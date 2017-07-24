/**
 * Viewer: a NGL wrapper
 */
require_js(["/js-find/ngl/ngl.min.js"], function() {
  $.widget('rbbt.ngl_viewer', {
    options: {
      pdb_url: null,
      pdb_url_download: null,
      pdb: null,
      component_promise: null,
      stage: null,
      colorScheme: [],
      seq2pdb: null,
      pdb2seq: null,
      default_colors: null,
      marks: {},
      layers_sequence: []
    },

    _create: function() {
      this.options.stage = new NGL.Stage("ngl_viewport", {backgroundColor: 'black'});
      var tooltip = document.createElement("div");
      Object.assign(tooltip.style, {
        display: "none",
        position: "absolute",
        zIndex: 10,
        pointerEvents: "none",
        backgroundColor: "rgba(0, 0, 0, 0.6)",
        color: "white",
        padding: "0.5em",
        fontFamily: "sans-serif"
      });

      this._init();

      // For use this inside mouse hooks
      var ths = this;

      this.options.stage.viewer.container.appendChild(tooltip);
      // remove default hoverPick mouse action
      this.options.stage.mouseControls.remove("hoverPick");

      // listen to `hovered` signal to move tooltip around and change its text
      this.options.stage.signals.hovered.add(function(pickingProxy){
        if(pickingProxy && (pickingProxy.atom || pickingProxy.bond)){
          var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
          var cp = pickingProxy.canvasPosition;
          var group = atom.qualifiedName().match(/\[(\w+)\](\d+)(.+)/);
          var aa = group[1];
          var pos_pdb = group[2];
          var first_pos_chain = Object.keys(ths.options.pdb2seq)[0];
          var chain_letters = new Set(Object.keys(ths.options.pdb2seq).map(pdb_pos => pdb_pos[0]));
          for (let chain_letter of chain_letters) {
            var pos_seq = ths.options.pdb2seq[chain_letter+':'+ pos_pdb];
            if (pos_seq) break;
          }
          var rest = group[3];
          tooltip.innerText = "AA: " + "[" + aa + "] " + (pos_seq === undefined ? ' not alignment here! ' : pos_seq + ' seq / '+ pos_pdb + ' pdb') + rest;
          tooltip.style.bottom = (cp.y + 18) + "px";
          tooltip.style.left = (cp.x + 15)  + "px";
          tooltip.style.display = "block";
        }else{
          tooltip.style.display = "none";
        }
      });

      this.options.stage.mouseControls.remove("clickPick");
      this.options.stage.signals.clicked.add(function(pickingProxy){
        if(pickingProxy && (pickingProxy.atom || pickingProxy.bond)){
          var buttons = ths.options.stage.mouseControls.mouse.buttons;
          // right click
          if (buttons == 2) {
            // implement here the context menu when right button is pressed
          } else {
          var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
          var cp = pickingProxy.canvasPosition;
          var group = atom.qualifiedName().match(/\[(\w+)\](\d+)(.+)/);
          var aa = group[1];
          var pos_pdb = group[2];
          var first_pos_chain = Object.keys(ths.options.pdb2seq)[0];
          var chain_letters = new Set(Object.keys(ths.options.pdb2seq).map(pdb_pos => pdb_pos[0]));
          for (let chain_letter of chain_letters) {
            var pos_seq = ths.options.pdb2seq[chain_letter+':'+ pos_pdb];
            if (pos_seq) break;
          }
          var rest = group[3]
          tooltip.innerText = "AA: " + "[" + aa + "] " + (pos_seq === undefined ? ' not alignment here! ' : pos_seq + ' seq / '+ pos_pdb + ' pdb') + rest;
          tooltip.style.bottom = (cp.y + 18) + "px";
          tooltip.style.left = (cp.x + 15)  + "px";
          tooltip.style.display = "block";
          if (pos_seq === undefined) return
          ths.mark_positions_by_colors('click', [pos_pdb], 'pdb', 'red')
          }
        }
      });
    },

    /**
     * Viewer constructor
     */
    _init: function() {
      this.options.component_promise = new Promise((resolve, reject) => {
        if (this.options.pdb_url.match(/interactome3d/)) {
          this.options.pdb_url_download ='/interactome3d?=pdb'+this.options.pdb_url;
        } else {
          this.options.pdb_url_download = this.options.pdb_url;
        }
        this.options.stage.removeAllComponents();
        this._reset();
        this._update_color_scheme();

        var ths = this;
        // only first model in NMR
        this.options.stage.loadFile(this.options.pdb_url_download, {ext: 'pdb', firstModelOnly: true}).then(function(component){
          component.addRepresentation("cartoon", {
            color: ths.options.colorSchemeNGL,
            sele: "protein and .CA",
            assembly: 'AU',
            multipleBond: true
          });
          component.autoView();
          resolve(component);
        });
      });
      this.options.stage.setSpin(false);
      this.options.stage.rockAnimation.pause(true)
    },

    /**
     * For drawing colors in order
     */
    _update_layers_sequence: function(layer) {
      var i = this.options.layers_sequence.indexOf(layer);
      if(i != -1) this.options.layers_sequence.splice(i, 1);
      if (layer != this.options.layers_sequence[0]) this.options.layers_sequence.unshift(layer);
    },

    /**
     * Generate the color scheme to be drawn by the viewer
     */
    _update_color_scheme: function(){
      this.options.colorScheme = [];

      for (let layer of this.options.layers_sequence) {
        for(let color of Object.keys(this.options.marks[layer])) {
          if (this.options.marks[layer][color].size == 0) continue;
          this.options.colorScheme.push([color, Array.from(this.options.marks[layer][color]).join(' or ')]);
        }
      }

      this.options.colorScheme.push([this.options.default_colors['default'], '*']);
      this.options.colorSchemeNGL = NGL.ColormakerRegistry.addSelectionScheme(this.options.colorScheme, "colorScheme");
    },

    draw: function(){
      this._update_color_scheme();
      var ths = this;
      this.options.component_promise.then(component => {component.addRepresentation("cartoon", {
        color: ths.options.colorSchemeNGL,
        sele: "protein and .CA",
        assembly: 'AU',
      });
        //component.autoView();
      });
    },

    _reset: function(){
      this.options.marks = {};
      this.options.marks['original_alignment'] = {}
      this.options.marks['original_alignment'][this.options.default_colors['original_alignment']] = new Set();
      this.options.layers_sequence = [];
      for (let pos of Object.keys(this.options.seq2pdb)) {
        this.options.marks['original_alignment'][this.options.default_colors['original_alignment']].add(this._get_pos(pos, 'seq'));
      }
      this._update_layers_sequence('original_alignment');
    },

    clear_all: function(){
      this._reset();
      this.draw();
    },

    clear_by_layers: function(layers){
      for (let layer of layers) {
        this.options.marks[layer] = new Set();
      }
      this.draw();
    },

    clear_by_layer: function(layer){
      this.clear_by_layers([layer])
    },

    /**
     * Get the PDB position from a seq/pdb position
     */
    _get_pos: function(pos, seq_or_pdb) {
      if (seq_or_pdb == 'seq') {
        chains = this.options.seq2pdb[pos];
        pos = chains[0].match(/[A-Z]:(\d+)/)[1];
      }
      return pos;
    },

    _mark_position: function(layer, pos, seq_or_pdb, color) {
      this.options.marks[layer][color].add(this._get_pos(pos, seq_or_pdb));
    },

    /*
     * Default colors: mark, align, original_alignment, feature, color
     */
    mark_positions_by_colors: function(layer, pos_list, seq_or_pdb, color_or_list, draw = true) {
      this._update_layers_sequence(layer);
      layer in this.options.marks || (this.options.marks[layer] = new Set());;
      if (Array.isArray(color_or_list)) {
        for (i = 0; i < pos_list.length; i++){
          var color = color_or_list[i];
          color in this.options.marks[layer] || (this.options.marks[layer][color] = new Set());;
          this._mark_position(layer, pos_list[i], seq_or_pdb, color);
        }
      } else {
        var color = color_or_list
        color in this.options.marks[layer] || (this.options.marks[layer][color] = new Set());;
        for (let pos of pos_list) {
          this._mark_position(layer, pos, seq_or_pdb, color);
        }
      }
      if (draw) this.draw();
    },

    color_mutation_density_subset: function(residues, residue_incidence) {
      var pos_list = [];
      var log10_counts = [];
      for (let pos of residues) {
        if (this.options.seq2pdb[pos]){
          var count = parseInt(residue_incidence[pos]) + 1;
          var log_count = Math.log(count) / Math.log(10);
          log10_counts.push(log_count);
          pos_list.push(pos);
        }
      }

      log10_counts.push(0);

      var colors = get_gradient(log10_counts, '#00FF00', '#FF0000');

      return {'pos_list': pos_list, 'colors': colors };
    },

    color_mutation_density: function(residue_incidence) {
      pos_colors = this.color_mutation_density_subset(Object.keys(residue_incidence), residue_incidence);
      this.mark_positions_by_colors('cosmic_mutation_density', pos_colors['pos_list'], 'seq', pos_colors['colors']);
    },

    screenshot: function() {
      this.options.stage.makeImage({factor: 1, antialias: true}).then(blob => NGL.download(blob, this.options.pdb + ".png"));
    },

    resize: function() {
      this.options.stage.autoView();
    },

    spin: function() {
      this.options.stage.toggleSpin();
    },

    fullscreen: function() {
      this.options.stage.toggleFullscreen();
    }
  })
})
