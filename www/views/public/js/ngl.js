require_js(["/js-find/ngl/ngl.min.js"], function() {y:
  $.widget('rbbt.ngl_tool', {
    options: {
      pdb_url: null,
      pdb_url_download: null,
      pdb: url,
      component_promise: null,
      stage: null,
      colorScheme: null,
      seq2pdb: null,
      pdb2seq: null,
      appris_features: null,
      marks: {},
      default_colors: {'mark': 'red', 'original_alignment': 'darkslategray', 'align': 'blue', 'feature': 'purple', 'default': 'lightgray'}
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

      this._init()

      // For use this inside mouse hooks
      var ths = this

      this.options.stage.viewer.container.appendChild(tooltip);
      // remove default hoverPick mouse action
      this.options.stage.mouseControls.remove("hoverPick");

      // listen to `hovered` signal to move tooltip around and change its text
      this.options.stage.signals.hovered.add(function(pickingProxy){
        if(pickingProxy && (pickingProxy.atom || pickingProxy.bond)){
          var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
          var cp = pickingProxy.canvasPosition;
          group = atom.qualifiedName().match(/\[(\w+)\](\d+)(.+)/)
          aa = group[1]
          pos_pdb = group[2]
          first_pos_chain = Object.keys(ths.options.pdb2seq)[0]
          chain_letters = new Set(Object.keys(ths.options.pdb2seq).map((pdb_pos) => pdb_pos[0]))
          for (let chain_letter of chain_letters) {
            pos_seq = ths.options.pdb2seq[chain_letter+':'+ pos_pdb]
            if (pos_seq) break
          }
          rest = group[3]
          tooltip.innerText = "AA: " + "[" + aa + "]" + (pos_seq === undefined ? ' not alignment here! ' : pos_seq + 'seq/'+ pos_pdb + 'pdb') + rest
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
          var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
          var cp = pickingProxy.canvasPosition;
          group = atom.qualifiedName().match(/\[(\w+)\](\d+)(.+)/)
          aa = group[1]
          pos_pdb = group[2]
          first_pos_chain = Object.keys(ths.options.pdb2seq)[0]
          chain_letters = new Set(Object.keys(ths.options.pdb2seq).map((pdb_pos) => pdb_pos[0]))
          for (let chain_letter of chain_letters) {
            pos_seq = ths.options.pdb2seq[chain_letter+':'+ pos_pdb]
            if (pos_seq) break
          }
          rest = group[3]
          tooltip.innerText = "AA: " + "[" + aa + "]" + (pos_seq === undefined ? ' not alignment here! ' : pos_seq + 'seq/'+ pos_pdb + 'pdb') + rest
          tooltip.style.bottom = (cp.y + 18) + "px";
          tooltip.style.left = (cp.x + 15)  + "px";
          tooltip.style.display = "block";
          ths.mark_position(pos_pdb, 'pdb')
        }
      });
    },

    _init: function() {
      this.options.component_promise = new Promise((resolve, reject) => {
        if (this.options.pdb_url.match(/interactome3d/)) {
          this.options.pdb_url_download ='/interactome3d?=pdb'+this.options.pdb_url
        } else {
          this.options.pdb_url_download = this.options.pdb_url
        }
        this.options.stage.removeAllComponents()
        this._reset()

        // only first model in NMR
        this.options.stage.loadFile(this.options.pdb_url_download, {ext: 'pdb', firstModelOnly: true}).then(function(component){
          //if(component.type !== "structure") return;
          component.addRepresentation("cartoon", {
            //color: "sstruc" ,
            sele: "protein and .CA",
            assembly: 'AU',
            multipleBond: true
          });
          component.autoView();
          resolve(component)
        })
      })
      this.clear()
    },

    _draw: function(){
      this.options.colorScheme_promise = new Promise((resolve, reject) => {
        var colorScheme = []
        if (this.options.marks['mark'].size > 0) {
          colorScheme.push([this.options.default_colors['mark'], Array.from(this.options.marks['mark']).join(' or ')])
        }
        if (this.options.marks['feature'].size > 0) {
          colorScheme.push([this.options.default_colors['feature'], Array.from(this.options.marks['feature']).join(' or ')])
        }
        if (Object.keys(this.options.marks).length > 4) {
          for (let key of Object.keys(this.options.marks)){
            if (!Object.keys(this.options.default_colors).includes(key)) {
              colorScheme.push([key, Array.from(this.options.marks[key]).join(' or ')])
            }
          }
        }
        if (this.options.marks['align'].size > 0) {
          colorScheme.push([this.options.default_colors['align'], Array.from(this.options.marks['align']).join(' or ')])
        }
        if (this.options.marks['original_alignment'].size > 0) {
          colorScheme.push([this.options.default_colors['original_alignment'], Array.from(this.options.marks['original_alignment']).join(' or ')])
        }
        colorScheme.push([this.options.default_colors['default'], '*'])
        resolve(colorScheme)
      })
      this.options.colorScheme_promise.then((colorScheme) => {
        this.options.colorScheme = NGL.ColormakerRegistry.addSelectionScheme(colorScheme, "colorScheme" );
        this.options.component_promise.then((component) => {component.addRepresentation("cartoon", {
          color: this.options.colorScheme,
          sele: "protein and .CA",
          assembly: 'AU',
        })
          component.autoView()
        })
      })
    },

    _reset: function(){
      this._reset_marks()
    },

    clear: function(){
      this._reset()
      this._draw()
    },

    _get_pos: function(pos, seq_or_pdb) {
      if (seq_or_pdb == 'seq') {
        chains = this.options.seq2pdb[pos]
        pos = chains[0].match(/[A-Z]:(\d+)/)[1]
      }
      return pos
    },

    _reset_marks: function() {
      this.options.marks = {}
      this.options.marks['mark'] = new Set()
      this.options.marks['feature'] = new Set()
      this.options.marks['align'] = new Set()
      this.options.marks['original_alignment'] = new Set()
      for (let pos of Object.keys(this.options.seq2pdb)) {
        this.options.marks['original_alignment'].add(this._get_pos(pos, 'seq'))
      }
    },

    // mark, align, original_alignment, feature, color
    _mark_position: function(pos, seq_or_pdb, color) {
      color in this.options.marks || (this.options.marks[color] = new Set());
      this.options.marks[color].add(this._get_pos(pos, seq_or_pdb))
    },

    mark_positions_by_colors: function(pos_list, seq_or_pdb, color_or_list) {
      if (Array.isArray(color_or_list)) {
        for (i = 0; i < pos_list.length; i++){
          this._mark_position(pos_list[i], 'seq', color_or_list[i])
        }
      } else {
        for (let pos of pos_list) {
          this._mark_position(pos, 'seq', color_or_list)
        }
      }
      this._draw()
    },

    mark_position: function(pos, seq_or_pdb){
      this._mark_position(pos, seq_or_pdb, 'mark')
      this._draw()
    },

    mark_positions: function(pos_list, seq_or_pdb) {
      this.mark_positions_by_colors(pos_list, 'seq', 'mark')
    },

    mark_appris_features: function() {
      pos_seq_features = []
      for (let feature of this.options.appris_features) {
        pos_seq = feature['start']
        pos_pdb = this.options.seq2pdb[pos_seq]
        if (feature['type'] == 'firestar' && pos_pdb) {
          pos_seq_features.push(pos_seq)
        }
      }
      this.mark_positions_by_colors(pos_seq_features, 'seq', 'feature')
    },

    align: function(){
      this._reset()
      this.mark_positions_by_colors(Object.keys(this.options.seq2pdb), 'seq', 'align')
    },

    color_mutation_density_subset: function(residues, residue_incidence) {
      pos_list = []
      var log10_counts = [];
      for (let pos of residues) {
        if (this.options.seq2pdb[pos]){
          var count = parseInt(residue_incidence[pos]) + 1;
          var log_count = Math.log(count) / Math.log(10);
          log10_counts.push(log_count)
          pos_list.push(pos)
        }
      }

      log10_counts.push(0)

      colors = get_gradient(log10_counts, '#00FF00', '#FF0000');

      return {'pos_list': pos_list, 'colors': colors }
    },

    color_mutation_density: function(residue_incidence) {
      console.log('color_mutation_density')
      this._reset()
      pos_colors = this.color_mutation_density_subset(Object.keys(residue_incidence), residue_incidence)
      this.mark_positions_by_colors(pos_colors['pos_list'], 'seq', pos_colors['colors'])
    },

    screenshot: function() {
      this.options.stage.makeImage({factor: 1, antialias: true}).then((blob) => NGL.download(blob, this.options.pdb + ".png"))
    },

    resize: function() {
      this.options.stage.autoView()
    },

    spin: function() {
      this.options.stage.toggleSpin()
    },

    fullscreen: function() {
      this.options.stage.toggleFullscreen()
    }
  })
})
