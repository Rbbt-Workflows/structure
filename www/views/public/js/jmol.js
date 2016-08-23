var jmol_tools = [];


$.widget("rbbt.jmol_tool", {

  options: {
    width: '100%',
    height: 480,
    debug: false,
    color: "0xFFFFFF",
    addSelectionOptions: true,
    use: "HTML5",   // JAVA HTML5 WEBGL are all options
    j2sPath: "/js-find/jsmol/j2s", // this needs to point to where the j2s directory is.
    script: null,
    serverURL: '/get_pdb',
    readyFunction: null,
    disableJ2SLoadMonitor: true,
    disableInitialConsole: true,
    allowJavaScript: true
  },

  _create: function() {
    var tool = this;

    require_js(["/js-find/jsmol/JSmol.min.nojq.js"], function(){
      Jmol._z=0
      Jmol._serverUrl='/get_pdb'
      tool.element.addClass('jmol_tool_init')
      tool.options.jmol_window = tool.element.find('.window');
      tool.options.jmol_window.html(Jmol.getAppletHtml("jmolApplet0", tool.options))
      tool.options.applet = jmolApplet0
    })

  },

  _wrapper: function() {
    var tool = this
    return {script: function(code){Jmol.script(tool.options.applet, code)}, getProperty: function(prop){return JSON.parse(Jmol.getPropertyAsJSON(tool.options.applet, prop))}};
  },

  background: function(color){
    this._wrapper().script('background ' + color);
  },

  clear: function(){
    var resetStyleScript = "select all; wireframe off; spacefill off; cartoon off; ribbon off; rocket off; strand off; trace off; halos off;select protein; backbone off; color pink;cartoons on;color structure; ";
    this._wrapper().script(resetStyleScript);
  },

  // PDB LOADING

  load_pdb: function(pdb) {
    this._loaded_pdb = pdb.replace(/^=/,'')
    this._wrapper().script("load " + pdb + "; wireframe off; restrict helix and sheet; select protein; backbone off; color pink;cartoons on;color structure;");
  },

  is_pdb_loaded: function(){
    var wrapper = this._wrapper();
    var filename = wrapper.getProperty("filename").filename;
    return (undefined != filename && filename != 'zapped' && filename != []);
  },

  pdbfile: function(){
    return this._wrapper().getProperty("filename");
  },

  //{{{ JOBS

  _mutation_positions: function(mutations, organism, watson, complete){
    var mutated_isoforms = rbbt_job("Sequence", "mutated_isoforms_fast", {organism: organism, watson: watson, mutations: mutations}, complete)
    return mutated_isoforms;
  },

  _sequence_positions_in_pdb: function(positions, complete){
    var job = new rbbt.Job("Structure", "sequence_position_in_pdb", {sequence: this.options.sequence, pdb: this._loaded_pdb, positions: positions.join("|")})
    job.run(true).then(complete)
  },

  alignment_map: function(complete){
    var job = new rbbt.Job("Structure", "pdb_alignment_map", {sequence: this.options.sequence, pdb: this._loaded_pdb})
    job.run(true).then(complete)
  },

  pdb_alignment_map: function(){
    var job = new rbbt.Job("Structure", "pdb_alignment_map", {sequence: this.options.sequence, pdb: this._loaded_pdb})
    var promise = job.run(true)
    return promise
  },

  //{{{ JMOL STUFF

  _select: function(position, chain){
    if (position instanceof Array){
      var tmp = $(position).map(function(){ return (parseInt(this)) }).toArray()
      position_str = tmp.join(", ");
    }else{
        position_str = position
    }

    this._wrapper().script("select protein and *.CA and " + position_str + ":" + chain + ";");
  },

  _color: function(color){
    this._wrapper().script("color " + color + ";");
  },

  _halos: function(color){
    this._wrapper().script("color halos " + color + ";");
    this._wrapper().script("halos on;");
  },

  _style: function(style){
    this._wrapper().script(style + ";");
  },

  //{{{ MANIPULATION

  mark_chain_position: function(chain, position, color){
    this._select(position, chain);
    this._style("spacefill")
    if (undefined === color){
      this._halos("color");
    }else{
      this._color(color)
    }
  },

  mark_chain_region: function(chain, start, end, color){
    if (start == end){
      this._select(start - 1, chain);
    }else{
      this._select((start - 1) + '-' + (end - 1), chain);
    }

    this._style("cartoon")
    if (undefined === color){
      this._halos("color");
    }else{
      this._color(color)
    }
  },

  //{{{ HIGH LEVEL

  mark_positions: function(positions, color){
    var tool = this;
    tool._sequence_positions_in_pdb(positions, function(pdb_positions){
     console.log(pdb_positions)
      for (var chain in pdb_positions){
        var position_list = pdb_positions[chain]
        position_list = $.grep(position_list,function(n){ return(n)});
        if (position_list != null && position_list.length > 0){
          tool.mark_chain_position(chain, position_list, color);
        }
      }
    })
  },

  mark_position: function(position, color){
    this.mark_positions([position], color)
  },


  mark_region: function(start, end, color){
    var tool = this;
    tool._sequence_positions_in_pdb([start, end], function(pdb_positions){
      for (var chain in pdb_positions){
        var positions = pdb_positions[chain]
        if (null != positions[0] && null != positions[1]){
          tool.mark_chain_region(chain, positions[0], positions[1], color)
        }
      }
    })
  },

  mark_genomic_mutations: function(list, color){
    var protein = this.options.protein;

    var info = list_info("GenomicMutation", list);
    var organism = info.organism
    var watson = info.watson
    var tool = this

    this._mutation_positions(list_array("GenomicMutation", list), organism, watson, function(mutated_isoforms){
      var all_positions = [];
      for (var mutation in mutated_isoforms){
        var mis = mutated_isoforms[mutation];

        var good_mis = $.grep(mis, function(mi){ return(null !== mi.match(protein))});

        var changes = $(good_mis).map(function(){
          return this.split(/:/)[1];
        })

        var good_changes = $.grep(changes, function(c){ var m = c.match(/([A-Z]*)\d+([A-Za-z*]*)/); return(m != null && m[1] != m[2])});

        var positions = $(good_changes).map(function(){
          return this.match(/[A-Z](\d+)[A-Za-z*]/)[1];
        }).toArray()

        all_positions = all_positions.concat(positions)
      }

      tool.mark_positions(all_positions, color)
    });
  },

  mark_region: function(color, map){
    var tool = this
    var chains = {};
    for(seq_pos in map){
      var chain_pos_list = map[seq_pos]
      for (var i = 0; i < chain_pos_list.length; i++) {
        var chain_pos = chain_pos_list[i];
        var chain = chain_pos.split(":")[0]
        var pos = chain_pos.split(":")[1]

        if (undefined === chains[chain]){
          chains[chain] = [];
        }
        chains[chain].push(parseInt(pos));
      }
    }
    for (chain in chains){
      var positions = chains[chain];
      positions = $(positions).map(function(){return parseInt(this);}).toArray().sort(function(a,b){return a-b});
      var last = -1
      var start = -1;
      for (var i = 0; i < positions.length; i++){
        if (positions[i] != last + 1){
          if (start != -1) { tool.mark_chain_region(chain, start, last, color); }
          start = positions[i]
        }
        last = positions[i]
      }
      if (start != -1) { tool.mark_chain_region(chain, start, last, color); }
    }
  },

  mark_aligned_region: function(color){
    var tool = this
    this.alignment_map(function(map){
      tool.mark_region(color, map)
    })
  },

});
