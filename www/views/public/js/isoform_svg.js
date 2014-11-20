$.widget("rbbt.isoform_svg_tool", {

  options: {},

  _create: function() {
    this.options.seq_len = this.options.sequence.length;
  },

  _svg: function(){
    return $(this.element.find('.window > svg')[0]);
  },

  clear: function(){
    this._svg().find('.rbbt-vline, .rbbt-region').remove();
  },

  //{{{ JOBS

  _mutation_positions: function(mutations, organism, watson, complete){
    var mutated_isoforms = rbbt_job("Sequence", "mutated_isoforms_fast", {organism: organism, watson: watson, mutations: mutations}, complete);
    return mutated_isoforms;
  },

  //{{{ SVG STUFF

  _position_in_svg: function(position, jitter){
    var seq_len = this.options.seq_len;
    var svg = this._svg();
    var width  = parseInt(svg.attr('width'));
    var start  = parseInt(svg.find('rect.ac').attr('x'));

    if (position.toArray){
      position = position.toArray()
    }

    if (position instanceof Array){ 
      var tool = this
      return $(position).map(function(){
        return tool._position_in_svg(this, jitter);
      }).toArray();
    }else{
      offset = ((width - start)/seq_len * parseInt(position)) + start
      if (jitter){ offset += Math.random() * 10 - 5; }
      return offset;
    }
  },

  _vline: function(x, color, options){
    var line =  document.createElementNS("http://www.w3.org/2000/svg", "line");
    line.setAttributeNS(null, "x1", x);
    line.setAttributeNS(null, "y1", 5);
    line.setAttributeNS(null, "x2", x);
    line.setAttributeNS(null, "y2", this._svg().attr('height') - 5);
    line.setAttributeNS(null, "class", 'rbbt-vline');
    line.setAttributeNS(null, "style", "stroke:" + color + ";opacity:0.5;stroke-width:1;");

    if (options !== null){
      for (option in options){
        line.setAttributeNS(null, option, options[option]);
      }
    }

    this.element.find('.window > svg')[0].appendChild(line);
  },

  _region: function(start, end, color, options){
    var rect =  document.createElementNS("http://www.w3.org/2000/svg", "rect");

    rect.setAttributeNS(null, "x", start);
    rect.setAttributeNS(null, "class", 'rbbt-region');
    rect.setAttributeNS(null, "y", 10);
    rect.setAttributeNS(null, "width",  end - start);
    rect.setAttributeNS(null, "height", this._svg().attr('height') - 30);
    rect.setAttributeNS(null, "style", "stroke:black;stroke-width: 2; opacity:0.3;fill:" + color + ";");

    if (options !== null){
      for (option in options){
        rect.setAttributeNS(null, option, options[option]);
      }
    }

    this.element.find('.window > svg')[0].appendChild(rect);

  },

  //{{{ MANIPULATION
  
  mark_position: function(position, color, jitter){
    var tool = this;
    var svg = tool._svg()
    var seq_len = tool.options.seq_len;
    var width  = parseInt(svg.attr('width'));
    var start  = parseInt(svg.find('rect.ac').attr('x'));

    $(this._position_in_svg($(position), jitter)).each(function(){
      tool._vline(this, color) 
    })
  },

  mark_positions: function(positions, color, jitter){
    var tool = this;
    $(positions).each(function(){ tool.mark_position(this, color, jitter); })
  },

  mark_region: function(r_start, r_end, color){
    var svg = this.element.find('.window > svg')
    var width  = parseInt(svg.attr('width'));
    var start  = parseInt(svg.find('rect.ac').attr('x'));
    var seq_len = this.options.seq_len;

    var offset_start = this._position_in_svg(r_start)
    var offset_end = this._position_in_svg(r_end)

    this._region(offset_start, offset_end, color)
  },


  //{{{ HIGH LEVEL
 
  mark_genomic_mutations: function(list, color, jitter){
    var info = list_info("GenomicMutation", list);
    var organism = info.organism
    var watson = info.watson

    var tool = this;
    this._mutation_positions(list_array("GenomicMutation", list), organism, watson, function(mutated_isoforms){

      var protein = tool.options.protein;
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

        if (positions.length > 0){
          tool.mark_position(positions, color, jitter)
        }
        }
      });
  },

  mark_aligned_region: function(map, color){
    var positions = []
    for(seq_pos in map){
      positions.push(parseInt(seq_pos))
    }
    positions = positions.sort();
    var last = -1
    var start = -1;
    for (var i = 0; i < positions.length; i++){
      if (positions[i] != last + 1){
        if (start != -1) { this.mark_region(start, last, color); }
        start = positions[i]
      }
      last = positions[i]
    }
    if (start != -1) { this.mark_region(start, last, color); }
  }

})

