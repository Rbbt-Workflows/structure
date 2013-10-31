$.widget("rbbt.sequence_tool", {

  options: {},

  _create: function() {
    this.options.seq_len = this.options.sequence.length;
  },

  _clear: function(){
    var seq = this.element.find('span.sequence');
    seq.html(this.options.sequence)
  },

  _mark_position: function(pos, color){
    pos = pos - 1
    var seq = this.element.find('span.sequence');
    var str = seq.html();
    str = str.slice(0, pos) + $('<span class="sequence_position" style="color:' +  color + '">' + str[pos] + '</span>')[0].outerHTML + str.slice(pos+1)
    seq.html(str)
  },

  _mark_region: function(start, end, color){
    start = start - 1
    end = end - 1
    var seq = this.element.find('span.sequence');
    var str = seq.html();
    str = str.slice(0, start) + $('<span class="sequence_range" style="color:' +  color + '">' + str.slice(start, end+1) + '</span>')[0].outerHTML + str.slice(end+1)
    seq.html(str)
  },

  //{{{ JOBS

  _mutation_positions: function(mutations, organism, watson, complete){
    var mutated_isoforms = rbbt_job("Sequence", "mutated_isoforms_for_genomic_mutations", {organism: organism, watson: watson, mutations: mutations}, complete);
    return mutated_isoforms;
  },

  // HIGH LEVEL

  place: function(left, size){
    this.element.css('margin-left', left).css('width', size)
  },

  mark_position: function(position, color){
    this._clear();
    this._mark_position(position, color);
  },

  mark_positions: function(positions, color){
    this._clear();
    positions = _.uniq(positions.sort(function(a,b){return a-b}).reverse());
    for(i in positions){
      var position = positions[i]
      this._mark_position(position, color);
    }
  },

  mark_region: function(start, end, color){
    this._clear();
    this._mark_region(start, end, color);
  },
 
  mark_genomic_mutations: function(list, color, jitter){
    var info = list_info("GenomicMutation", list);
    var organism = info.organism
    var watson = info.watson

    var tool = this;
    this._mutation_positions(list_array("GenomicMutation", list), organism, watson, function(mutated_isoforms){

      var protein = tool.options.protein;
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

      if (all_positions.length > 0){
        tool.mark_positions(all_positions, color, jitter)
      }
    });
  },

  clear: function(){
    this._clear()
  }
})

