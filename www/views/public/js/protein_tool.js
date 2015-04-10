//$.widget("rbbt.protein_tool", {
//
//  options: { },
//
//  _create: function() {
//    this.element.addClass('protein_tool_init');
//
//    this.options.secondary_structure = this.element.find('.secondary_structure_tool');
//    this.options.jmol = this.element.find('.jmol_tool');
//  },
//
//  _redirect_function: function(name){
//    var args = arguments;
//
//    this.options.secondary_structure.secondary_structure_tool(name, args[1], args[2]);
//    this.options.secondary_structure.secondary_structure_tool(name, args[1], args[2], args[3], args[4]);
//
//    if (this.options.jmol !== undefined && this.options.jmol.jmol_tool('is_pdb_loaded')){ 
//      console.log(this.options.jmol.jmol_tool('loaded_pdb'))
//      this.options.jmol.jmol_tool(name, args[1], args[2], args[3], args[4])
//    }
//  },
//
//  clear: function(){
//    this._redirect_function('clear');
//  },
//
//  mark_position: function(position, color){
//    if (typeof position == 'string') position = parseInt(position)
//    this._redirect_function('mark_position', position, color);
//  },
//
//  mark_positions: function(positions, color){
//    this._redirect_function('mark_positions', positions, color);
//  },
//
//  mark_region: function(start, end, color){
//    this._redirect_function('mark_region', start, end, color);
//  },
//});


