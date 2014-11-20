$.widget("rbbt.secondary_structure_tool", {

  options: { },

  _create: function() {
    this.element.addClass('secondary_structure_tool_init');
    this.options.sequence = this.element.find('.sequence_tool');
    this.options.svg = this.element.find('.isoform_svg_tool');
  },

  _redirect_function: function(name){
    var args = arguments;
    this.options.sequence.sequence_tool(name, args[1], args[2], args[3]);
    this.options.svg.isoform_svg_tool(name, args[1], args[2], args[3]);
  },

  clear: function(){
    this._redirect_function('clear');
  },

  mark_position: function(position, color){
    if (typeof position == 'string') position = parseInt(position)
    this._redirect_function('mark_position', position, color);
  },

  mark_positions: function(positions, color){
    this._redirect_function('mark_positions', positions, color);
  },

  mark_region: function(start, end, color){
    this._redirect_function('mark_region', start, end, color);
  },
});
