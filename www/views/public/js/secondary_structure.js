$.widget("rbbt.secondary_structure_tool", {

  options: { },

  _create: function() {
    var tool = this;
    tool.element.addClass('secondary_structure_tool_init');
    tool.options.sequence = tool.element.find('.sequence_tool');
    tool.options.svg = tool.element.find('.isoform_svg_tool');
  },

  _redirect_function: function(name){
    var args = arguments;
    var tool = this;
    require_js(["/js/sequence.js", "/js/isoform_svg.js"], function(){
      tool.options.sequence.sequence_tool(name, args[1], args[2], args[3]);
      tool.options.svg.isoform_svg_tool(name, args[1], args[2], args[3]);
    })
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
