rbbt.sequence = {}
rbbt.sequence.clear = function(element){
  var seq = element.find('span.sequence');
  var marked_chars = seq.find('.sequence_char_position')
  marked_chars.each(function(marked_char){ marked_char = $(this); marked_char.replaceWith(marked_char.html()) })
}

rbbt.sequence.mark_position = function(element, position, color){
  if (undefined === color) color = 'red';
  var seq = element.find('span.sequence');
  var str = seq.html();
  var char_pos = position - 1;
  str = str.slice(0, char_pos) + $('<span class="sequence_char_position" data-sequence_position=' + position + ' style="color:' +  color + '">' + str[char_pos] + '</span>')[0].outerHTML + str.slice(char_pos+1)
  seq.html(str)
}

rbbt.svg = {}

rbbt.svg.clear = function(element, layer=''){
  var svg = element.find('svg');
  var vlines;
  var vregions;
  if (layer) {
    vlines = svg.find('.rbbt-vline.' + layer)
    vregions = svg.find('.rbbt-region.' + layer)
  } else {
    vlines = svg.find('.rbbt-vline')
    vregions = svg.find('.rbbt-region')
  }
  vlines.remove()
  vregions.remove()
}

rbbt.svg.position_offset = function(element, position){
  var svg = element.find('svg');
  var width  = parseInt(svg.attr('width'));
  var start  = parseInt(svg.find('rect.ac').attr('x'));
  var seq_len  = parseInt(element.attr('data-sequence_length'));
  return start + position * (width - start)/seq_len
}

rbbt.svg.mark_position = function(element, position, color, layer=''){
  if (undefined === color) color = 'red';
  var svg = element.find('svg');

  var width  = parseInt(svg.attr('width'));
  var height  = parseInt(svg.attr('height'));
  var offset = rbbt.svg.position_offset(element, position)

  var line =  document.createElementNS("http://www.w3.org/2000/svg", "line");
  line.setAttributeNS(null, "x1", offset);
  line.setAttributeNS(null, "y1", 5);
  line.setAttributeNS(null, "x2", offset);
  line.setAttributeNS(null, "y2", height - 5);
  line.setAttributeNS(null, "class", 'rbbt-vline' + (layer ? ' ' + layer : ''));
  line.setAttributeNS(null, "style", "stroke:" + color + ";opacity:0.5;stroke-width:1;");

  svg.append(line);
}

rbbt.svg.mark_positions = function(element, positions, color, layer=''){
  for (var i = 0; i < positions.length; i++){
    var position = positions[i];
    var seq_len  = parseInt(element.attr('data-sequence_length'));
    jitter = Math.ceil(seq_len / 1000) + 2;
    position = position + Math.floor(Math.random() * 2 * jitter) -jitter + 1
    if (Array.isArray(color)) {
      rbbt.svg.mark_position(element, position, color[i], layer);
    } else {
      rbbt.svg.mark_position(element, position, color, layer);
    }
  }
}

rbbt.svg.mark_region = function(element, first, last, color, layer=''){
  var svg = element.find('svg')
  var width = parseInt(svg.attr('width'));
  var height  = parseInt(svg.attr('height'));
  var start = parseInt(svg.find('rect.ac').attr('x'));
  var seq_len  = parseInt(element.attr('data-sequence_length'));

  var offset_start = rbbt.svg.position_offset(element, first);
  var offset_end = rbbt.svg.position_offset(element, last);

  var rect =  document.createElementNS("http://www.w3.org/2000/svg", "rect");

  rect.setAttributeNS(null, "x", offset_start);
  rect.setAttributeNS(null, "class", 'rbbt-region' + (layer ? ' ' + layer : ''));
  rect.setAttributeNS(null, "y", 10);
  rect.setAttributeNS(null, "width",  offset_end - offset_start);
  rect.setAttributeNS(null, "height", height - 30);
  rect.setAttributeNS(null, "style", "stroke:black;stroke-width: 2; opacity:0.3;fill:" + color + ";");

  svg.append(rect);
}

rbbt.svg.mark_aligned_region = function(element, map, color, layer=''){
  if (undefined === color) color = 'blue'

    var positions = []
    for(seq_pos in map){
      positions.push(parseInt(seq_pos));
    }

    positions = positions.sort();

    var last = -1;
    var start = -1;
    for (var i = 0; i < positions.length; i++){
      if (positions[i] != last + 1){
        if (start != -1) { rbbt.svg.mark_region(element, start, last, color, layer); }
        start = positions[i];
      }
      last = positions[i];
    }
    if (start != -1) { rbbt.svg.mark_region(element, start, last, color, layer); }
}

