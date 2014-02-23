require 'structure/ssw'

module Structure

  def self.alignment_map(alignment_source, alignment_target)
    map = {}

    offset_source, alignment_source = alignment_source.match(/^(_*)(.*)/).values_at( 1, 2)
    offset_target, alignment_target = alignment_target.match(/^(_*)(.*)/).values_at( 1, 2)
 
    gaps_source = 0 
    gaps_target = 0
    alignment_source.chars.zip(alignment_target.chars).each_with_index do |p,i|
      char_source, char_target = p
      gaps_source += 1 if char_source == '-'
      gaps_target += 1 if char_target == '-'
      map[i + 1 + offset_source.length - gaps_source] = i + 1 + offset_target.length - gaps_target if char_source == char_source  and char_source != "-"
    end
    
    map
  end

  def self.sequence_map(source_sequence, target_sequence)
    source_alignment, target_alignment = SmithWaterman.align(source_sequence, target_sequence)
    Structure.alignment_map(source_alignment, target_alignment)
  end

  def self.match_position(protein_position, protein_alignment, chain_alignment)
    map = alignment_map(protein_alignment, chain_alignment)
    if Array === protein_position
      map.chunked_values_at(protein_position.collect{|p| p.to_i})
    else
      map[protein_position]
    end
  end
end
