require 'structure/alignment'

module Structure

  def self.sequence_position_in_pdb(protein_sequence, protein_positions, pdb, pdbfile)
    chains = PDBHelper.pdb_chain_sequences(pdb, pdbfile)

    alignments = {}
    chains.each do |chain,chain_sequence|

      chain_alignment, protein_alignment = SmithWaterman.align(chain_sequence, protein_sequence)

      alignments[chain] = Structure.match_position(protein_positions, protein_alignment, chain_alignment).compact
    end

    alignments.delete_if{|c,p| p.nil? or p.empty?}

    alignments
  end


  def self.pdb_chain_position_in_sequence(pdb, pdbfile, chain, positions, protein_sequence)
    chains = PDBHelper.pdb_chain_sequences(pdb, pdbfile)

    chain_sequence = chains[chain] 

    protein_alignment, chain_alignment = SmithWaterman.align(protein_sequence, chain_sequence)

    Structure.match_position(positions, chain_alignment, protein_alignment)
  end


  def self.pdb_alignment_map(protein_sequence, pdb, pdbfile)
    chains = PDBHelper.pdb_chain_sequences(pdb, pdbfile)

    result = TSV.setup({}, :key_field => "Sequence position", :fields => ["Chain:Position in PDB"], :type => :flat)
    chains.each do |chain,chain_sequence|
      chain_alignment, protein_alignment = SmithWaterman.align(chain_sequence, protein_sequence)

      map = Structure.alignment_map(protein_alignment, chain_alignment)

      map.each do |seq_pos, chain_pos|
        if result[seq_pos].nil?
          result[seq_pos] = [[chain, chain_pos] * ":"]
        else
          result[seq_pos] << [chain, chain_pos] * ":"
        end
      end
    end

    result
  end

  def self.pdb_positions_to_sequence(pdb_positions, sequence, target_chain, pdb = nil, pdbfile = nil)
    chain_positions = {}
    pdb_positions.collect do |cp|
      chain, position = cp.split(":")
      chain_positions[chain] ||= []
      chain_positions[chain] << position
    end
    return [] unless chain_positions.include? target_chain
    Structure.job(:pdb_chain_position_in_sequence, "TEST", :pdb => pdb, :pdbfile => pdbfile, :sequence => sequence, :chain => target_chain, :positions => chain_positions[target_chain]).exec
  end

  def self.neighbours_in_pdb(sequence, positions, pdb = nil, pdbfile = nil, chain = nil, distance = 5)

    positions_in_pdb = Structure.job(:sequence_position_in_pdb, "TEST", :pdb => pdb, :pdbfile => pdbfile, :sequence => sequence, :positions => positions).exec

    Log.debug "Position in PDB: #{Misc.fingerprint positions_in_pdb}"

    return [] if positions_in_pdb[chain].nil?

    chain ||=  positions_in_pdb.sort{|c,p| p.length}.last.first

    neighbour_map = Structure.job(:neighbour_map, "PDB Neighbours", :pdb => pdb, :pdbfile => pdbfile, :distance => distance).run

    return [] if neighbour_map.nil?

    inverse_neighbour_map = {}
    neighbour_map.each do |k,vs|
      vs.each do |v|
        inverse_neighbour_map[v] ||= []
        inverse_neighbour_map[v] << k
      end
    end

    neighbours_in_pdb = positions_in_pdb[chain].collect do |position|
      position_in_chain = [chain, position] * ":"
      neigh = neighbour_map[position_in_chain]
      ineigh = inverse_neighbour_map[position_in_chain]
      (neigh || []) + (ineigh || [])
    end.compact.flatten
  end

  def self.neighbours_in_pdb(sequence, positions, pdb = nil, pdbfile = nil, chain = nil, distance = 5)

    positions_in_pdb = Structure.job(:sequence_position_in_pdb, "TEST", :pdb => pdb, :pdbfile => pdbfile, :sequence => sequence, :positions => positions).exec

    Log.debug "Position in PDB: #{Misc.fingerprint positions_in_pdb}"

    chain ||=  positions_in_pdb.sort{|c,p| p.length}.last.first

    return {} if positions_in_pdb[chain].nil?

    neighbour_map = Structure.job(:neighbour_map, "PDB Neighbours", :pdb => pdb, :pdbfile => pdbfile, :distance => distance).run

    return {} if neighbour_map.nil?

    inverse_neighbour_map = {}
    neighbour_map.each do |k,vs|
      vs.each do |v|
        inverse_neighbour_map[v] ||= []
        inverse_neighbour_map[v] << k
      end
    end

    neighbours_in_pdb = TSV.setup({}, :key_field => "Original", :fields => ["Neighbours"], :type => :flat)
    positions_in_pdb[chain].each do |position|
      position_in_chain = [chain, position] * ":"
      neigh = neighbour_map[position_in_chain]
      ineigh = inverse_neighbour_map[position_in_chain]
      neighbours_in_pdb[position_in_chain] = (neigh || []) + (ineigh || [])
    end

    neighbours_in_pdb
  end

end

