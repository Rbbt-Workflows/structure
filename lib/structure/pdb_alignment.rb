require 'structure/alignment'

module Structure

  SSW_ALIGNMENT_REPO = cache_dir.ssw_alignments.find
  Open.repository_dirs << SSW_ALIGNMENT_REPO unless Open.repository_dirs.include? SSW_ALIGNMENT_REPO

  ALIGNMENT_REPO = cache_dir.alignments.find
  Open.repository_dirs << ALIGNMENT_REPO unless Open.repository_dirs.include? ALIGNMENT_REPO

  def self.pdb_chain_position_in_sequence(pdb, pdbfile, chain, positions, protein_sequence)
    begin
      raise "Protein sequence missing" if protein_sequence.nil? or protein_sequence.empty?

      protein_sequence.gsub!(/\s/,'')
      protein_alignment, chain_alignment = Misc.insist do
        Persist.persist("SW PDB Alignment", :array,
                        :dir => ALIGNMENT_REPO, :persist => true,
                        :other => {:pdb => pdb, :pdbfile => pdbfile, :chain => chain, :protein_sequence => protein_sequence}) do
          chains = PDBHelper.pdb_chain_sequences(pdb, pdbfile)
          chain_sequence = chains[chain] 

          protein_alignment, chain_alignment =  Persist.persist("SmithWaterman alignment", :marshal, :dir => SSW_ALIGNMENT_REPO, :other => {:source => protein_sequence, :target => chain_sequence}) do 
            SmithWaterman.align(protein_sequence, chain_sequence)
          end
                        end
      end

      seq_pos = Structure.match_position(positions, chain_alignment, protein_alignment)
    rescue Exception
      Log.exception $!
      seq_pos = [nil] * positions.length
    end
    res = Hash[*positions.zip(seq_pos).flatten]
    TSV.setup(res, :key_field => "Chain position", :fields => ["Sequence position"], :type => :single, :unnamed => true)
  end

  def self.pdb_alignment_map(protein_sequence, pdb, pdbfile)
    protein_sequence.gsub!(/\s/,'')
    chains = PDBHelper.pdb_chain_sequences(pdb, pdbfile)

    result = TSV.setup({}, :key_field => "Sequence position", :fields => ["Chain:Position in PDB"], :type => :flat)
    chains.each do |chain,chain_sequence|
      chain_alignment, protein_alignment =  Persist.persist("SmithWaterman alignment", :marshal, :dir => SSW_ALIGNMENT_REPO, :other => {:source => chain_sequence, :target => protein_sequence}) do 
        SmithWaterman.align(chain_sequence, protein_sequence)
      end

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
    sequence.gsub!(/\s/,'')
    chain_positions = {}
    pdb_positions.collect do |cp|
      chain, position = cp.split(":")
      chain_positions[chain] ||= []
      chain_positions[chain] << position
    end
    return [] unless chain_positions.include? target_chain
    Structure.job(:pdb_chain_position_in_sequence, "TEST", :pdb => pdb, :pdbfile => pdbfile, :sequence => sequence, :chain => target_chain, :positions => chain_positions[target_chain]).exec
  end

  def self.sequence_position_in_pdb(protein_sequence, protein_positions, pdb, pdbfile)
    protein_sequence.gsub!(/\s/,'')
    chains = PDBHelper.pdb_chain_sequences(pdb, pdbfile)

    protein_positions = [protein_positions] unless Array === protein_positions
    alignments = {}
    chains.each do |chain,chain_sequence|
      alignments[chain] = begin
                           chain_alignment, protein_alignment =  Persist.persist("SmithWaterman alignment", :marshal, :update => true, :dir => SSW_ALIGNMENT_REPO, :other => {:source => chain_sequence, :target => protein_sequence}) do 
                             SmithWaterman.align(chain_sequence, protein_sequence)
                           end
                           Structure.match_position(protein_positions, protein_alignment, chain_alignment)
                          end
    end

    TSV.setup(alignments, :key_field => "PDB Chain", :fields => protein_positions, :type => :list, :cast => :to_i, :unnamed => true)
  end

  def self.neighbour_map_job(pdb, pdbfile, distance)
    Misc.insist do
      begin
        Persist.persist("Neighbour map", :marshal, :dir => NEIGHBOUR_MAP, :other => {:pdb => pdb, :pdbfile => pdbfile, :distance => distance}) do  |filename|
          job = Structure.job(:neighbour_map, "PDB Neighbours", :pdb => pdb, :pdbfile => pdbfile, :distance => distance)
          job.run
        end
      rescue Exception
        Log.warn "Exception calculating neighbour map: #{$!.message}. Retrying"
        persist = :update
        raise $!
      end
    end
  end

  def self.neighbours_in_pdb(sequence, positions, pdb = nil, pdbfile = nil, chain = nil, distance = 5)
    sequence.gsub!(/\s/,'')
    raise "No sequence positions specified" if positions.nil?
    neighbours_in_pdb = TSV.setup({}, :key_field => "Sequence position", :fields => ["Neighbours"], :type => :flat)

    positions_in_pdb = Structure.sequence_position_in_pdb(sequence, positions, pdb, pdbfile)

    Log.debug "Position in PDB: #{Misc.fingerprint positions_in_pdb}"

    chain ||=  positions_in_pdb.collect.reject{|c,p| p.nil? }.sort_by{|c,p| p.length}.first.first

    return neighbours_in_pdb if positions_in_pdb.nil? or positions_in_pdb[chain].nil?

    neighbour_map = neighbour_map_job(pdb,pdbfile,distance)
    return neighbours_in_pdb if neighbour_map.nil?

    inverse_neighbour_map = {}
    neighbour_map.each do |k,vs|
      vs.each do |v|
        inverse_neighbour_map[v] ||= []
        inverse_neighbour_map[v] << k
      end
    end

    positions_in_pdb[chain].each do |position|
      position_in_chain = [chain, position] * ":"
      neigh = neighbour_map[position_in_chain]
      ineigh = inverse_neighbour_map[position_in_chain]
      all_neighbours = (neigh || []) + (ineigh || []).uniq

      position_in_sequence = Structure.pdb_chain_position_in_sequence(pdb, pdbfile, chain, [position], sequence).values.compact.flatten.first
      next if position_in_sequence.nil?
      all_neighbours_in_sequence = all_neighbours.collect{|n| c,p = n.split(":"); Structure.pdb_chain_position_in_sequence(pdb, pdbfile, c, [p], sequence).values.compact.flatten }.flatten.uniq
      neighbours_in_pdb[position_in_sequence] = all_neighbours_in_sequence
    end

    neighbours_in_pdb
  end

  def self.neighbour_map(distance, pdb = nil, pdbfile = nil)
    return TSV.setup({}, :key_field => "Residue", :fields => ["Neighbours"], :type => :flat) if pdb == '1cw3'
    close_residues = PDBHelper.pdb_close_residues(distance, pdb, pdbfile)
    tsv = TSV.setup close_residues, :key_field => "Residue", :fields => ["Neighbours"], :type => :flat
    new = {}
    tsv.each do |p,ns|
      ns.each{|n| new[n] ||= []; new[n] << p }
    end
    new.each{|p,ns| tsv[p] ||= []; tsv[p].concat ns }
    tsv
  end

end

