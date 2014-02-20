
module Structure
  ISO2UNI = Organism.protein_identifiers("Hsa").index :target => "UniProt/SwissProt Accession", :persist => true
  I3D_PROTEINS = Interactome3d.proteins_tsv.tsv :merge => true

  def self.pdb_position_to_sequence(neighbours, sequence, pdb = nil, pdbfile = nil)
    neighbours.collect do |cp|
      chain, position = cp.split(":")
      Structure.job(:pdb_chain_position_in_sequence, "TEST", :pdb => pdb, :pdbfile => pdbfile, :sequence => sequence, :chain => chain, :positions => [position.to_i]).run
    end.compact.flatten
  end

  def self.neighbours_i3d(protein, positions)
    Log.info("PROCESSING #{Term::ANSIColor.red(protein)} -- #{Misc.fingerprint positions}")
    uniprot = ISO2UNI[protein]
    return nil if uniprot.nil?
    return nil unless I3D_PROTEINS.include? uniprot
    sequence = protein.sequence

    I3D_PROTEINS[uniprot].zip_fields.each do |values|
      #ToDo: we are looking only for chains specified by the source file. Some
      #are missing and some a wrong
      
      chain, filename = values.values_at "CHAIN", "FILENAME"
      next if chain.strip.empty?

      type = filename =~ /EXP/ ? :pdb : :model
      url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=proteins&type2=#{ type }&pdb=#{ filename }"
      pdbfile = Open.read(url)

      positions_in_pdb = Structure.job(:sequence_position_in_pdb, "TEST", :pdbfile => pdbfile, :sequence => sequence, :positions => positions).run
      next if positions_in_pdb.empty?

      neighbour_map = Persist.persist("PDB Neighbours: #{url}", :marshal, :pdb => url, :dir => Rbbt.var.persist.find(:lib)) do
        PDBHelper.pdb_close_residues(4, nil, pdbfile)
      end


      next if positions_in_pdb[chain].nil?
      neighbours_in_pdb = positions_in_pdb[chain].collect do |position|
        position_in_chain = [chain, position] * ":"
        neighbour_map[position_in_chain]
      end.compact.flatten

      #Try another PDB unless at least one neighbour is found
      next if neighbours_in_pdb.empty?

      sequence_positions = pdb_position_to_sequence(neighbours_in_pdb, sequence, nil, pdbfile) 

      return sequence_positions 
    end
    nil
  end

  def self.neighbours(*args)
    neighbours_i3d(*args)
  end
end
