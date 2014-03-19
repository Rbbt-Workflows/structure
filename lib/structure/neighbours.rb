require 'rbbt/tsv'
require 'rbbt/sources/organism'
require 'structure/interactome_3d'
require 'structure/alignment'

module Structure
  #ISO2UNI = Organism.protein_identifiers("Hsa").index :target => "UniProt/SwissProt Accession", :persist => true, :unnamed => true
  #ISO2SEQ = Organism.protein_sequence("Hsa").tsv :persist => true, :unnamed => true
  I3D_PROTEINS = Interactome3d.proteins_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS = Interactome3d.interactions_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS_REVERSE = Interactome3d.interactions_tsv.tsv :merge => true, :key_field => "PROT2", :zipped => true, :unnamed => true, :persist => true

  def self.iso2uni(organism = "Hsa")
    @@iso2uni ||= {}
    @@iso2uni[organism] ||= Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true, :unnamed => true
  end

  def self.iso2seq(organism = "Hsa")
    @@iso2seq ||= {}
    @@iso2seq[organism] ||= Organism.protein_sequence(organism).tsv :persist => true, :unnamed => true
  end

  def self.neighbours_i3d(protein, positions, organism = "Hsa", only_pdb = false)

    uniprot = iso2uni(organism)[protein]
    sequence = iso2seq(organism)[protein]

    tsv = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)
    return tsv if sequence.nil?

    if uniprot and  I3D_PROTEINS.include? uniprot

      field_positions = ["CHAIN", "FILENAME"].collect{|f| I3D_PROTEINS.identify_field f}
      Misc.zip_fields(I3D_PROTEINS[uniprot]).each do |values|

        # Be more careful with identifying positions in the wrong chain do to
        # aberrant alignments
        chain, filename = values.values_at *field_positions
        next if chain.strip.empty?

        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=proteins&type2=#{ type }&pdb=#{ filename }"

        begin
          neighbours_in_pdb = self.neighbours_in_pdb(sequence, positions, url, nil, chain, 5) 
        rescue 
          Log.warn "Error processing #{ url }: #{$!.message}"
          next
        end

        #Try another PDB unless at least one neighbour is found
        next if neighbours_in_pdb.nil? or neighbours_in_pdb.empty?

        neighbours_in_pdb.each do |seq_pos, seq_neigh|
          tsv[[protein, seq_pos] * ":"] = [protein, seq_pos, url, seq_neigh * ";"]
        end

        return tsv 
      end
    else
      return tsv if only_pdb
    end

    positions.each do |p| 
      new = []
      new << p-1 if p > 1
      new << p+1 if p < sequence.length 
      tsv[[protein,p]*":"] = [protein, p, [nil], new * ";"]
    end

    tsv
  end

  def self.interface_neighbours_i3d(protein, positions, organism = "Hsa")
    tsv = TSV.setup({}, :key_field => "Isoform", :fields => ["Position", "Partner Ensembl Protein ID", "PDB", "Partner residues"], :type => :double)

    uniprot = iso2uni(organism)[protein]
    return tsv if uniprot.nil?
    sequence = iso2seq(organism)[protein]
    return tsv if sequence.nil?

    Protein.setup(protein, "Ensembl Protein ID", organism)

    forward_positions = ["PROT2", "CHAIN1", "CHAIN2", "FILENAME"].collect{|f| I3D_INTERACTIONS.identify_field f}
    reverse_positions = ["PROT1", "CHAIN2", "CHAIN1", "FILENAME"].collect{|f| I3D_INTERACTIONS_REVERSE.identify_field f}

    {:forward => I3D_INTERACTIONS,
      :reverse => I3D_INTERACTIONS_REVERSE}.each do |db_direction, db|

      next unless db.include? uniprot
      Misc.zip_fields(db[uniprot]).each do |values|

        if db_direction == :forward
          partner, chain, partner_chain, filename = values.values_at *forward_positions
          chain, partner_chain = "A", "B"
        else
          partner, chain, partner_chain, filename = values.values_at *reverse_positions
          chain, partner_chain = "B", "A"
        end

        next if chain.strip.empty? or partner_chain.strip.empty?

        Protein.setup(partner, "UniProt/SwissProt Accession", protein.organism)
        partner_sequence = partner.sequence
        partner_ensembl =  partner.ensembl

        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=interactions&type2=#{ type }&pdb=#{ filename }"
        Log.debug "Processing: #{ url }"

        positions_in_pdb = sequence_position_in_pdb(sequence, positions, url, nil)[chain]
        next if positions_in_pdb.nil? or positions_in_pdb.empty?

        begin
          map = Structure.job(:neighbour_map, protein, :distance =>  8, :pdb => url).run 
        rescue 
          Log.warn "Error processing #{ url }: #{$!.message}"
          next
        end

        positions_in_pdb.zip(positions).each do |pdb_position, position|
          code = [chain,pdb_position]*":"
          neighbours = map[code]
          next if neighbours.nil? or neighbours.empty?
          partner_neighbours = neighbours.select{|c| c.split(":").first == partner_chain }.collect{|c| c.split(":").last}
          next if partner_neighbours.empty?
          partner_neighbours_in_sequence = pdb_chain_position_in_sequence(url, nil, partner_chain, partner_neighbours, partner_sequence).values.compact.flatten
          next if partner_neighbours_in_sequence.empty?
          tsv.zip_new(protein, [position, partner_ensembl, url, partner_neighbours_in_sequence * ";"])
        end
      end
    end
    tsv
  end
end
