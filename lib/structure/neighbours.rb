require 'rbbt/tsv'
require 'rbbt/sources/organism'
require 'structure/identifiers'
require 'structure/interactome_3d'
require 'structure/alignment'

module Structure
  I3D_PROTEINS = Interactome3d.proteins_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS = Interactome3d.interactions_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS_REVERSE = Interactome3d.interactions_tsv.tsv :merge => true, :key_field => "PROT2", :zipped => true, :unnamed => true, :persist => true

  def self.neighbours_uniprot(protein, positions, organism = Organism.default_code("Hsa"), only_pdb = false)
    uniprot = iso2uni(organism)[protein]
    sequence = iso2seq(organism)[protein]

    tsv = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)
    return tsv if sequence.nil?

    uni_pdbs = UniProt.pdbs(uniprot) unless uniprot.nil?

    if uniprot and  uni_pdbs.any?

      uni_pdbs.each do |pdb, info|

        # Be more careful with identifying positions in the wrong chain do to
        # aberrant alignments

        neighbours_in_pdb = nil
        begin
          neighbours_in_pdb = self.neighbours_in_pdb(sequence, positions, pdb, nil, nil, 5) 
        rescue 
          Log.warn "Error processing #{ pdb }: #{$!.message}"
          next
        end

        #Try another PDB unless at least one neighbour is found
        next if neighbours_in_pdb.nil? or neighbours_in_pdb.empty?

        neighbours_in_pdb.each do |seq_pos, seq_neigh|
          next if seq_pos.nil?
          tsv[[protein, seq_pos] * ":"] = [protein, seq_pos, pdb, seq_neigh * ";"]
        end

        return tsv if tsv.any?
      end
    else
      return tsv if only_pdb
    end

    # TODO Add contiguous also when PDB is present
    positions.each do |p| 
      new = []
      new << p-1 if p > 1
      new << p+1 if p < sequence.length 
      tsv[[protein,p]*":"] = [protein, p, nil, new * ";"]
    end

    tsv
  end

  def self.neighbours_i3d(protein, positions, organism = Organism.default_code("Hsa"), only_pdb = false)

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

        neighbours_in_pdb = nil
        begin
          neighbours_in_pdb = self.neighbours_in_pdb(sequence, positions, url, nil, chain, 5) 
        rescue 
          Log.warn "Error processing #{ url }: #{$!.message}"
          next
        end

        #Try another PDB unless at least one neighbour is found
        next if neighbours_in_pdb.nil? or neighbours_in_pdb.empty?

        neighbours_in_pdb.each do |seq_pos, seq_neigh|
          next if seq_pos.nil?
          tsv[[protein, seq_pos] * ":"] = [protein, seq_pos, url, seq_neigh * ";"]
        end

        return tsv if tsv.any?
      end
    end

    return tsv if only_pdb

    # TODO Add contiguous also when PDB is present
    positions.each do |p| 
      new = []
      new << p-1 if p > 1
      new << p+1 if p < sequence.length 
      tsv[[protein,p]*":"] = [protein, p, nil, new * ";"]
    end

    tsv
  end

  def self.neighbours(protein, positions, organism = Organism.default_code("Hsa"), only_pdb = false)
    tsv = self.neighbours_i3d(protein, positions, organism, true)
    return tsv unless tsv.empty?
    self.neighbours_uniprot(protein, positions, organism, only_pdb)
  end

  def self.interface_neighbours_i3d(protein, positions, organism = Organism.default_code("Hsa"), distance = 8)
    tsv = TSV.setup({}, :key_field => "Isoform", :fields => ["Position", "Partner Ensembl Protein ID", "PDB", "Partner residues"], :type => :double)

    uniprot = iso2uni(organism)[protein]
    return tsv if uniprot.nil?
    sequence = iso2seq(organism)[protein]
    return tsv if sequence.nil?

    uni2iso = uni2iso(organism)

    forward_positions = ["PDB_ID", "PROT2", "CHAIN1", "CHAIN2", "FILENAME"].collect{|f| I3D_INTERACTIONS.identify_field f}
    reverse_positions = ["PDB_ID", "PROT1", "CHAIN2", "CHAIN1", "FILENAME"].collect{|f| I3D_INTERACTIONS_REVERSE.identify_field f}

    {:forward => I3D_INTERACTIONS,
      :reverse => I3D_INTERACTIONS_REVERSE}.each do |db_direction, db|

      next unless db.include? uniprot
      values_list = db[uniprot]
      seen_pbs = []
      Misc.zip_fields(values_list).each do |values|
        if db_direction == :forward
          pdb, partner, orig_chain, orig_partner_chain, filename = values.values_at *forward_positions
          chain, partner_chain = "A", "B"
        else
          pdb, partner, orig_chain, orig_partner_chain, filename = values.values_at *reverse_positions
          chain, partner_chain = "B", "A"
        end

        next if chain.strip.empty? or partner_chain.strip.empty?

        if uniprot == partner
          partner_ensembl = protein
        else
          partner_ensembl =  uni2iso[partner]
        end

        if partner_ensembl.nil?
          Log.warn "Could not translate partner to Ensembl: #{ partner }"
          next
        end

        partner_sequence = iso2seq[partner_ensembl]
        if partner_sequence.nil?
          Log.warn "Could get partner sequence: #{ partner } (#{partner_ensembl})"
          next
        end


        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=interactions&type2=#{ type }&pdb=#{ filename }"
        Log.debug "Processing: #{ url }"

        positions_in_pdb = sequence_position_in_pdb(sequence, positions, url, nil)[chain]
        next if positions_in_pdb.nil? or positions_in_pdb.empty?

        map = neighbour_map_job url, nil, distance
        #map.unnamed = true

        next if map.nil? or map.empty?

        positions_in_pdb.zip(positions).each do |pdb_position, position|
          code = [chain,pdb_position]*":"
          begin
            neighbours = map[code]
          rescue
            Log.exception $!
            next
          end
          next if neighbours.nil? or neighbours.empty?
          partner_neighbours = neighbours.select{|c| c.split(":").first == partner_chain }.collect{|c| c.split(":").last}
          next if partner_neighbours.empty?
          partner_neighbours_in_sequence = pdb_chain_position_in_sequence(url, nil, partner_chain, partner_neighbours, partner_sequence).values.compact.flatten
          next if partner_neighbours_in_sequence.empty?
          tsv.zip_new(protein, [position, partner_ensembl, url, partner_neighbours_in_sequence * ";"])
        end

        if not seen_pbs.include? pdb
          next if orig_chain == orig_partner_chain
          seen_pbs << pdb
          url = pdb
          Log.debug "Processing: #{ url }"

          positions_in_pdb = sequence_position_in_pdb(sequence, positions, url, nil)[orig_chain]
          next if positions_in_pdb.nil? or positions_in_pdb.empty?

          map = neighbour_map_job url, nil, distance
          #map.unnamed = true

          next if map.nil? or map.empty?

          positions_in_pdb.zip(positions).each do |pdb_position, position|
            code = [orig_chain,pdb_position]*":"
            begin
              neighbours = map[code]
            rescue
              Log.exception $!
              next
            end
            next if neighbours.nil? or neighbours.empty?
            partner_neighbours = neighbours.select{|c| c.split(":").first == orig_partner_chain }.collect{|c| c.split(":").last}
            next if partner_neighbours.empty?
            partner_neighbours_in_sequence = pdb_chain_position_in_sequence(url, nil, orig_partner_chain, partner_neighbours, partner_sequence).values.compact.flatten
            next if partner_neighbours_in_sequence.empty?
            tsv.zip_new(protein, [position, partner_ensembl, url, partner_neighbours_in_sequence * ";"])
          end
        end
      end
    end

    tsv
  end
end
