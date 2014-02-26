require 'rbbt/tsv'
require 'rbbt/sources/organism'
require 'structure/interactome_3d'
require 'structure/alignment'

module Structure
  ISO2UNI = Organism.protein_identifiers("Hsa").index :target => "UniProt/SwissProt Accession", :persist => true, :unnamed => true
  I3D_PROTEINS = Interactome3d.proteins_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS = Interactome3d.interactions_tsv.tsv :merge => true, :unnamed => true, :persist => true
  I3D_INTERACTIONS_REVERSE = Interactome3d.interactions_tsv.tsv :merge => true, :key_field => "PROT2", :zipped => true, :unnamed => true, :persist => true

  def self.neighbours_i3d(protein, positions, only_pdb = false)
    Log.info("PROCESSING #{Log.color :red, protein} -- #{Misc.fingerprint positions}")

    uniprot = ISO2UNI[protein]
    sequence = protein.sequence

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
        rescue Exception
          Log.warn "Error processing #{ url }: #{$!.message}"
          next
        end

        #Try another PDB unless at least one neighbour is found
        next if neighbours_in_pdb.nil? or neighbours_in_pdb.empty?

        tsv = TSV.setup({}, :key_field => "Isoform residue", :fields => ["PDB", "Neighbours"], :type => :list)
        neighbours_in_pdb.each do |pos, neigh|
          seq_pos = pdb_positions_to_sequence([pos], sequence, chain, url, nil)
          next if seq_pos.nil? or seq_pos.empty?
          seq_pos = seq_pos.first 
          seq_neigh = pdb_positions_to_sequence(neigh, sequence, chain, url, nil) 
          tsv[[protein, seq_pos] * ":"] = [url, seq_neigh * ";"]
        end

        return tsv 
      end
    else
      return nil if only_pdb
    end

    tsv = TSV.setup({}, :key_field => "Isoform residue", :fields => ["PDB", "Neighbours"])
    positions.each do |p| 
      new = []
      new << p-1 if p > 1
      new << p+1 if p< sequence.length 
      tsv[p] = [[nil], new * ";"]
    end

    tsv
  end












  def self.interface_neighbours_i3d(protein, positions)

    Log.info("PROCESSING #{Term::ANSIColor.red(protein)} -- #{Misc.fingerprint positions}")

    uniprot = ISO2UNI[protein]
    sequence = protein.sequence
 
    if uniprot and  I3D_INTERACTIONS.include? uniprot
      partner_interface = {}

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

          begin
            neighbours_in_pdb = self.neighbours_in_pdb(sequence, positions, url, nil, chain, 8).values.compact.flatten.uniq
          rescue Exception
            Log.exception $!
            next
          end

          next if neighbours_in_pdb.nil? or neighbours_in_pdb.empty?

          sequence_positions = pdb_positions_to_sequence(neighbours_in_pdb, partner_sequence, partner_chain, url, nil) 

          next unless (sequence_positions and sequence_positions.any?)

          partner_interface[partner_ensembl] ||= []
          partner_interface[partner_ensembl].concat sequence_positions
          partner_interface[partner_ensembl].uniq!

        end
      end
      partner_interface
    else
      return {} 
    end
  end
end
