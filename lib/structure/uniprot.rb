module Structure
  def self.uniprot_sequence_map(uniprot, sequence)
    uniprot_sequence = UniProt.sequence(uniprot)
    map = sequence_map(uniprot_sequence, sequence)
  end

  def self.corrected_uniprot_features(uniprot, sequence)
    uniprot_sequence = UniProt.sequence(uniprot)

    map = uniprot_sequence_map(uniprot, sequence)

    features = UniProt.features(uniprot)
    corrected_features = []
    features.each do |info|
      corrected_start = map[info[:start]]
      corrected_end = map[info[:end]]
      next if corrected_start.nil? or corrected_end.nil?
      corrected_info = info.dup
      corrected_info[:start] = corrected_start
      corrected_info[:end] = corrected_end
      corrected_features << corrected_info
    end

    corrected_features
  end


  def self.UniProt_residues
    @UniProt_residues ||= Persist.persist_tsv(UniProt.annotated_variants, "UniProt::residues", {}, :persist => true, :serializer => :list, :dir => Rbbt.var.persistence.find(:lib)) do |data|
                           isoform_residue_mutations = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["UniProt Variant ID"], :type => :flat)

                           uni2ensp = Organism.protein_identifiers(Organism.default_code("Hsa")).tsv :fields => ["Ensembl Protein ID"], :key_field => "UniProt/SwissProt Accession", :persist => true, :type => :flat, :merge => true, :unnamed => true
                           ensp2sequence = Organism.protein_sequence(Organism.default_code("Hsa")).tsv :persist => true, :unnamed => true

                           db = UniProt.annotated_variants.tsv(:fields => ["Amino Acid Mutation", "UniProt Variant ID"], :persist => true, :type => :double, :unnamed => true)
                           db.monitor = {:desc => "Processing UniProt", :step => 1000}

                           db.with_unnamed do
                             db.through do |uniprot,values|
                               Log.debug Log.color :red, uniprot
                               begin
                                 ensps = uni2ensp[uniprot]
                                 raise "No translation to Ensembl: #{ uniprot }" if ensps.nil? or ensps.empty?
                                 found = false
                                 ensps.each do |ensp|
                                   Log.debug Log.color :blue, ensp
                                   begin
                                     ensp_sequence = ensp2sequence[ensp]
                                     raise "No sequence: #{ ensp } " if ensp_sequence.nil?
                                     uniprot_sequence = UniProt.sequence(uniprot)
                                     map = Structure.sequence_map(uniprot_sequence, ensp_sequence)

                                     Misc.zip_fields(values).each do |change,vid|
                                       match = change.match(/^([A-Z])(\d+)([A-Z])$/)
                                       raise "Unknown change: #{ ensp } #{change}" if match.nil?
                                       ref, _pos, mut = match.values_at 1,2,3
                                       pos = map[_pos.to_i]
                                       raise "Unmapped position: #{ ensp } #{_pos}" if pos.nil?
                                       isoform_residue_mutations[[ensp,pos]*":"] ||= []
                                       isoform_residue_mutations[[ensp,pos]*":"] << vid
                                       found = true
                                     end
                                   rescue
                                     Log.debug $!.message
                                     next
                                   end
                                 end
                                 raise "No suitable mapings for #{ uniprot }" unless found
                               rescue
                                 Log.warn $!.message
                                 next
                               end
                             end
                           end

                           Log.info "Merging data"
                           data.merge! isoform_residue_mutations
                           Log.info "Annotating database"
                           isoform_residue_mutations.annotate data

                           data
    end
  end

  def self.UniProt_mutation_annotations
    @UniProt_mutation_annotations ||= begin
                                        fields = [
                                          'Amino Acid Mutation',
                                          'Type of Variant',
                                          'SNP ID',
                                          'Disease'
                                        ]
                                       UniProt.annotated_variants.tsv(:key_field => "UniProt Variant ID", :fields => fields, :persist => true, :type => :double, :zipped => true, :unnamed => true).to_list
                                     end
  end
end
