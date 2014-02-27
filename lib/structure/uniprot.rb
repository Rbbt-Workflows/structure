module Structure
  def self.corrected_uniprot_features(uniprot, sequence)
    uniprot_sequence = UniProt.sequence(uniprot)

    map = sequence_map(uniprot_sequence, sequence)

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

                           uni2ensp = Organism.protein_identifiers("Hsa").index :target => "Ensembl Protein ID", :fields => ["UniProt/SwissProt Accession"], :persist => true
                           ensp2sequence = Organism.protein_sequence("Hsa").tsv :persist => true

                           db = UniProt.annotated_variants.tsv(:fields => ["Amino Acid Mutation", "UniProt Variant ID"], :persist => true, :type => :double)
                           db.monitor = {:desc => "Processing UniProt", :step => 10000}

                           db.with_unnamed do
                             db.through do |uniprot,values|
                               ensp = uni2ensp[uniprot]
                               next if ensp.nil?
                               ensp_sequence = ensp2sequence[ensp]
                               next if ensp_sequence.nil?
                               uniprot_sequence = UniProt.sequence(uniprot)
                               map = Structure.sequence_map(uniprot_sequence, ensp_sequence)

                               Misc.zip_fields(values).each do |change,vid|
                                 match = change.match(/^([A-Z])(\d+)([A-Z])$/)
                                 next if match.nil?
                                 ref, pos, mut = match.values_at 1,2,3
                                 pos = map[pos.to_i]
                                 next if pos.nil?
                                 isoform_residue_mutations[[ensp,pos]*":"] ||= []
                                 isoform_residue_mutations[[ensp,pos]*":"] << vid
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
                                       UniProt.annotated_variants.tsv(:key_field => "UniProt Variant ID", :fields => fields, :persist => true, :type => :double).to_list
                                     end
  end
end
