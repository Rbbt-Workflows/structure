
Workflow.require_workflow "COSMIC"

module Structure
  def self.COSMIC_residues
    @COSMIC_residues ||= Persist.persist_tsv(nil, "COSMIC::residues", {}, :persist => true, :serializer => :list, :dir => Rbbt.var.persistence.find(:lib)) do |data|
                           isoform_residue_mutations = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Genomic Mutations"], :type => :flat)

                           db = COSMIC.knowledge_base.get_database('mutation_isoforms')

                           db.monitor = {:desc => "Processing COSMIC", :step => 10000}
                           db.with_unnamed do
                            db.through  do |mutation, mis|
                              protein_residues = {}
                              mis.flatten.each do |mi|
                                next unless mi =~ /(ENSP\d+):([A-Z])(\d+)([A-Z])$/ and $2 != $4
                                residue = $3.to_i
                                protein = $1
                                isoform_residue_mutations[[protein, residue] * ":"] ||= []
                                isoform_residue_mutations[[protein, residue] * ":"] << mutation
                              end
                            end
                           end

                           data.merge! isoform_residue_mutations
                           isoform_residue_mutations.annotate data

                           data
                         end
  end

  def self.COSMIC_mutation_annotations
    @COSMIC_mutation_annotations ||= begin
                                       fields = [
                                         'Sample name',
                                         'Primary site',
                                         'Site subtype',
                                         'Primary histology',
                                         'Histology subtype',
                                         'Pubmed_PMID',
                                         ]
                                       COSMIC.mutations.tsv :key_field => "Genomic Mutation", :fields => fields, :persist => true, :unamed => true, :type => :double, :zipped => true
                                     end
  end
end
