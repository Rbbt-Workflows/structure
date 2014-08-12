Workflow.require_workflow "InterPro"
require 'rbbt/sources/InterPro'

module Structure
  def self.interpro_protein_domains
    @interpro_protein_domains ||= InterPro.protein_domains.tsv :persist => true, :unnamed => true
  end

  def self.corrected_interpro_features(uniprot, sequence)
    uniprot_sequence = UniProt.sequence(uniprot)

    map = sequence_map(uniprot_sequence, sequence)

    features = interpro_protein_domains[uniprot]
    return [] if features.nil? 
    corrected_features = []
    Misc.zip_fields(features).each do |code,start,eend|
      corrected_start = map[start.to_i]
      corrected_end = map[eend.to_i]
      next if corrected_start.nil? or corrected_end.nil?
      corrected_info = {}
      corrected_info[:start] = corrected_start
      corrected_info[:end] = corrected_end
      corrected_info[:code] = code
      corrected_features << corrected_info
    end

    corrected_features
  end
end
