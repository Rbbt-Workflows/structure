require 'rbbt'

module Interactome3d
  extend Resource
  self.subdir = 'share/databases/interactome3d'

  self.claim self.proteins_tsv, :proc do |filename|
    tsv = TSV.open('http://interactome3d.irbbarcelona.org/user_data/human/download/complete/proteins.dat', 
                   :header_hash => '',
                   :merge => true)
    tsv.to_s
  end

  self.claim self.interactions_tsv, :proc do |filename|
    tsv = TSV.open('http://interactome3d.irbbarcelona.org/user_data/human/download/complete/interactions.dat', 
                   :header_hash => '',
                   :merge => true)
    tsv.to_s
  end
end
