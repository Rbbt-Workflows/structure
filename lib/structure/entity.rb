if defined? Entity and defined? MutatedIsoform and Entity === MutatedIsoform
  module MutatedIsoform
    property :pdbs_and_positions => :single do
      return [] if pdbs.nil?
      pdbs.collect do |pdb, info|
        [pdb, Structure.job(:sequence_position_in_pdb, "Protein: #{ self }", :sequence => protein.sequence, :organism => organism, :position => position, :pdb => pdb).run]
      end
    end
  end
end
