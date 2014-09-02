require 'rbbt/entity'
Workflow.require_workflow 'Genomics'
require 'rbbt/entity/mutated_isoform'
require 'rbbt/entity/protein'

module MutatedIsoform
  property :pdbs_and_positions => :single do
    return [] if pdbs.nil?
    pdbs.collect do |pdb, info|
      [pdb, Structure.job(:sequence_position_in_pdb, "Protein: #{ self }", :sequence => protein.sequence, :organism => organism, :position => position, :pdb => pdb).run]
    end
  end
end

module Protein
  extend Entity
  self.format = ["Partner Ensembl Protein ID"]

  property :name => :single do
    "(#{gene.name}) " << self
  end
end
