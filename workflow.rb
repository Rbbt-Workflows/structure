require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'

require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/InterPro'

require 'structure/alignment'
require 'structure/pdb_alignment'

require 'structure/pdb_helper'
require 'structure/neighbours'
require 'structure/interactome_3d'

require 'structure/uniprot'
require 'structure/appris'
require 'structure/COSMIC'
require 'structure/interpro'

Workflow.require_workflow 'Genomics'
Workflow.require_workflow 'Translation'
Workflow.require_workflow 'PdbTools'
Workflow.require_workflow "Sequence"

require 'rbbt/entity'
require 'rbbt/entity/mutated_isoform'
require 'structure/entity'

require 'rbbt/util/simpleopt'
$cpus ||= SOPT.get("--cpus* CPUs to use in map-reduce (Structure workflow)")[:cpus]
$cpus = $cpus.to_i if $cpus

Log.info "Loading Structure with #{ $cpus.inspect }" unless $cpus.nil?

module Structure
  extend Workflow
end

require 'structure/workflow/alignments'
require 'structure/workflow/helpers'
require 'structure/workflow/annotate'
