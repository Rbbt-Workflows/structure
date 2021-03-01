require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'

require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'

require 'structure/uniprot'
require 'structure/appris'
require 'structure/COSMIC'
require 'structure/interpro'


Workflow.require_workflow 'Translation'
#Workflow.require_workflow 'PdbTools'
Workflow.require_workflow "Sequence"
#Workflow.require_workflow "SphereClustering"
#
module Structure
  extend Workflow

  class << self
    attr_accessor :cache_dir
    def cache_dir
      @cache_dir ||= Rbbt.var.Structure
    end
  end

  extend Resource
  self.subdir = Structure.cache_dir
  Structure.claim Structure.root, :proc do |dirname|
    FileUtils.mkdir_p dirname
    nil
  end
end

require 'structure/alignment'
require 'structure/pdb_alignment'

require 'structure/pdb_helper'
require 'structure/neighbours'
require 'structure/interactome_3d'

require 'rbbt/util/simpleopt'
$cpus ||= SOPT.get("--cpus* CPUs to use in map-reduce (Structure workflow)")[:cpus]
$cpus = $cpus.to_i if $cpus

Log.info "Loading Structure with #{ $cpus.inspect }" unless $cpus.nil?

require 'structure/workflow/alignments'
require 'structure/workflow/helpers'
require 'structure/workflow/annotate'
require 'structure/workflow/wizard'
