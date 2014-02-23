$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"
require 'structure/ssw'
require 'structure/alignment'
require 'awesome_print'


class TestClass < Test::Unit::TestCase
  def test_alignment_map
    ensp = Protein.setup("ENSP00000308495", "Ensembl Protein ID", "Hsa")
    uniprot = "G3V5T7"

    uniprot_sequence = UniProt.sequence(uniprot)
    ensp_sequence = ensp.sequence
    uniprot_alignment, ensp_alignment = SmithWaterman.align(uniprot_sequence, ensp_sequence)

    map = Structure.alignment_map(uniprot_alignment, ensp_alignment)
    
    map.each do |u,e|
      assert_equal ensp_sequence[e-1], uniprot_sequence[u-1]
    end
  end
end

