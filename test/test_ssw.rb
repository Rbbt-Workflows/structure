$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"
require 'ssw'


class TestClass < Test::Unit::TestCase
  def test_alignment_map
    protein = Protein.setup("ENSP00000308495", "Ensembl Protein ID", "Hsa")
    uniprot = "G3V5T7"

    uniprot_sequence = UniProt.sequence(uniprot)
    uniprot_alignment, original_alignment = SmithWaterman.align(uniprot_sequence, protein.sequence)
    map = Structure.alignment_map(uniprot_alignment, original_alignment)

    assert_equal 0, map[0]
    assert_equal protein.sequence.sub(/\*$/,'').length-1, map[uniprot_sequence.length-1]
  end
end

