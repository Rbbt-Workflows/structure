$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"
require 'structure/uniprot'


class TestClass < Test::Unit::TestCase
  def test_corrected_features
    protein = Protein.setup("ENSP00000308495", "Ensembl Protein ID", "Hsa")
    uniprot = "G3V5T7"
    assert Structure.corrected_uniprot_features(uniprot, protein.sequence)
  end
end

