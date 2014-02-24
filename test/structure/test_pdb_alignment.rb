$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"
require 'structure/ssw'
require 'structure/alignment'
require 'awesome_print'


class TestClass < Test::Unit::TestCase
  def test_sequence_position_in_pdb
    protein = "ENSP00000419692"
    Protein.setup(protein, "Ensembl Protein ID", "Hsa")

    protein_sequence = protein.sequence
    position = 427
    pdb_url =  "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=interactions&type2=model&pdb=O75469-P19793-MDL-1xvp.pdb3-D-0-C-0.pdb"

    assert_equal [487], Structure.sequence_position_in_pdb(protein_sequence, [position], pdb_url, nil)["B"]
  end
end

