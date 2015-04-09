$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '../..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"
require 'structure/pdb_helper'

class TestClass < Test::Unit::TestCase
  def _test_pdb_atom_distance
    pdb = '3dzy'
    assert PDBHelper.pdb_atom_distance(5, pdb).length > 1
  end

  def test_missing_chain
    iii PDBHelper.pdb_chain_sequences(nil, Open.read('/home/mvazquezg/modbase_UNC5C_3g5bA.pdb'))
  end
end

