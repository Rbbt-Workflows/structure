$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '../..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'structure/pdb_helper'

class TestClass < Test::Unit::TestCase
  def test_pdb_atom_distance
    pdb = '3dzy'
    assert PDBHelper.pdb_atom_distance(5, pdb).length > 1
  end
end

