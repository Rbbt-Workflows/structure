$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '../..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))

require 'rbbt-util'
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Genomics"
Workflow.require_workflow "Structure"
require 'structure/neighbours'


class TestClass < Test::Unit::TestCase
  def test_close_by
    protein = "ENSP00000386794"
    Protein.setup(protein, "Ensembl Protein ID", "Hsa")
    assert_equal [12,14,15,17,189,191], Structure.neighbours_i3d(protein, [13,16,190])
  end

  def test_neighbours_i3d
    protein = "ENSP00000251849"
    Protein.setup(protein, "Ensembl Protein ID", "Hsa")

    assert_equal [258,260], Structure.neighbours_i3d(protein, [259])
  end


  def test_interface_i3d
    protein = "ENSP00000419692"
    Protein.setup(protein, "Ensembl Protein ID", "Hsa")
    position = 427

    partner =  "ENSP00000336528"
    Protein.setup(partner, "Ensembl Protein ID", "Hsa")


    Log.debug "AA at #{ position }: #{protein.sequence[position-1]}"

    require 'pp'
    res = Structure.interface_neighbours_i3d(protein, [position])
    assert res.include? partner
    assert_equal ["D", "R", "R"], res[partner].sort.collect{|pos| partner.sequence[pos-1] }
  end
end

