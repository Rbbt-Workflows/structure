$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), 'lib'))

require 'rbbt'
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"

class TestStructure < Test::Unit::TestCase
  def test_neighbours_annotations
    mis = ["ENSP00000419692:S427V"]
    job =  Structure.job(:mutated_isoform_neighbour_annotation, nil, :mutated_isoforms => mis)
    job.clean
    job.run

    assert job.load.include? "ENSP00000419692:427"

    assert job.file(:neighbour_uniprot).tsv["ENSP00000419692"]["Uniprot Features"].include? "HELIX"
    assert job.file(:neighbour_appris).tsv["ENSP00000419692"]["Appris Features"].include? "firestar"

    assert job.file(:neighbour_appris).tsv.unzip(0, true)["ENSP00000419692:432"]["Appris Features"].include? 'firestar'
  end
  
  def _test_uniprot_variants
    r = TSV.setup({})
    r["ENSP00000236709"] = [218]
    job =  Structure.job(:annotate_variants_UNIPROT, nil, :residues => r)
    job.clean.run
    assert job.load.include? "ENSP00000236709"
    assert job.load["ENSP00000236709"]["Residue"].include? "218"
    assert job.load["ENSP00000236709"]["SNP ID"].include? "rs2246945"
  end
end
