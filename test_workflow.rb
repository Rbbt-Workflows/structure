$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), 'lib'))

require 'rbbt'
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"

class TestStructure < Test::Unit::TestCase
  def test_neighbours_annotations
    mis = ["ENSP00000419692:S427V"]

    job =  Structure.job(:annotated_variant_neighbours, nil, :mutated_isoforms => mis, :database => "UniProt"); job.clean.run
    assert job.load.include? mis.first
    assert job.file(:residue_annotations).tsv["ENSP00000419692"]["Uniprot Features"].include? "HELIX"

    job =  Structure.job(:annotated_variant_neighbours, nil, :mutated_isoforms => mis, :database => "Appris"); job.clean.run
    assert job.load.include? mis.first
    assert job.file(:residue_annotations).tsv["ENSP00000419692"]["Appris Features"].include? "firestar"
  end
  
  def test_uniprot_residues
    r = TSV.setup({})
    r["ENSP00000391127"] = [110]
    job =  Structure.job(:annotate_residues_UniProt, nil, :residues => r); job.clean.run
    assert job.load.include? "ENSP00000391127"
  end

  def test_uniprot_variants
    r = TSV.setup({})
    r["ENSP00000236709"] = [218]
    job =  Structure.job(:annotate_residues_UniProt_variants, nil, :residues => r)
    assert job.load.include? "ENSP00000236709"
    assert job.load["ENSP00000236709"]["Residue"].include? 218
    assert job.load["ENSP00000236709"]["SNP ID"].include? "rs2246945"
  end

  def test_AAAS_residue_overlap
    job = Structure.example_step(:annotate_residues_COSMIC, "AAAS")
    job.clean.run
    unzipped = job.load.unzip
    assert unzipped.include? "ENSP00000446885:1"
    assert unzipped["ENSP00000446885:1"]["Genomic Mutation"].include? "12:53709535:T"
    assert(! unzipped["ENSP00000446885:95"]["Genomic Mutation"].include?("12:53709535:T"))
  end
end
