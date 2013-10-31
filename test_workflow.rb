$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), 'lib'))

require 'rbbt'
require 'test/unit'
require 'rbbt/util/log'
require 'rbbt/util/tmpfile'
require 'rbbt/resource/path'
require 'fileutils'

class Test::Unit::TestCase
  include FileUtils

  #def setup
  #  Persist.cachedir = Rbbt.tmp.test.persistence.find :user
  #end

  #def teardown
  #  FileUtils.rm_rf Path.setup("", 'rbbt').tmp.test.find :user
  #  Persist::TC_CONNECTIONS.values.each do |c| c.close end
  #  Persist::TC_CONNECTIONS.clear
  #end

  def datafile_test(file)
    File.join(File.dirname(__FILE__), 'data', file)
  end
end

require 'rbbt/workflow'
Workflow.require_workflow 'workflow'

Rbbt.share.databases.CATH.CathRegions.produce

class TestFixWidthTable < Test::Unit::TestCase

  def _test_protein_aa_pos_to_pdb
    pdbs = Structure.pdbs_covering_aa_position('P04637', 350).collect{|pdb,info| pdb}

    cath = pdbs.collect do |pdb|
      pdb = pdb.downcase
      structure =  Uniprot.cath(pdb)
      next if structure.nil?
      structure["CATH Code"]
    end.compact

    assert_equal  ["4.10.170"], cath
  end

  def _test_protein_cath_domains
    assert_equal ["1olgA00", "2ioiA00"], Uniprot.cath_domains('P04637')
  end

  def test_protein_cath_domains_at_aa_position
    ddd Structure.cath_domains_at_aa_position('P04637', 350)
  end
end
