$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), 'lib'))

require 'rbbt'
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Structure"

class TestStructure < Test::Unit::TestCase
  def test_interfaces
    file = "/home/mvazquezg/test/mis.10000"
    mis = Open.read(file).split("\n")
    Log.severity = 4
    bar = Log::ProgressBar.new
    TSV.traverse mis, :threads => 3 do |mi|
      bar.tick
      next unless mi =~ /(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
      isoform = $1
      residue = $3.to_i
      Structure.interface_neighbours_i3d(isoform, [residue], "Hsa")
    end
  end
end
