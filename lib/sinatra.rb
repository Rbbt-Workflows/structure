require 'rbbt/rest/web_tool'
include Sinatra::RbbtToolHelper

require 'rbbt-util'
require 'rbbt/workflow'

require 'rbbt/entity/structure'

require 'rbbt/entity/gene'
require 'rbbt/entity/mutated_isoform'
require 'rbbt/entity/InterPro'

[Gene, Protein, MutatedIsoform, InterPro].each do |mod|
  mod.instance_eval do
    include Entity::REST
  end
end

Workflow.require_workflow "Appris"

require 'rbbt/entity/appris'

$title = "Structure PPI"
