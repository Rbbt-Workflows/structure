require 'rbbt/rest/web_tool'
include Sinatra::RbbtToolHelper

require 'rbbt-util'
require 'rbbt/workflow'

require 'rbbt/entity/structure'

require 'rbbt/entity/gene'
require 'rbbt/entity/mutated_isoform'
require 'rbbt/entity/InterPro'
require 'rbbt/sources/pfam'

[Gene, Protein, MutatedIsoform, InterProDomain, PfamDomain].each do |mod|
  mod.instance_eval do
    include Entity::REST
  end
end

Workflow.require_workflow "Appris"

require 'rbbt/entity/appris'

post '/wizard' do
  template_render('wizard', params, "Wizard", :cache_type => :sync)
end

get '/get_pdb' do
  query = consume_parameter :query
  content_type 'text/plain'
  send_file Open.open(query)
end
