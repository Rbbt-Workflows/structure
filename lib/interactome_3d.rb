require 'rbbt'

module Interactome3d
  extend Resource
  self.subdir = 'share/databases/interactome3d'

  self.claim self.proteins_tsv, :proc do |filename|
    tsv = TSV.open('http://interactome3d.irbbarcelona.org/user_data/human/download/complete/proteins.dat', 
                   :header_hash => '',
                   :merge => true)
    tsv.to_s
  end

  self.claim self.interactions_tsv, :proc do |filename|
    tsv = TSV.open('http://interactome3d.irbbarcelona.org/user_data/human/download/complete/interactions.dat', 
                   :header_hash => '',
                   :merge => true)
    tsv.to_s
  end

  #self.claim self.proteins, :proc do |filename|
  #  FileUtils.mkdir_p filename

  #  html =  Open.read("http://interactome3d.irbbarcelona.org/downloadset.php?queryid=human&release=current&path=complete")
  #  links = html.scan(%r{http://interactome3d.irbbarcelona.org/user_data/human/download/complete/proteins_\d+.tgz})
  #  links.each do |link|
  #    TmpFile.with_dir do |tmpdir|
  #      FileUtils.chdir tmpdir
  #      begin
  #        CMD.cmd('tar xvfz -', :in => Open.open(link))
  #      rescue
  #        Log.error("Error downloading #{ link }: #{$!.message}")
  #      end
  #      FileUtils.mv Dir.glob("*.pdb"), filename
  #    end
  #  end
  #  nil
  #end

  #self.claim self.interactions, :proc do |filename|
  #  FileUtils.mkdir_p filename

  #  html =  Open.read("http://interactome3d.irbbarcelona.org/downloadset.php?queryid=human&release=current&path=complete")
  #  links = html.scan(%r{http://interactome3d.irbbarcelona.org/user_data/human/download/complete/interactions_\d+.tgz})
  #  links.each do |link|
  #    TmpFile.with_dir do |tmpdir|
  #      FileUtils.chdir tmpdir
  #      begin
  #        CMD.cmd('tar xvfz -', :in => Open.open(link))
  #      rescue
  #        Log.error("Error downloading #{ link }: #{$!.message}")
  #      end
  #      FileUtils.mv Dir.glob("*.pdb"), filename
  #    end
  #  end
  #  nil
  #end
end

if __FILE__ == $0
  Interactome3d.proteins_tsv.produce
  Interactome3d.proteins.produce
  Interactome3d.interactions_tsv.produce
  Interactome3d.interactions.produce
end
