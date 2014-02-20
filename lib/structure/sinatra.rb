require 'rbbt-util'
require 'rbbt/rest/web_tool'

include Sinatra::RbbtToolHelper

Rbbt.claim Rbbt.www.views.public.js.jmol.find(:lib), :proc do |dir|
  url = "http://sourceforge.net/projects/jmol/files/Jmol/Version%2013.0/Version%2013.0.13/Jmol-13.0.13-binary.tar.gz/download"
  TmpFile.with_file do |zip_file|
    Open.write(zip_file, Open.read(url, :mode => 'rb'), :mode => 'wb')
    TmpFile.with_file do |unzip_dir|
      FileUtils.mkdir_p unzip_dir unless File.exists? unzip_dir
      CMD.cmd("tar xvfz '#{zip_file}' -C '#{unzip_dir}'")
      FileUtils.mv(Dir.glob(File.join(unzip_dir, '*')).first, dir)
    end
  end
  nil
end


Rbbt.www.views.public.js.jmol.find(:lib).produce
