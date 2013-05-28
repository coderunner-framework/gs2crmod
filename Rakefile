# encoding: utf-8

require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "gs2crmod"
  gem.homepage = "http://gs2crmod.sourceforge.net"
  gem.license = "GSLv3"
  gem.summary = %Q{Module to allow CodeRunner to run and analyse the GS2 and AstroGK codes.}
  gem.description = %Q{GS2 is a gyrokinetic flux tube initial value turbulence code which can be used for fusion or astrophysical plasmas. CodeRunner is a framework for the automated running and analysis of large simulations. This module allows GS2 (and its sister code AstroGK) to harness the power of the CodeRunner framework.}
  gem.email = "edmundhighcock@sourceforge.net"
  gem.authors = ["Edmund Highcock", "Ferdinand van Wyk"]
	gem.extensions = "ext/extconf.rb"
	gem.files.include('ext/*.c', 'include/*.h', 'ext/*.rb')
	gem.required_ruby_version = '>= 1.9.1'
  # dependencies defined in Gemfile
end
Jeweler::RubygemsDotOrgTasks.new

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/**/test_*.rb'
  test.verbose = true
end

#require 'rcov/rcovtask'
#Rcov::RcovTask.new do |test|
  #test.libs << 'test'
  #test.pattern = 'test/**/test_*.rb'
  #test.verbose = true
  #test.rcov_opts << '--exclude "gems/*"'
#end

task :default => :test

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "gs2crmod #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end