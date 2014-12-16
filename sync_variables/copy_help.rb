require 'pp'
namelists = eval(File.read('lib/gs2crmod/namelists.rb'))
p 'read namelists'
namelists[:diagnostics_config][:variables].each do |var, varhash|
  varhash[:help] ||= namelists[:gs2_diagnostics_knobs][:variables][var][:help] rescue ""
  varhash[:description] ||= namelists[:gs2_diagnostics_knobs][:variables][var][:description] rescue ""
end
File.open('lib/gs2crmod/namelists.rb', 'w'){|file| file.puts namelists.pretty_inspect}
