##########################################
# Calculations for AGK Code Runner Module
# (Only if they differ from GS2)
#
# This module contains any methods that
# begin with the word calculate. 
#
# These methods calculate results and 
# quantities that are not directly 
# obtainable from the GS2 output files.
#
##########################################


class CodeRunner
class Gs2::Astrogk 

	def read_transfers
		Dir.chdir(@directory) do
			raise "No ktrans file" unless FileTest.exist? "#@run_name.ktrans"
		file = File.open("#@run_name.ktrans", "r")
		wtr_hash = {}
		etr_hash = {}
		while line = file.gets
			#p line;
			arr =   line.split(/\s+/).map{|w| w.to_f}
			#p arr
			dud, t, kf, kt, w, e = arr
			next unless t 
			etr_hash[t] ||= {}
			etr_hash[t][kf] ||= {}
			etr_hash[t][kf][kt] = e

			wtr_hash[t] ||= {}
			wtr_hash[t][kf] ||= {}
			wtr_hash[t][kf][kt] = w




		end

		
		#p etr_hash

		File.open('energy_transfer.rb', 'w'){|file| file.puts etr_hash.pretty_inspect}
		File.open('w_transfer.rb', 'w'){|file| file.puts wtr_hash.pretty_inspect}

		end # Dir.chdir(@directory)

	end


end
end


