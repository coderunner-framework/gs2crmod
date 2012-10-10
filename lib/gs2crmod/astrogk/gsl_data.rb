class CodeRunner
	class Astrogk
		module AstrogkGSLVectors
			def etrans_by_kfrom_by_kto_over_time_gsl_vector(options)
				read_transfers unless FileTest.exist? "energy_transfer.rb"
				etr_hash = eval(File.read("energy_transfer.rb"))
				etr_hash.values.map{|hash| hash[options[:kf]][options[:kt]]}.to_gslv

			end
		end
		include AstrogkGSLVectors
	end
end
