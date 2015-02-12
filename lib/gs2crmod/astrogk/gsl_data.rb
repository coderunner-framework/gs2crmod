class CodeRunner
	class Gs2::Astrogk
		module AstrogkGSLVectors
			def etrans_by_kfrom_by_kto_over_time_gsl_vector(options)
				read_transfers unless FileTest.exist? "energy_transfer.rb"
				etr_hash = eval(File.read("energy_transfer.rb"))
				etr_hash.values.map{|hash| hash[options[:kf]][options[:kt]]}.to_gslv

			end
		end
		include AstrogkGSLVectors
		def geometric_factors_gsl_tensor(options)
			#ops = options.dup; ops.delete :phi
		#ep ops; gets
			theta_vec = gsl_vector(:theta, options)
			factors = GSL::Tensor.alloc(6,theta_vec.size)
			factors[true, true] = 1.0
			factors
		end
		def correct_3d_options(options)
			eputs "Info: setting options[:gs2_coordinate_factor] to 1.0 (for slab geometry)"
		  options[:gs2_coordinate_factor] = 1.0
			#raise "Please specify options[:rho_star]" unless options[:rho_star]
			options[:rho_star_actual] = 1.0
			options[:q_actual] = 1.0
			options[:rhoc_actual] = 1.0
		end
	end
end
