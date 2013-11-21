#########################
# GS2 CodeRunner module
# Data Analysis
#########################
#

#

class NumRu::NetCDF
	aliold :var
	def var(*args)
		if args[0].kind_of? Symbol
			args[0] = args[0].to_s
		end
		return old_var(*args)
	end
end
class CodeRunner
class Gs2 


eval(File.read(File.dirname(__FILE__) + '/gsl_tools.rb'), GLOBAL_BINDING, File.dirname(__FILE__) + '/gsl_tools.rb')

# def gsl_vector(name, options={})
# 	if options[:t_index] or options[:frame_index] or not [:Failed, :Complete].include? status
# 		return get_gsl_vector(name, options)
# 	else
# 		return cache[[:gsl_vector, name, options]] ||= get_gsl_vector(name, options)
# 	end
# end

def netcdf_file
	#if @runner.cache[:runs] and (open = @runner.cache[:runs].keys.find_all{|id| @runner.cache[:runs][id][:netcdf_file]}).size > 200
	#ep "my id", id
	if (open = @runner.run_list.keys.find_all{|id|  @runner.run_list[id].cache[:netcdf_file]}).size > 200
		open = open.sort_by{|id| @runner.run_list[id].cache[:netcdf_file_otime]}
		@runner.run_list[open[0]].ncclose
	end

	if cache[:netcdf_file] and not [:Complete, :Failed].include? @status
		ncclose
	end
	cache[:netcdf_file_otime] = Time.now.to_i
	cache[:netcdf_file] ||= NumRu::NetCDF.open(netcdf_filename)
	cache[:netcdf_file].sync
	cache[:netcdf_file]
end

def netcdf_filename
	@directory + '/' +  @run_name + '.out.nc'
end


def ncclose
	cache[:netcdf_file].close
	cache.delete(:netcdf_file)
end



module FixNormOption
#class << self
		def fix_norm_action(options)
			case options[:set_norm_option]
			when "t_over_m", "bd"
				if ["t_over_m", "bd"].include? @norm_option
					return :none
				else
					return :from_root_2
				end
			else
				#eputs "else", norm_option
				if ["t_over_m", "bd"].include? @norm_option
					#eputs "norm old"
					return :to_root_2
				else
					return :none
				end
			end
		end

		def fix_heat_flux_norm(tensor, options={})
			case fix_norm_action(options)
			when :none
				#eputs "none"
				return tensor
			when :to_root_2
				eputs "to_root_2"
				return tensor / 2.0**1.5
			when :from_root_2
				return tensor * 2.0**1.5
			end
		end

		# Return tensor normalised according to options[:set_norm_option] 
		# (which may be "t_over_m", "bd" or "with_root_2", "mtk" etc)
		# regardless of the original normalisation.
		#
		# <tt>power</tt> should be the power to which the reference thermal 
		# velocity is raised in the normalising quantity. For example, 
		# <tt>t</tt> is normalised to a / v_thr, so for times, power should
		# be set equal to -1.
		
		def fix_norm(tensor, power, options={})
			case fix_norm_action(options)
			when :none
				#eputs "none"
				return tensor
			when :to_root_2
				eputs "to_root_2"
				return tensor / 2.0**(0.5 * power)
			when :from_root_2
				return tensor * 2.0**(0.5 * power)
			end
		end

#end # class << self
end # module FixNormOption

include FixNormOption

def gsl_vector(name, options={})
  Dir.chdir(@directory) do
		options[:t_index_window] ||= @scan_index_window
		options.setup_time_window
		if [:ky, :kx].include? name.to_sym
			vec = fix_norm(
				GSL::Vector.alloc(netcdf_file.var(name.to_s).get.to_a.sort),
				-1, options
			) # ky, ky are normalised to 1 / rho_i
			if i = options[:interpolate_ + name.to_s.sub(/k/, '').to_sym]
				if name.to_sym == :ky
					s = (vec.size - 1)*i + 1
					#return vec.connect(GSL::Vector.alloc((vec.size-1)*(i-1)) * 0.0)
					return (0...s).map{|k| k.to_f * vec[1]}.to_gslv
				else
					size = vec.size
					#vec = vec.to_box_order
					raise "Hmmm, kx.size should be odd" unless size%2 == 1
					s = (size-1)/2 * i
					return (-s..s).to_a.map{|i| i.to_f * vec.to_box_order[1]}.to_gslv
					#new_vec = GSL::Vector.alloc((s-1)*i + 1)
					#new_vec *= 0.0
					#for j in 0...((s-1)/2+1)
						#new_vec[j] = vec[j]
					#end
					#for j in 0...((s-1)/2)
						#new_vec[-j-1] = vec[-j-1]
					#end
					#return new_vec.from_box_order
				end


			else
				return vec
			end
	  elsif [:theta].include? name.to_sym
			#ep options; gets
			#vec = GSL::Vector.alloc(netcdf_file.var(name.to_s).get({'start' => [options[:thetamin]||0], 'end' => [options[:thetamax]||-1]}).to_a)
			vec = GSL::Vector.alloc(netcdf_file.var(name.to_s).get.to_a)
			if gryfx? and options[:periodic]
				#vec = vec.connect([2.0*vec[-1] - vec[-2]].to_gslv)
				vec = vec.connect([-vec[0]].to_gslv)
			end
			if ith = options[:interpolate_theta]
				osize = vec.size
				newsize = (osize-1)*ith+1
				newvec = GSL::Vector.alloc(newsize)
				newvec[newsize-1] = vec[osize-1]# * ith.to_f
				for i in 0...(newsize-1)
					im = i%ith
					frac = im.to_f/ith.to_f
					#iold = (i-im)/(new_shape[-1]-1)*(shape[-1]-1)
					iold = (i-im)/ith
					newvec[i] =  (vec[iold] * (1.0-frac) + vec[iold+1] * frac)
				end
				vec = newvec
			end
			start = options[:thetamin]||0
			endv = options[:thetamax]||vec.size-1
			#ep ['options', options, 'vec.size', vec.size]
			vec = vec.subvector(start, (endv-start+1)).dup
			return vec
		elsif name.to_sym == :t
			#options.setup_time_window
			t = GSL::Vector.alloc(netcdf_file.var(name.to_s).get('start' => [options[:begin_element]], 'end' => [options[:end_element]]).to_a)
			t = t - t[0] if options[:sync_time]
			return fix_norm(t, -1, options) # t is normalised to a/v_thi
		end
		options = eval(options) if options.class == String
		if options[:saturated_time_average] or options[:sta]
			raise "Not Saturated" unless @saturation_time_index
			tmax = list(:t).keys.max
			return ((@saturation_time_index..tmax).to_a.map do |t_index|
				gsl_vector(name, options.dup.absorb({t_index: t_index, saturated_time_average: nil, sta: nil}))
			end).sum / (list(:t).values.max - list(:t)[@saturation_time_index])
		elsif options[:time_average] or options[:ta]
			tmax = list(:t).keys.max
			start_t = 2
			return ((start_t..tmax).to_a.map do |t_index|
				gsl_vector(name, options.dup.absorb({t_index: t_index, time_average: nil, ta: nil}))
			end).sum / (list(:t).values.max - list(:t)[start_t])
		end
		if method = self.class.instance_methods.find{|meth| (name + '_gsl_vector').to_sym == meth}
			options[:graphkit_name] = name
			return send(method, options)
		end
	end
	raise "GSL Vector #{name} not found"
end

module GSLVectors

	# The square of the potential summed over all wave numbers, indexed by time, normalised to (e/T)(rho_1/a). 
	
	def phi2tot_over_time_gsl_vector(options)

		Dir.chdir(@directory) do			#Necessary options: ky
			#log 'about to open netcdf file'
			#options.setup_time_window
			phis = netcdf_file.var('phi2').get('start'=>[options[:begin_element]], 'end'=>[options[:end_element]] ).to_a
			log 'about to allocate gsl vector'
			vec = GSL::Vector.alloc(phis)
			log 'finished'
			return fix_norm(vec, 1, options) 
		end 
	end
	def apar2_over_time_gsl_vector(options)

		Dir.chdir(@directory) do			#Necessary options: ky
			#log 'about to open netcdf file'
			#options.setup_time_window
			phis = netcdf_file.var('apar2').get('start'=>[options[:begin_element]], 'end'=>[options[:end_element]] ).to_a
			log 'about to allocate gsl vector'
			vec = GSL::Vector.alloc(phis)
			log 'finished'
			return fix_norm(vec, 1, options) 
		end 
	end

	def transient_es_heat_flux_amplification_over_kx_gsl_vector(options)
		options[:direction] = :kx
		transient_es_heat_flux_amplification_over_kxy_gsl_vector(options)
	end

	def transient_es_heat_flux_amplification_over_ky_gsl_vector(options)
		options[:direction] = :ky
		transient_es_heat_flux_amplification_over_kxy_gsl_vector(options)
	end
	def transient_es_heat_flux_amplification_over_kxy_gsl_vector(options)
		Dir.chdir(@directory) do			# i.e. phi2_by_ky_vs_time or phi2_by_kx_vs_time
			kxy = options[:direction].to_sym
			
# 			ep :growth_rate_at_ + kxy
			p send(:transient_es_heat_flux_amplification_at_species_at_ + kxy)
			return GSL::Vector.alloc(send(:transient_es_heat_flux_amplification_at_species_at_ + kxy)[options[:species_index]-1].values)

		end
	end

	def transient_amplification_over_kx_gsl_vector(options)
		options[:direction] = :kx
		transient_amplification_over_kxy_gsl_vector(options)
	end
	def transient_amplification_over_kxy_gsl_vector(options)
		options[:direction] = :ky
		transient_amplification_over_kxy_gsl_vector(options)
	end
	def transient_amplification_over_kxy_gsl_vector(options)
		Dir.chdir(@directory) do			# i.e. phi2_by_ky_vs_time or phi2_by_kx_vs_time
			kxy = options[:direction]
# 			ep :growth_rate_at_ + kxy
			return GSL::Vector.alloc(send(:transient_amplification_at_ + kxy).values)

		end
	end
	private :transient_amplification_over_kxy_gsl_vector

		# The growth rate of the fluctuations, calculated from the potential, indexed by time and normalised to vth_1/a.
	  # :kx or :kx_index must be specified in options
	#
		def growth_rate_by_kx_over_time_gsl_vector(options)
			options[:direction] = :kx
			growth_rate_by_kxy_over_time_gsl_vector(options)
		end

		# The growth rate of the fluctuations, calculated from the potential, indexed by time and normalised to vth_1/a.
	  # :ky or :ky_index must be specified in options
	
		def growth_rate_by_ky_over_time_gsl_vector(options)
			options[:direction] = :ky
			growth_rate_by_kxy_over_time_gsl_vector(options)
		end
		def growth_rate_by_kxy_over_time_gsl_vector(options)
			# i.e. time_dependent_gr_by_ky_vs_time or phi2_by_kx_vs_time
			
			kxy = options[:direction]

			phi = gsl_vector("phi2_by_#{kxy}_over_time", options).log / 2.0

			size = phi.size
			dphi = phi.subvector(1, size - 1) - phi.subvector(0, size-1)
			# NB dt already has norm fixed, dphi is dimensionless
			return fix_norm(dphi/gsl_vector('dt'), 0, options)
		end

# <MJL edits on 2013-09-19>
                  # The real frequency of the fluctuations, read from the .out file, indexed by time and normalised to vth_1/a.
                  # :ky_index or :kx_index must be specified in options.

                def frequency_by_kx_over_time_gsl_vector(options)
                        options[:direction] = :kx
                        frequency_by_kxy_over_time_gsl_vector(options)
                end

                def frequency_by_ky_over_time_gsl_vector(options)
                        options[:direction] = :ky
                        frequency_by_kxy_over_time_gsl_vector(options)
                end

                def frequency_by_kxy_over_time_gsl_vector(options)
                    kxy = options[:direction]
                    kxy_index = kxy + :_index
                    kxys = get_list_of(kxy)
                    desired_kxy = kxys[options[kxy_index]]
                    raise "No k found at the desired index" if desired_kxy.nil?

                    omega_reals = []
                    File.open(@run_name+".out",'r') do |fileHandle|
                      fileHandle.each_line do |fileLine|
                        if fileLine.include?('aky=')  # Only examine the lines of the .out file that contain frequency information.

                          index = fileLine.index('akx=')
                          raise "akx wasn't found where it was expected in the .out file." if index.nil?
                          akx = fileLine[(index+4)..-1].to_f

                          index = fileLine.index('aky=')
                          raise "aky wasn't found where it was expected in the .out file." if index.nil?
                          aky = fileLine[(index+4)..-1].to_f

                          index = fileLine.index('om=')
                          raise "om wasn't found where it was expected in the .out file." if index.nil?
                          omr = fileLine[(index+3)..-1].to_f
                          if kxy == :kx
                            # You need to be careful when testing equality of the desired k with the k in the .out file
                            # since the .out file is only written to ~ 5 significant digits:
                            omega_reals << omr if ((desired_kxy - akx).abs/(desired_kxy.abs + 1e-7) < 1e-4)
                          else
                            omega_reals << omr if ((desired_kxy - aky).abs/(desired_kxy.abs + 1e-7) < 1e-4)
                          end
                        end
                      end
                     end
                     raise "No real frequencies found in the .out file for the desired k" if (omega_reals.size==0)
                     GSL::Vector.alloc(omega_reals)
                end
# </MJL>


		# The size of each time step,  indexed by time, normalised to a/v_th1.
		
	 	def dt_gsl_vector(options)
			t = gsl_vector('t', options)
			size = t.size
			# NB t already has norm fixed
		  return t.subvector(1, size - 1) - t.subvector(0, size-1)
		end

		# The growth rate, calculated from the potential, indexed by kx. Only makes sense in linear calculations. 
	def growth_rate_over_kx_gsl_vector(options)
		options[:direction] = :kx
		growth_rate_over_kxy_gsl_vector(options)
	end
		# The growth rate, calculated from the potential, indexed by ky. Only makes sense in linear calculations. 
	def growth_rate_over_ky_gsl_vector(options)
		options[:direction] = :ky
		growth_rate_over_kxy_gsl_vector(options)
	end

	def growth_rate_over_kxy_gsl_vector(options)
		Dir.chdir(@directory) do			# i.e. phi2_by_ky_vs_time or phi2_by_kx_vs_time
			kxy = options[:direction]
# 			ep :growth_rate_at_ + kxy
			return GSL::Vector.alloc(send(:growth_rate_at_ + kxy).values)

		end
	end
	private :growth_rate_over_kxy_gsl_vector

	# The growth rate, calculated from the potential, indexed by kx. Only makes sense in linear calculations. 
	def growth_rate_over_kx_slice_gsl_vector(options)
		Dir.chdir(@directory) do
			slice_of_growth_rates = send(:growth_rate_at_ky_at_kx)[options[:ky]].values
			raise "Something went wrong: slice of growth rates seems empty" if slice_of_growth_rates.nil?
			return GSL::Vector.alloc(slice_of_growth_rates)
			#return GSL::Vector.alloc(send(:growth_rate_at_ky_at_kx[ky]).values)
		end
	end

	# The growth rate, calculated from the potential, indexed by ky. Only makes sense in linear calculations. 
	def growth_rate_over_ky_slice_gsl_vector(options)
		Dir.chdir(@directory) do
			slice_of_growth_rates = send(:growth_rate_at_ky_at_kx).values.map{|h| h[options[:kx]]}
			raise "Something went wrong: slice of growth rates seems empty" if slice_of_growth_rates.nil?
			return GSL::Vector.alloc(slice_of_growth_rates)
		end
	end

	# Frequency, indexed over ky, taken direct from the gs2 output file
	def frequency_over_ky_gsl_vector(options)
		  options.convert_to_index(self, :kx)
			return GSL::Vector.alloc(gsl_vector('ky').to_a.map{|ky| frequency_at_ky_at_kx[ky].values[options[:kx_index]-1]})
	end

	def es_heat_by_kx_over_time_gsl_vector(options)
		options[:direction] = :kx
		es_heat_by_kxy_over_time_gsl_vector(options)
	end
	def es_heat_by_ky_over_time_gsl_vector(options)
		options[:direction] = :ky
		es_heat_by_kxy_over_time_gsl_vector(options)
	end

	def es_heat_by_kxy_over_time_gsl_vector(options)
		Dir.chdir(@directory) do
			kxy = options[:direction]
			kxy_index = kxy + :_index
			options.convert_to_index(self, kxy)
			if kxy==:ky
				lkx = list(:kx)
				es_heat_av = (lkx.keys.map do |kx_index|		
					es_heat =  netcdf_file.var('es_heat_by_k').get({'start' => [kx_index-1,options[:ky_index]-1,options[:species_index]-1, 0], 'end' => [kx_index-1,options[:ky_index]-1,options[:species_index]-1, -1]})
					#ep phi.shape
					es_heat.reshape(*es_heat.shape.values_at(3))
				end).sum / lkx.size
				return es_heat_av.to_gslv
			else
				lky = list(:ky)
				es_heat_av = (lky.keys.map do |ky_index|		
					es_heat =  netcdf_file.var('es_heat_by_k').get({'start' => [options[:kx_index]-1,ky_index-1,options[:species_index]-1, 0], 'end' => [options[:kx_index]-1,ky_index-1,options[:species_index]-1, -1]})
					#ep phi.shape
					es_heat.reshape(*es_heat.shape.values_at(3))
				end).sum / lky.size
				return es_heat_av.to_gslv
			end

		end
	end
	private :es_heat_by_kxy_over_time_gsl_vector

	def phi2_by_kx_over_time_gsl_vector(options)
		options[:direction] = :kx
	 	phi2_by_kxy_over_time_gsl_vector(options)
	end
	def phi2_by_ky_over_time_gsl_vector(options)
		options[:direction] = :ky
	 	phi2_by_kxy_over_time_gsl_vector(options)
	end
	def phi2_by_kxy_over_time_gsl_vector(options)
		Dir.chdir(@directory) do	
			# i.e. phi2_by_ky_vs_time or phi2_by_kx_vs_time
			
			kxy = options[:direction]
			if list(kxy).size == 1
				return phi2tot_over_time_gsl_vector(options)
			end
			kxy_index = kxy + :_index
			
			
			#Necessary options: :ky or :kx
			#Optional options: :t_index_window
			# 		eputs "got here"
			#options[:begin_element], options[:end_element] = (options[:t_index_window] ? options[:t_index_window].map{|ind| ind -1} : [0, -1])
			phi_t_array=nil
			if @grid_option == "single"
				phi_t_array = netcdf_file.var('phi2').get('start' => [options[:begin_element]], 'end' => [options[:end_element]]).to_a.flatten
			else
# 				value = options[:ky]
# 				eputs value
# 				get_list_of(:ky)
# 				index = @ky_list.find{|index,val| (val-value).abs < Float::EPSILON}[0]
#         ep options
				options.convert_to_index(self, kxy)
 				#ep options
				phi_t_array = netcdf_file.var("phi2_by_#{kxy}").get('start' => [options[kxy_index] - 1, options[:begin_element]], 'end' => [options[kxy_index] - 1, options[:end_element]]).to_a.flatten
# 				eputs 'phi_t_array.size', phi_t_array.size
			end
			return GSL::Vector.alloc(phi_t_array)

		end
	end
	private :phi2_by_kxy_over_time_gsl_vector


	def phi2_by_mode_over_time_gsl_vector(options)
		Dir.chdir(@directory) do			#Necessary options: :ky and :kx
			#Optional options: :t_index_window
			# 		eputs "got here"
			#options[:begin_element], options[:end_element] = (options[:t_index_window] ? options[:t_index_window].map{|ind| ind -1} : [0, -1])
			options.setup_time_window			
			phi_t_array=nil
			if @grid_option == "single"
				phi_t_array = netcdf_file.var('phi2').get('start' => [options[:begin_element]], 'end' => [options[:end_element]]).to_a.flatten
			else
# 				value = options[:ky]
# 				eputs value
# 				get_list_of(:ky)
# 				index = @ky_list.find{|index,val| (val-value).abs < Float::EPSILON}[0]
				options.convert_to_index(self, :kx, :ky)
# 				p options
				phi_t_array = netcdf_file.var("phi2_by_mode").get('start' => [options[:kx_index] - 1, options[:ky_index] - 1, options[:begin_element]], 'end' => [options[:kx_index] - 1, options[:ky_index] - 1, options[:end_element]]).to_a.flatten
# 				eputs 'phi_t_array.size', phi_t_array.size
			end
			return GSL::Vector.alloc(phi_t_array)

		end
	end

	def tpar2_by_mode_over_time_gsl_vector(options)
		Dir.chdir(@directory) do			#Necessary options: :ky and :kx
			#Optional options: :t_index_window
			# 		eputs "got here"
			#options[:begin_element], options[:end_element] = (options[:t_index_window] ? options[:t_index_window].map{|ind| ind -1} : [0, -1])
			options.setup_time_window			
			tpar_t_array=nil
			if @grid_option == "single"
				tpar_t_array = netcdf_file.var('tpar2').get('start' => [options[:begin_element]], 'end' => [options[:end_element]]).to_a.flatten
			else
# 				value = options[:ky]
# 				eputs value
# 				get_list_of(:ky)
# 				index = @ky_list.find{|index,val| (val-value).abs < Float::EPSILON}[0]
				options.convert_to_index(self, :kx, :ky, :species)
# 				p options
				tpar_t_array = netcdf_file.var("tpar2_by_mode").get('start' => [options[:kx_index] - 1, options[:ky_index] - 1, options[:species_index] - 1, options[:begin_element]], 'end' => [options[:kx_index] - 1, options[:ky_index] - 1, options[:species_index] - 1, options[:end_element]]).to_a.flatten
# 				eputs 'tpar_t_array.size', tpar_t_array.size
			end
			return GSL::Vector.alloc(tpar_t_array)

		end
	end

	def tperp2_by_mode_over_time_gsl_vector(options)
		Dir.chdir(@directory) do			#Necessary options: :ky and :kx
			#Optional options: :t_index_window
			# 		eputs "got here"
			#options[:begin_element], options[:end_element] = (options[:t_index_window] ? options[:t_index_window].map{|ind| ind -1} : [0, -1])
			options.setup_time_window			
			tperp_t_array=nil
			if @grid_option == "single"
				tperp_t_array = netcdf_file.var('tperp2').get('start' => [options[:begin_element]], 'end' => [options[:end_element]]).to_a.flatten
			else
# 				value = options[:ky]
# 				eputs value
# 				get_list_of(:ky)
# 				index = @ky_list.find{|index,val| (val-value).abs < Float::EPSILON}[0]
				options.convert_to_index(self, :kx, :ky, :species)
# 				p options
				tperp_t_array = netcdf_file.var("tperp2_by_mode").get('start' => [options[:kx_index] - 1, options[:ky_index] - 1, options[:species_index] - 1, options[:begin_element]], 'end' => [options[:kx_index] - 1, options[:ky_index] - 1, options[:species_index] - 1, options[:end_element]]).to_a.flatten
# 				eputs 'tperp_t_array.size', tperp_t_array.size
			end
			return GSL::Vector.alloc(tperp_t_array)

		end
	end

	def phi0_by_kx_by_ky_over_time_gsl_vector(options)
		Dir.chdir(@directory) do		
			options.convert_to_index(self, :kx, :ky)
			phi0_array = netcdf_file.var('phi0').get.to_a.map{|arr| arr[options[:kx_index] - 1][options[:ky_index] - 1][options[:ri]]}
			return GSL::Vector.alloc(phi0_array)

		end
	end

	def linked_kx_elements_gsl_vector(options)
		Dir.chdir(@directory) do		
			return GSL::Vector.alloc([0]) if @grid_option == "single" or agk?
			if agk? or (@s_hat_input or @shat).abs < 1.0e-5
				#p 'op1', options
				options.convert_to_index(self, :ky, :kx)
				#p 'op2', options
				#eputs "No Magnetic Shear"

# 				begin
# 					options.convert_to_index(:kx)
# 				rescue
# 					raise "Must specify kx or kx_index if no magnetics shear"
# 				end
# # 				theta0 = (options[:theta0] || 0)
# # 				theta0 += jump(options) if @g_exb

				#theta0 = (options[:kx_index])
				#if @g_exb and @g_exb.abs > 0.0
          #theta0 += jump(options)
          #theta0 = theta0%((list(:kx).size-1)/2) if list(:kx).size > 1
        #end

				return GSL::Vector.alloc([options[:kx_index] - 1])
			end
			
			options.convert_to_index(self, :ky, :kx)
			nkx = netcdf_file.var('kx').dims[0].length
# 			p nkx
			stride = @jtwist * (options[:ky_index] -1 )
			#stride = 3
			nlinks = [(nkx / stride).floor, 1].max 
			theta0 = options[:kx_index] % @jtwist  #(options[:theta0] || 0)
			log 'stride', stride, 'nlinks', nlinks, 'theta0', theta0
			#if @g_exb and @jtwist > 1 #and options[:t_index]
# 				kx_shift = list(:ky)[options[:ky_index]]  * @g_exb
# 				p list(:t)[options[:t_index]], options[:t_index], kx_shift
				
# 				kx_shift *=  list(:t)[(options[:t_index] or list(:t).keys.max)]
# 				jump = (kx_shift / list(:kx)[2]).round
				#theta0  += (@jtwist - jump(options) % @jtwist) % @jtwist
				
# 				else
# 					jump = 0	
			#end
			ep 'stride', stride, 'nlinks', nlinks, 'theta0', theta0
			ep GSL::Vector.indgen(nlinks / 2,  nkx + theta0 - nlinks / 2 * stride, stride).connect(GSL::Vector.indgen(nlinks / 2, theta0, stride)).reverse if nlinks > 1
			#return [7,5,3,1,34].to_gslv
			return GSL::Vector.alloc([theta0 % jtwist]) if nlinks ==1
			return GSL::Vector.indgen(nlinks / 2,  nkx + theta0 - nlinks / 2 * stride, stride).connect(GSL::Vector.indgen(nlinks / 2, theta0, stride)).reverse

		end
	end

		def spectrum_over_kpar_gsl_vector(options)
		Dir.chdir(@directory) do
 # , /kpar_spectrum/
			#ep 'zero?', (@s_hat_input||@shat)==0.0
			unless agk? or (@s_hat_input||@shat||0.0).abs<1.0e-5
				phi = gsl_vector_complex('phi_along_field_line', options)
        phi = phi.subvector(0,phi.size-1)
				#i = 0
				#phi = phi.collect{|re,im| 
					#i+=1; GSL::Complex.alloc(Math.sin(0.1*i), Math.cos(0.1*i))+ 
					#GSL::Complex.alloc(Math.sin(0.4*i), Math.cos(0.4*i)) 

				#}
				##GraphKit.quick_create([phi.square]).gnuplot
				phi_k = phi.forward
				phi_kr = phi_k.square
				case phi_kr.size%2
				when 0
					spec = phi_kr.subvector((phi_kr.size+2)/2, (phi_kr.size-2)/2).connect(phi_kr.subvector(0, (phi_kr.size+2)/2))
				when 1
					spec = phi_kr.subvector((phi_kr.size + 1)/2, (phi_kr.size-1)/2).connect(phi_kr.subvector(0, (phi_kr.size+1)/2))
				end
				##spec = phi_kr
				#ep 'spec.class', spec.class
				return spec
			else

				gm = gsl_matrix('spectrum_over_ky_over_kpar', options)
				vec = GSL::Vector.alloc(gm.shape[1])
				vec.set_all(0.0)
				for ky_element in 0...gm.shape[0]
					vec+= gm.row(ky_element)
				end
				vec = vec/gm.shape[0]
				return vec
			end
		end
	end
	def kpar_gsl_vector(options)

		Dir.chdir(@directory) do	
			if agk? or (@s_hat_input||@shat).abs	< 1.0e-5
				dk = 1
				phi = list(:theta).values
			else
				kxe = gsl_vector('linked_kx_elements', options)
				dk = 1.0/kxe.size
			phi = gsl_vector_complex('phi_along_field_line', options)
			end
			case phi.size%2
			when 0
				kpar = GSL::Vector.indgen(phi.size-1, -((phi.size-3)/2))*dk
			when 1
				kpar = GSL::Vector.indgen(phi.size-1, -((phi.size-2)/2))*dk
			end
      #ep 'kpar', kpar, 'phi.size', phi.size

			#ep 'kpar.class', kpar.class
			return kpar

		end
	end

	def phi_along_field_line_gsl_vector(options)
		Dir.chdir(@directory) do	
			complex_phi_vector= gsl_vector_complex('phi_along_field_line', options)
			case options[:imrc]
			when :im
				phi_vector = complex_phi_vector.imag
			when :mag
				mag = true
				phi_vector = complex_phi_vector.abs2
			when :corr
				thetas = gsl_vector('theta_along_field_line', options)
				min = thetas.abs.to_a.index(thetas.abs.min)
				at_0 = complex_phi_vector[min]
# 				ep at_0.class
				phi_vector = (complex_phi_vector * (at_0 / at_0.mag).conj).real
# 				gsl_complex('correcting_phase', options)).real
			when :real
				phi_vector = complex_phi_vector.real
			else
				raise CRError.new("options[:imrc] was: #{options[:irmc]}")
			end
			phi_vector *= -1.0 if options[:flip]
			(phi_vector /= phi_vector.abs.max; phi_vector *= (options[:height] || 1.0)) if options[:norm]
			phi_vector = phi_vector.reverse if options[:rev]
			return phi_vector

		end
	end

	def theta_along_field_line_gsl_vector(options)
		Dir.chdir(@directory) do
			case @grid_option
			when "single", "range"
				theta_vector = gsl_vector(:theta)
			when "box"
				#eputs "Start theta_along_field_line"

				kx_elements = gsl_vector('linked_kx_elements', options).to_a
				#if @grid_option == "range"
					#kx_elements = kx_elements.to_gslv.from_box_order.to_a
				#end
				ep 'kx_elements', kx_elements.to_a
# 				ep list(:kx).keys.max
# 				ep kx_elements[0], list(:kx)[(kx_elements[0] + 1).to_i]
# 				ep kx_elements[-1], list(:kx)[(kx_elements[-1] + 1).to_i]
				thetas = gsl_vector(:theta)
# 				ep thetas
				#eputs "End theta_along_field_line"
				return thetas if agk? or (@s_hat_input or @shat).abs < 1.0e-5
				if gryfx?
					theta_list = ((1..kx_elements.size).to_a.map do |i|
						thetas * i
					end)
					thetas = theta_list.inject{|o,n| o.connect(n)}
					thetas -= Math::PI*(kx_elements.size-1)
					return thetas

				end
				theta_list = (kx_elements.map do |element|
				        
					kx = list(:kx)[(element + 1).to_i]
# 				             ep element
 				             #ep 'kx', kx, 'shat', (@s_hat_input or @shat), 'ky',   list(:ky)[options[:ky_index]]
					thetas - 1.0 / (@s_hat_input or @shat) / list(:ky)[options[:ky_index]] * kx
				end).inject{|old, new| old.connect(new)}
# 				thetas = gsl_vector(:theta) - 1.0 / @shat / list(:ky)[options[:ky_index]] * list(:kx)[(kx_elements[0] + 1).to_i] #- Math::PI*(kx_elements.size  - 1)
# 				get_list_of(:ky, :t)
# 				if @g_exb #and options[:t_index]

					if options[:moving]
						theta_list = theta_list  -  Math::PI * 2.0 * (jump(options) / @jtwist)
					else
# 						ep 'jump % jtwist is!!', jump(options) % @jtwist
						theta_list = theta_list - Math::PI * 2.0 / @nx.to_f * ((jump(options) % @jtwist).to_f / @jtwist.to_f)
					end
# 					jump = 0	
# 				end
# 				theta_list = thetas.dup #gsl_vector(:theta) - Math::PI*kx_elements.size
# 				(kx_elements.size - 1).times do 
# 					thetas = thetas + Math::PI * 2.0
# 					theta_list = theta_list.connect(thetas)
# 				end
# 				pp theta_list.to_a.values_at(0, theta_list.size - 1)
# 				pp theta_list.to_a.max
				theta_vector = theta_list
			end
# 			theta_vector = theta_vector.reverse if options[:rev]
			theta_vector *= (@shat) if options[:z]
			return theta_vector

		end
	end

	def phi_for_eab_movie_gsl_vector(options)
		Dir.chdir(@directory) do			#options required are x_index, y_index and tm_index (Time)
			mvf_name = @run_name + '.movie.nc'
			raise CRError.new("cannot find file #{mvf_name}") unless FileTest.exist? mvf_name
			ncf = NumRu::NetCDF.open(mvf_name)
# 			p ncf.var('phi_by_xmode').get.to_a[0][0][0]
			return GSL::Vector.alloc(ncf.var('phi_by_xmode').get.to_a[options[:tm_index] - 1].map{|xy_arr| xy_arr[options[:x_index] - 1][options[:y_index] - 1]})

		end
	end

		def hflux_tot_over_time_gsl_vector(options)
			Dir.chdir(@directory) do
				options.setup_time_window
				narr = netcdf_file.var('hflux_tot').get('start' => [options[:begin_element]], 'end' => [options[:end_element]])
	 			#eputs 'Got narr'
				#ep 'hflux_tot', hflux
				#eputs "fixing norm"
				return fix_heat_flux_norm(GSL::Vector.alloc(narr.to_a), options)
			end
		end
		alias :hflux_tot_gsl_vector :hflux_tot_over_time_gsl_vector
		def es_heat_flux_over_time_gsl_vector(options)
			Dir.chdir(@directory) do

				options.setup_time_window
				return GSL::Vector.alloc(netcdf_file.var('es_heat_flux').get('start' => [options[:species_index].to_i - 1, options[:begin_element]], 'end' => [options[:species_index].to_i - 1, options[:end_element]]).to_a.flatten)
			end
		end
		def es_heat_par_over_time_gsl_vector(options)
			Dir.chdir(@directory) do

				options.setup_time_window
				return GSL::Vector.alloc(netcdf_file.var('es_heat_par').get('start' => [options[:species_index].to_i - 1, options[:begin_element]], 'end' => [options[:species_index].to_i - 1, options[:end_element]]).to_a.flatten)
			end
		end
		alias :es_heat_par_gsl_vector :es_heat_par_over_time_gsl_vector
		def es_heat_perp_over_time_gsl_vector(options)
			Dir.chdir(@directory) do

				options.setup_time_window
				return GSL::Vector.alloc(netcdf_file.var('es_heat_perp').get('start' => [options[:species_index].to_i - 1, options[:begin_element]], 'end' => [options[:species_index].to_i - 1, options[:end_element]]).to_a.flatten)
			end
		end
		alias :es_heat_perp_gsl_vector :es_heat_perp_over_time_gsl_vector
		def es_heat_flux_over_time_gsl_vector(options)
			Dir.chdir(@directory) do

				options.setup_time_window
				return GSL::Vector.alloc(netcdf_file.var('es_heat_flux').get('start' => [options[:species_index].to_i - 1, options[:begin_element]], 'end' => [options[:species_index].to_i - 1, options[:end_element]]).to_a.flatten)
			end
		end
		def es_mom_flux_over_time_gsl_vector(options)
			Dir.chdir(@directory) do

				options.setup_time_window
				return GSL::Vector.alloc(netcdf_file.var('es_mom_flux').get('start' => [options[:species_index].to_i - 1, options[:begin_element]], 'end' => [options[:species_index].to_i - 1, options[:end_element]]).to_a.flatten)
			end
		end
		# Velocity space diagnostics: fraction of dist func in higher 
		# pitch angle harmonics
		def lpc_pitch_angle_gsl_vector(options)
			raise "Velocity space lpc diagnostics not found" unless FileTest.exist? "#@directory/#@run_name.lpc"
			lpc = GSL::Vector.filescan("#@directory/#@run_name.lpc")
			return lpc[1]
		end
		# Velocity space diagnostics: fraction of dist func in higher 
		# energy harmonics
		def lpc_energy_gsl_vector(options)
			raise "Velocity space lpc diagnostics not found" unless FileTest.exist? "#@directory/#@run_name.lpc"
			lpc = GSL::Vector.filescan("#@directory/#@run_name.lpc")
			return lpc[2]
		end
		# Velocity space diagnostics: integral error due to 
		# pitch angle resolution
		def vres_pitch_angle_gsl_vector(options)
			raise "Velocity space vres diagnostics not found" unless FileTest.exist? "#@directory/#@run_name.vres"
			vres = GSL::Vector.filescan("#@directory/#@run_name.vres")
			return vres[1]
		end
		# Velocity space diagnostics: integral error due to 
		# energy resolution
		def vres_energy_gsl_vector(options)
			raise "Velocity space vres diagnostics not found" unless FileTest.exist? "#@directory/#@run_name.vres"
			vres = GSL::Vector.filescan("#@directory/#@run_name.vres")
			return vres[2]
		end
		def par_mom_flux_over_time_gsl_vector(options)
		Dir.chdir(@directory) do

			options.setup_time_window
			# This is a hack... one day some one will put it in the NetCDF file (haha).
			momlines = `grep parmom #@run_name.out`
			mom = []
			momlines.scan(Regexp.new("#{LongRegexen::FLOAT.to_s}$")) do 
				mom.push $~[:float].to_f
			end
			options[:end_element] = (mom.size + options[:end_element]) if options[:end_element] < 0
# 			p options
			return GSL::Vector.alloc(mom).subvector(options[:begin_element], options[:end_element] - options[:begin_element] + 1)
		end
		end

		def perp_mom_flux_over_time_gsl_vector(options)

			Dir.chdir(@directory) do
				options.setup_time_window
				# This is a hack... one day some one will put it in the NetCDF file (haha).
				momlines = `grep perpmom #@run_name.out`
				mom = []
				momlines.scan(Regexp.new("#{LongRegexen::FLOAT.to_s}$")) do 
					mom.push $~[:float].to_f
				end
				options[:end_element] = (mom.size + options[:end_element]) if options[:end_element] < 0
	# 			p options
				return GSL::Vector.alloc(mom).subvector(options[:begin_element], options[:end_element] - options[:begin_element] + 1)
			end
		end

		def scan_parameter_value_gsl_vector(options)
			return GSL::Vector.alloc(netcdf_file.var('scan_parameter_value').get.to_a)
		end
		def spectrum_over_kx_gsl_vector(options)
			options[:direction] = :kx
			spectrum_over_kxy_gsl_vector(options)
		end
		def spectrum_over_ky_gsl_vector(options)
			options[:direction] = :ky
			spectrum_over_kxy_gsl_vector(options)
		end
		def spectrum_over_kxy_gsl_vector(options)
			Dir.chdir(@directory) do
				# i.e. spectrum_over_ky or spectrum_over_kx
				kxy = options[:direction]
	# 			eputs options[:t_index]
				raise "Spectrum makes no sense for single modes" if @grid_option == "single"

				options.convert_to_index(:t) if options[:t] or options[:t_element]
	# 			eputs options[:t_index]

				options[:t_index] ||= list(:t).keys.max
	# 			eputs options[:t_index]
				phi_array = netcdf_file.var("phi2_by_#{kxy}").get('start' => [0, options[:t_index] - 1], 'end' => [-1, options[:t_index] - 1]).to_a.flatten
				v = GSL::Vector.alloc(phi_array)
				v = v.from_box_order if kxy == :kx
				v = v.mul(gsl_vector(kxy).square) unless options[:phi2_only]
				return v
			end
		end
		def x_gsl_vector(options)
			raise "options nakx and interpolate_x are incompatible" if options[:nakx] and options[:interpolate_x]
			kx = gsl_vector(:kx, options)
			lx = 2*Math::PI/kx.to_box_order[1]
			#ep 'lx', lx
			nx = options[:nakx]||kx.size
			GSL::Vector.indgen(nx, 0, lx/nx)
		end
		def y_gsl_vector(options)
			raise "options naky and interpolate_y are incompatible" if options[:naky] and options[:interpolate_y]
			ky = gsl_vector(:ky, options)
			ly = 2*Math::PI/ky[1]
			ny = options[:naky]||ky.size
			ysize = ny*2-2+ny%2
			GSL::Vector.indgen(ysize, 0, ly/ysize)
		end
		def zonal_spectrum_gsl_vector(options)
			Dir.chdir(@directory) do

		     gmzf = gsl_matrix('spectrum_over_ky_over_kx',options)
                     veczf = GSL::Vector.alloc(gmzf.shape[1])
                   # p gmzf.get_row(0).size
                   # p gmzf.get_row(0)
		     gmzf.shape[1].times{|i| veczf[i] = gmzf[0,i]}
		     return veczf
		#else
			#raise CRError.new("Unknown gsl_vector requested: #{name}")
			end
	#  			eputs data; gets
		end
  
end	# module GSLVectors	
include GSLVectors

def gsl_vector_complex(name, options={})
	options = eval(options) if options.class == String

		if method = self.class.instance_methods.find{|meth| (name + '_gsl_vector_complex').to_sym == meth}
			options[:graphkit_name] = name
			return send(method, options)
		end
end

module GSLVectorComplexes

	def phi_along_field_line_gsl_vector_complex(options)
	Dir.chdir(@directory) do
# 		eputs options[:ky]
# 		eputs Dir.pwd
			#eputs "Start phi_along_field_line"
			options.convert_to_index(self, :ky)
			if options[:t_index] or options[:t]
				#extra option required is t_index
				raise CRFatal.new("write_phi_over_time is not enabled so this function won't work") unless @write_phi_over_time
				
				options.convert_to_index(self, :t)
				case @grid_option
				when "single"
					temp = GSL::Vector.alloc(netcdf_file.var('phi_t').get({'start' => [0,0,0,0, options[:t_index] - 1], 'end' => [-1,-1,0,0, options[:t_index] - 1]}).to_a[0][0][0].flatten)
				when "range"
					a = netcdf_file.var('phi_t').get({'start' => [0, 0, options[:kx_index]-1, options[:ky_index] - 1, options[:t_index] - 1], 'end' => [-1, -1, options[:kx_index]-1, options[:ky_index] - 1, options[:t_index]-1]})
					#temp =  GSL::Vector.alloc(a.to_a[0].values_at(*kx_elements).flatten)
					temp =  GSL::Vector.alloc(a.to_a[0][0].flatten)
				when "box"
					options.convert_to_index(self, :ky, :kx)
					kx_elements = gsl_vector('linked_kx_elements', options).to_a
	#  				pp kx_elements
					a = netcdf_file.var('phi_t').get({
						'start' => [0,0,0,options[:ky_index] - 1, options[:t_index] - 1], 
						'end' => [-1,-1,-1, options[:ky_index] - 1, options[:t_index] - 1]
					}).to_a[0][0].values_at(*kx_elements).flatten
	# 				pp a.index(nil)
	# 				temp = GSL::Vector.alloc(netcdf_file.var('phi').get.to_a[options[:ky_index] - 1 ].values_at(*kx_elements).flatten)
					#ep a
					temp = GSL::Vector.alloc(a)
				end

				#eputs "End phi_along_field_line"
				return GSL::Vector::Complex.alloc(temp.subvector_with_stride(0, 2), temp.subvector_with_stride(1, 2)) 
			else
				case @grid_option
				when "single"
					temp = GSL::Vector.alloc(netcdf_file.var('phi').get({'start' => [0,0, 0, 0], 'end' => [-1,-1,0,0]}).to_a.flatten)
				when "range"
					a = netcdf_file.var('phi').get({'start' => [0, 0, 0, options[:ky_index] - 1], 'end' => [-1, -1, -1, options[:ky_index] - 1]})
					#temp =  GSL::Vector.alloc(a.to_a[0].values_at(*kx_elements).flatten)
					temp =  GSL::Vector.alloc(a.to_a[0][0].flatten)
				when "box"
					ep 'kx_elements', kx_elements = gsl_vector('linked_kx_elements', options).to_a
					a = netcdf_file.var('phi').get({'start' => [0, 0, 0, options[:ky_index] - 1], 'end' => [-1, -1, -1, options[:ky_index] - 1]})
					temp =  GSL::Vector.alloc(a.to_a[0].values_at(*kx_elements).flatten)
				else
					raise "invalid grid option"
				end
				
				vector = GSL::Vector::Complex.alloc(temp.subvector_with_stride(0, 2), temp.subvector_with_stride(1, 2))
				#ep 'vector', vector.real
				return vector
			end
		end
	#  			eputs data; gets
	end

end		
include GSLVectorComplexes

def gsl_matrix(name, options={})
	options = eval(options) if options.class == String
	if options[:saturated_time_average] or options[:sta]
		raise "Not Saturated" unless @saturation_time_index
		tmax = list(:t).keys.max
		return ((@saturation_time_index..tmax).to_a.map do |t_index|
			gsl_matrix(name, options.dup.absorb({t_index: t_index, saturated_time_average: nil, sta: nil}))
		end).sum / (list(:t).values.max - list(:t)[@saturation_time_index])
	end
	if method = self.class.instance_methods.find{|meth| (name + '_gsl_matrix').to_sym == meth}
			options[:graphkit_name] = name
			return send(method, options)
	end
end

module GSLMatrices
	def growth_rate_over_ky_over_kx_gsl_matrix(options)
		if @growth_rate_at_ky_at_kx.nil?
		   raise("The CodeRunner variable growth_rate_at_ky_at_kx does not seem to have been calculated for this run. This may result when the environment variable GS2_CALCULATE_ALL is not set when the run was analyzed. Try setting GS2_CALCULATE_ALL and then re-analyze the run using, e.g. from the command line,\n $ coderunner rc 'cgrf\' -j #{@id}")
      		end
		array = @growth_rate_at_ky_at_kx.values.map{|h| h.values}
		return GSL::Matrix.alloc(array.flatten, array.size, array[0].size)
	end
	def transient_amplification_over_ky_over_kx_gsl_matrix(options)
			array = @transient_amplification_at_ky_at_kx.values.map{|h| h.values}
			return GSL::Matrix.alloc(array.flatten, array.size, array[0].size)
	end
	def es_heat_flux_over_ky_over_kx_gsl_matrix(options)
	Dir.chdir(@directory) do
			raise "Heat flux spectrum makes no sense for single modes" if @grid_option == "single"
			options.convert_to_index(:t) if options[:t] or options[:t_element]
			options[:t_index] ||= list(:t).keys.max
			#es_heat_by_k index order (in Fortran) is kx, ky, t
			es_heat_narray = netcdf_file.var("es_heat_by_k").get('start' => [0, 0, 0, options[:t_index] - 1], 'end' => [-1, -1, 0, options[:t_index] - 1])
			es_heat_narray.reshape!(*es_heat_narray.shape.slice(0..1))
			
			gm =  es_heat_narray.to_gm.move_cols_from_box_order
			if options[:limit]
				for i in 0...gm.shape[0]
					for j in 0...gm.shape[1]
# 						j+= extra if 
						gm[i, j] = [[gm[i,j], (options[:limit][0] or gm[i,j])].max, (options[:limit][1] or gm[i,j])].min
# 						mat[i, j+extra] = gm[i,-j] unless j==0
					end
				end
			end
			return gm
	end
	end
	def spectrum_over_ky_over_kx_gsl_matrix(options)
	Dir.chdir(@directory) do
			raise "Spectrum makes no sense for single modes" if @grid_option == "single"
			options.convert_to_index(:t) if options[:t] or options[:t_element] 
			options[:t_index] ||= list(:t).keys.max
			#phi2_by_mode index order (in Fortran) is kx, ky, t
			phi_narray = netcdf_file.var("phi2_by_mode").get('start' => [0, 0, options[:t_index] - 1], 'end' => [-1, -1, options[:t_index] - 1])
			phi_narray.reshape!(*phi_narray.shape.slice(0..1))
			
			gm =  phi_narray.to_gm.move_cols_from_box_order
			if options[:times_kx4] or options[:times_kx2]
# 				puts 'normalising'
				vals = list(:kx).values.sort
				for i in 0...gm.shape[0]
					for j in 0...gm.shape[1]
# 						p vals[j]
						gm[i,j] =  gm[i,j] * (vals[j])**4 if options[:times_kx4]
						gm[i,j] =  gm[i,j] * (vals[j])**2 if options[:times_kx2]
					end
				end
			end
						if options[:no_zonal]
				
				for i in 0...gm.shape[1]
					gm[0,i] = 0.0
				end
			end
			if options[:log]
				gm = gm.log
			end

			return gm
	end
	end
	def spectrum_over_ky_over_kpar_gsl_matrix(options)
	Dir.chdir(@directory) do

			#:re, :theta, :kx, :ky
			lkx = list(:kx)
		
			if options[:t_index] or options[:t]
				#extra option required is t_index
				raise CRFatal.new("write_phi_over_time is not enabled so this function won't work") unless @write_phi_over_time
				options.convert_to_index(self, :t)
			end
			temp = phi_av = (lkx.keys.map do |kx_index|		
				if options[:t_index]
					phi =  netcdf_file.var('phi_t').get({'start' => [0,0,kx_index-1,0, options[:t_index] - 1], 'end' => [-1,-2,kx_index-1,-1, options[:t_index] - 1]})
				else
					phi = netcdf_file.var('phi').get({'start' => [0, 0, kx_index - 1, 0], 'end' => [-1, -2, kx_index-1, -1]})
				end
				#ep phi.shape
				phi.reshape(*phi.shape.values_at(0,1,3))
			end).sum / lkx.size

			phi_t = phi_av.to_a #.map{|arr| arr.transpose}.transpose.map{|a| a.transpose}
			#ep 'phi_t', phi_t.size, phi_t[0].size, phi_t[0][0].size
			gvky = gsl_vector('ky')
			gm = GSL::Matrix.alloc(gvky.size, gsl_vector('theta').size-1)
			for ky_element in 0...gm.shape[0]
				#p phi_t[ky_element].transpose[0]
				spectrum = GSL::Vector::Complex.alloc(phi_t[ky_element]).forward.square
				if options[:no_kpar0]
					spectrum[0]=0.0
				end
        #ep spectrum.size
				spectrum = spectrum.from_box_order
				#ep spectrum.shape
				spectrum = spectrum*gvky[ky_element]**2 unless options[:phi2_only]
        #ep gm.size
        #ep spectrum.size
				gm.set_row(ky_element, spectrum)
			end
			if options[:no_zonal]
				gm.row(0).set_all(0.0)
			end
			if options[:log]
				gm = gm.log
			end
			
			return gm
	end
	end

	def phi0_over_x_over_y_gsl_matrix(options)
	Dir.chdir(@directory) do

			options.convert_to_index(:t) if options[:t] or options[:t_element]
			options[:t_index] ||= list(:t).keys.max
			#phi2_by_mode index order (in Fortran) is kx, ky, t
			phi_re_narray = netcdf_file.var("phi0").get('start' => [0, 0, 0, options[:t_index] - 1], 'end' => [0, -1, -1, options[:t_index] - 1])
# 			ep phi_re_narray.shape
			phi_re_narray.reshape!(*phi_re_narray.shape.slice(1..2))
			# The narray has index order ky, kx, but we want kx, ky for historical reasons, hence the transpose. 
			gm_re = phi_re_narray.to_gm
			phi_im_narray = netcdf_file.var("phi0").get('start' => [1, 0, 0, options[:t_index] - 1], 'end' => [1, -1, -1, options[:t_index] - 1])
			phi_im_narray.reshape!(*phi_im_narray.shape.slice(1..2))
			# The narray has index order ky, kx, but we want kx, ky for historical imasons, hence the transpose. 
			gm_im = phi_im_narray.to_gm
			gm = GSL::Matrix::Complex.re_im(gm_re, gm_im)
# 			ep gm.shape

			if options[:no_zonal]
				
				for i in 0...gm.shape[1]
					gm[0,i] = GSL::Complex.alloc([0,0])
				end
			end
			if xres = (options[:xres] or options[:x_resolution])
				mat = GSL::Matrix::Complex.calloc(gm.shape[0], xres)
				extra = ((xres - gm.shape[1])).floor
				for i in 0...gm.shape[0]
					for j in 0...((gm.shape[1] + 1) / 2 )
# 						j+= extra if 
						mat[i, j] = gm[i,j]
						mat[i, j+extra] = gm[i,-j] unless j==0
					end
				end
				gm = mat
				
				
# 				gm = mat.vertcat(gm).vertcat(mat)
			end
			if yres = (options[:yres] or options[:y_resolution])
				mat = GSL::Matrix::Complex.calloc(yres, gm.shape[1])
				extra = ((yres - gm.shape[0])).floor
				for i in 0...gm.shape[0]
					for j in 0...gm.shape[1]
# 						j+= extra if 
						mat[i, j] = gm[i,j]
# 						mat[i, j+extra] = gm[i,-j] unless j==0
					end
				end
				gm = mat
				
				
# 				gm = mat.vertcat(gm).vertcat(mat)
			end
# 			ep gm_re, gm_im
# 			re = GSL::Complex.alloc([1.0, 0.0])
# 			gm = GSL::Matrix::Complex.calloc(*gm_re.shape)
# 			gm = gm_re * re  + gm_im * GSL::Complex.alloc([0.0, 1.0])
# 			gm_re, gm_im = fourier_transform_gm_matrix_complex_rows(gm_re, gm_im)
			
			gm = gm.backward_cols_c2c(true).backward_rows_cc2r(true)  
			if options[:limit]
				for i in 0...gm.shape[0]
					for j in 0...gm.shape[1]
# 						j+= extra if 
						gm[i, j] = [[gm[i,j], options[:limit][0]].max, options[:limit][1]].min
# 						mat[i, j+extra] = gm[i,-j] unless j==0
					end
				end
			end
			return gm
	end
	end
end

include GSLMatrices

def kx_shift(options)
#	ep options
	return 0 unless @g_exb and @g_exb.abs > 0.0
	#p options
	return - list(:ky)[options[:ky_index]] * list(:t)[(options[:t_index] or list(:t).keys.max)] * @g_exb
end

def jump(options)
#	ep 'kx_shift',  kx_shift(options)
	jump =  ((kx_shift(options) / list(:kx)[2]).round)
	case options[:t_index]
	when 1
		return jump
	else
		if @g_exb and @g_exb.abs > 0
			return jump + 1
		else
			return 0
		end
	end
end

# This function is used in the presence of perpendicular flow shear. It returns the (Eulerian) GS2
# kx_index as a function of the Lagrangian kx,  which is the kx_index of the mode in a shearing 
# coordinate system, I.e. if you give it an Lagrangian kx (which is the same as the Eulerian
# kx at t=0) it will tell you where it has now got to. It may have left the box, in which case
# this function will return an error.
#
# A given Lagrangian kx moves through the GS2 box, and thus for such a kx the response matrix varies
# in time. This is done because the effect of flow shear can be reduced by a shearing coordinate 
# transformation to become merely a time varying kx.
#
# At each timestep, phi(ikx_indexed(it)) is set equal to phi(ikx_indexed(it - jump(iky))
# kx_indexed is defined in the following way. 
#   do it=itmin(1), ntheta0
#			ikx_indexed (it+1-itmin(1)) = it
#		end do
#          
#		do it=1,itmin(1)-1
#			ikx_indexed (ntheta0 - itmin(1) + 1 + it)= it
#		end do
#
# In other words, what this means is that akx(ikx_indexed(0)) is the minimum kx, 
# and that akx(ikx_indexed(ntheta0)) gives the maximum kx, kx_indexed moves the
# kxs out of box order.
#
# So. remembering that jump is negative, phi(kx) is set equal phi(kx - jump * dkx)
# so the Lagrangian mode has moved to a lower kx. So get the Eulerian index, one
# starts with the Lagrangian index, and adds jump (which is negative!). This, however,
# must be done with indexes that are in the physical (not box) order. So this function
# first moves the indexes out of box order, then adds jump, then moves them back 
# into box order so that the index returned will give the correct kx from the GS2
# array.

def eulerian_kx_index(options)
	#eputs "Start eulerian_kx_index"
	lagrangian_kx_index = options[:kx_index]
	phys = physical_kx_index(lagrangian_kx_index)
	#ep 'jump', jump(options)
	index = phys + jump(options)
	raise ArgumentError.new("Lagrangian kx out of range") if index <= 0
	box= box_kx_index(index)
	#eputs "End eulerian_kx_index"
	return box
end

def kx_indexed
	return cache[:kx_indexed] if cache[:kx_indexed]
	#kx = cache[:kx_array] ||= gsl_vector('kx').to_a
	#kxphys = kx.from_box_order
	#min_index = kx.min_index + 1
	#cache[:kx_indexed] ||= kx.size.times.inject({}) do |hash, kx_element|
		#hash[kx_element + 1] = kxphs
	kx = gsl_vector('kx')
	size = kx.size
	box =  GSL::Vector::Int.indgen(size) + 1
	zero_element = kx.abs.min_index
	phys = box.subvector(zero_element, size-zero_element).connect(box.subvector(0, zero_element))
	cache[:kx_indexed] = [phys.to_a, box.to_a].transpose.inject({}){|hash, (phys, box)| hash[phys] = box; hash}
end

def box_kx_index(physical_kx_index)

	return kx_indexed[physical_kx_index]
end

def physical_kx_index(box_kx_index)
	return kx_indexed.key(box_kx_index)
	kx = cache[:kx_gslv] ||= gsl_vector('kx')
	return kx.from_box_order.to_a.index(kx[box_kx_index-1]) + 1
	#kx = cache[:kx_gslv] ||= gsl_vector('kx')
	#index_of_min_kx = cache[:index_of_min_kx] ||= kx.min_index + 1 # kx.min_index returns a 0-based index
	#if box_kx_index < index_of_min_kx
		#box_kx_index + (1 + kx.size - index_of_min_kx)
	#else
		#box_kx_index - (index_of_min_kx - 1)
	#end
end


	

def gsl_complex(name, options={})
	options = eval(options) if options.class == String
# 	p @directory
	Dir.chdir(@directory) do
# 		eputs Dir.pwd
		case name
		when /correcting_phase/
# 			options.convert_to_index(self, :ky)
# 			theta0 = (options[:theta0] or 0)
# # 			p 'options[:ky_index]', options[:ky_index]
# 			phase_array = NumRu::NetCDF.open("#@directory/#@run_name.out.nc").var('phase').get({"start" => [0, options[:ky_index] - 1, theta0], 'end' => [1, options[:ky_index] - 1, theta0] }).to_a.flatten
# 			p 'phase_array', phase_array
# 			thetaelement0 = (list(:theta).key(0.0) - 1).to_i
# 			p 'list(:theta)[thetaelement0 + 1]', list(:theta)[thetaelement0 + 1]
# 			p 'thetaelement0', thetaelement0
# 			p 'theta0 - jump(options)', theta0 - jump(options) % @jtwist
# 			p 'list(:kx)[2] * (theta0 - jump(options)%@jtwist)', list(:kx)[2] * (theta0 - jump(options)%@jtwist)
# 			kx_element = list(:kx).key(list(:kx)[2] * (theta0 - jump(options)%@jtwist)) - 1  
# 			at_0 = NumRu::NetCDF.open("#@directory/#@run_name.out.nc").var('phi').get({"start" => [0, thetaelement0, kx_element, options[:ky_index] - 1], 'end' => [1, thetaelement0, kx_element, options[:ky_index] - 1] }).to_a.flatten
# 			p 'at_0', at_0
# 			at_0 = GSL::Complex.alloc(at_0)
# 			p 'at_0', at_0
# 			return (at_0 / at_0.mag).conj
# # 			pp 'theta0', theta0
# # 			pp phase_array[5][theta0]
# 			return GSL::Complex.alloc(phase_array)
# # 			new_options = options.dup
# 			new_options[:imrc] = :real
# 			thetas = gsl_vector('theta_along_field_line', new_options)
# 			at_0 = gsl_vector_complex('phi_along_field_line', new_options)[.to_a.index(0.0)]
# 			p at_0
			exit
		else
			raise CRError.new("Unknown gsl_complex requested: #{name}")
		end
	#  			eputs data; gets
	end
end		

# def gsl_matrix(name, options={})
# 	if options[:t_index] or options[:frame_index]
# 		return get_gsl_matrix(name, options)
# 	else
# 		return cache[[:gsl_vector, name, options]] ||= get_gsl_matrix(name, options)
# 	end
# end



end
end
