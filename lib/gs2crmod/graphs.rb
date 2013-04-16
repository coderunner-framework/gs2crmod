

class CodeRunner
class Gs2 

# for backwards compatibility, and because I'm lazy
	# one day someone should get rid of this!
AxisKit = GraphKit::AxisKit
DataKit = GraphKit::DataKit

def auto_axiskits(name, options)
	hash = cache[:auto_axiskits] ||= {'t' => ['Time', ''],
                'phi2tot_over_time' => ['Phi^2 Total', ''],
                'apar2_over_time' => ['Apar^2 Total', ''],
                'growth_rate_by_ky_over_time' => ['Growth Rate by ky', ''],
                 'growth_rate_by_kx_over_time' => ['Growth Rate by ky', ''],  
								 'growth_rate_by_mode_over_time' => ["Growth Rate by mode", ''],
                'phi2_by_ky_over_time' => ['Phi^2 by ky', ''],
                 'phi2_by_kx_over_time' => ['Phi^2 by ky', ''],  
                'es_heat_by_ky_over_time' => ['Phi^2 by ky', ''],
                 'es_heat_by_kx_over_time' => ['Phi^2 by ky', ''],  
								 'phi2_by_mode_over_time' => ["Phi^2 by mode", ''],
						 		 'hflux_tot' => ['Total Heat Flux', ''],
                'ky' => ['ky', "1/rho_#{species_letter}"],
                'kx' => ['kx', "1/rho_#{species_letter}"],
	        'kpar' => ['kpar', "2 pi/qR"],
	        'growth_rate_over_kx' => ['Growth Rate', "v_th#{species_letter}/a", 1],
	        'growth_rate_over_ky' => ['Growth Rate', "v_th#{species_letter}/a", 1],
	        'transient_es_heat_flux_amplification_over_kx' => ['Transient Electrostatic Heat Amplification', "", 1],
	        'transient_es_heat_flux_amplification_over_ky' => ['Transient Electrostatic Heat Amplification', "", 1],
	        'transient_amplification_over_kx' => ['Transient Amplification', "", 1],
	        'transient_amplification_over_ky' => ['Transient Amplification', "", 1],
	        'spectrum_over_kx' => ["Spectrum at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 1],
	        'zonal_spectrum' => ["Zonal spectrum at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 1],
	        'spectrum_over_ky' => ["Spectrum at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 1],
	        'growth_rate_over_ky_over_kx' => ["Growth Rate", "v_th#{species_letter}/a", 2],
	       	'es_heat_flux_over_ky_over_kx' => ["Heat flux at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 2],
	       	'spectrum_over_kpar' => ["Spectrum at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 1],
	       	'spectrum_over_ky_over_kx' => ["Spectrum at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 2],
	       	'spectrum_over_ky_over_kpar' => ["Spectrum at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 2],
	        'phi0_over_x_over_y' => ["Phi at t = #{sprintf("%.3f" ,(options[:t] or list(:t)[options[:t_index]] or list(:t).values.max))}", '', 2],
	        'es_mom_flux_over_time' => ["#{species_type((options[:species_index] or 1)).capitalize} Momentum Flux", '', 1]

                  
                   }
	return hash[name]
end

def axiskit(name, options={})
	logf :axiskit
	if info = auto_axiskits(name, options)
		if info[2] and info[2] == 2
			axis =  GraphKit::AxisKit.autocreate({data: gsl_matrix(name, options), title: info[0], units: info[1]})
		elsif !info[2] or info[2] == 1
			axis =  GraphKit::AxisKit.autocreate({data: gsl_vector(name, options), title: info[0], units: info[1]})
			log 'successfully created axis'
		end
		return axis
	end
	case name
	when 'phi_along_field_line'
		title = options[:imrc].to_s.capitalize + " Phi"
		units = ""
		return GraphKit::AxisKit.autocreate(data: gsl_vector(name, options), title: title, units: units)
	when 'theta_along_field_line'
		title =  options[:z] ? "z/l_B" : 'Theta' 
		units = options[:z] ? '' : 'radians'
		return GraphKit::AxisKit.autocreate(data: gsl_vector(name, options), title: title, units: units)
	when 'es_heat_flux'
		type = species_type(options[:species_index]).capitalize
		units = ''
		return GraphKit::AxisKit.autocreate(data: gsl_vector('es_heat_flux_over_time', options), title: "#{type} Heat Flux", units: units)
# 	when 'spectrum_by_ky'
# 		return AxisKit.autocreate(data: gsl_vector('spectrum_by_ky', options), title: "Phi^2 at t = #{list(:t)[options[:t_index]]}", units: '')
	end
	raise CRError.new("Unknown axis kit: #{name}")
end

def self.cache
	@cache ||= {}
	@cache
end
	
def self.generate_graphs_rdoc_file
	File.open('graphs_rdoc.rb', 'w') do |file|
	graphs = self.instance_methods.find_all{|m| m.to_s =~ /_graphkit$/}.sort_by{|m| m.to_s}
	run = new(nil)
	file.puts "class #{self.to_s}::GraphKits\n"
	graphs.each do |graph|
		help = run.send(graph, command: :help)
		options = run.send(graph, command: :options)
		file.puts "# #{help}"
		if options and options.size > 0
			file.puts "# Options:"
			options.each do |op|
				file.puts "#\n# #{op}: #{GRAPHKIT_OPTIONS_HELP[op]}"
			end
		end
		file.puts "def #{graph}\nend"
	end
	file.puts "end"
	end
end
def self.help_graphs
# 	@@runner ||= CodeRunner.fetch_runner(U: true, 
	graphs = self.instance_methods.find_all{|m| m.to_s =~ /_graphkit$/}.sort_by{|m| m.to_s}
	run = new(nil)
	puts "-------------------------------------------\n    Available Graphs For #{self.to_s}\n-------------------------------------------\n"
	graphs.each do |graph|
		help = run.send(graph, command: :help)
		options = run.send(graph, command: :options)
		puts "\n------------------------------------\n#{graph.to_s.sub(/_graphkit/, '')}\n------------------------------------\n\n#{help}"
		if options and options.size > 0
			puts "\n\tOptions:"
			options.each do |op|
				puts "\t\t#{op}: #{GRAPHKIT_OPTIONS_HELP[op]}"
			end
		end
		
	end
end

GRAPHKIT_OPTIONS_HELP = {
	t_index_window: "[begin, end], window of time indices to plot (e.g. t_index_window: [0,10])",
	t_index: "integer, index of time at which to plot (e.g. t_index: 20)",
	t: "float, value of time at which to plot (e.g. t: 2.45)",
	ky_index: "integer, index of ky at which to plot (e.g. ky_index: 20)",
	ky: "float, value of ky at which to plot (e.g. ky: 0.1)",
	kx_index: "integer, index of kx at which to plot (e.g. kx_index: 20)",
	kx: "float, value of kx at which to plot (e.g. kx: 0.1)",
	with: "Gnuplot Option (may not apply when using other packages), e.g. with: 'lp' or with 'pm3d palette'",
	rgbformulae: "Gnuplot Option (may not apply when using other packages), sets colour mapping. See gnuplot help set rgbformulae",
	limit: "Limit the range of quantity begin plotted - any values of the quantity outside the limits will be set to the limit: eg. limit: [0,80]",
	flip: 'Flip the y axis,  e.g. flip: true',
	rev: 'Reverse the x axis, e.g. rev: true',
	z: 'Plot quantities vs z = theta/shat rather than theta. See Beer, Cowley Hammet 1996, eg. z: true',
	norm: 'Normalise the graph so that its maximum is 1, e.g. norm: true',
	mag: 'Plot the magnitude, e.g. mag: true',
	species_index: "Which GS2 species to plot the graph for (1-based).",
  strongest_non_zonal_mode: "Plot the graph requested for the mode with the highest value of phi^2. Overrides ky, kx, ky_index, kx_index. Can be set true or false; e.g. strongest_non_zonal_mode: true",
	no_zonal: "Don't plot the ky=0 part (boolean, e.g. no_zonal: true)",
	no_kpar0: "Don't plot the kpar=0 part (boolean, e.g. no_kpar0: true)",
	log: "Plot the log of a given quantity (exact meaning varies). boolean",
	Rmaj: "The major radius in metres. This has no effect on the shape of the graph: it merely multiplies every length",
 n0: " The toroidal mode number of the longest y mode. In effect it is the number of periodic copies of the flux tube that will fit in the torus. Periodicity requires that n0 q  is also an integer. If you specify :n0 where this is not the case, q will automatically be adjusted until it is",
 rho_star: " The ratio of the reference Lamour radius to the GS2 normalising length a. Cannot be specified at the same time as n0. If specified, both n0 and q will be adjusted to ensure periodicity",
 t_index: "The (1-based) time index",
 nakx: "The number of radial wave numbers to include in the plot. In effect, it is a low pass filter which reduces the resolution in the radial direction without changing the shape of the final surface. Minimum value is 4",
 naky: "The number of kys to include in the plot. In effect, it is a low pass filter which reduces the resolution in the y direction without changing the shape of the final surface. Minimum value is 4",
 gs2_coordinate_factor: "When set to 1, plot the graph in GS2 coordinates. When set to  0 plot the graph in real space. Can be set at any value between 0 and 1: the graph will smoothly distort between the two limits",
 xmax: "The (0-based) index of the maximum value of x to include in the plot",
 xmin: "The (0-based) index of the minimum value of x to include in the plot",
 ymax: "The (0-based) index of the maximum value of y to include in the plot",
 ymin: "The (0-based) index of the minimum value of y to include in the plot",
 thetamax: "The (0-based) index of the maximum value of theta to include in the plot",
 thetamin: "The (0-based) index of the minimum value of theta to include in the plot",
 ncopies: " The number of periodic copies of the flux tube to include",
 torphi_values: "An array of two values of the toroidal angle. The graph will be plotted in between those two values with poloidal cross sections at either end",
 magnify: " The magnification factor of the small section. It can take any value greater than or equal to 1",

}
	

# def graphkit(name, options={})
# 	unless [:Failed, :Complete].include? status
# 		return get_graphkit(name, options)
# 	else
# 		return cache[[:graphkit, name, options]] ||= get_graphkit(name, options)
# 	end
# end


def graphkit(name, options={})
	logf :graphkit
	# If an array of t, kx or ky values is provided, plot one graph for each value and then sum the graphs together
	[:t, :kx, :ky].each do |var|
		#ep 'index', var
		if options[var].class == Symbol and options[var] == :all
			options[var] = list(var).values
		elsif options[var+:_index].class == Symbol and options[var+:_index] == :all
			#ep 'Symbol'
			options[var+:_index] = list(var).keys
		end
		if options[var].class == Array
			return options[var].map{|value| graphkit(name, options.dup.absorb({var =>  value}))}.sum
		elsif options[var+:_index].class == Array
			#ep 'Array'
			return options[var+:_index].map{|value| graphkit(name, options.dup.absorb({var+:_index =>  value}))}.sum
		end
		if options[var].class == Symbol and options[var] == :max
			options[var] = list(var).values.max
		elsif options[var+:_index].class == Symbol and options[var+:_index] == :max
			ep 'Symbol'
			options[var+:_index] = list(var).keys.max
		end
	end
	options[:t_index] ||= options[:frame_index]  if options[:frame_index]

	


	# If a method from the new GraphKits module can generate this graphkit use it 
	#ep name + '_graphkit'
 	#ep self.class.instance_methods.find{|meth| (name + '_graphkit').to_sym == meth}

	if method = self.class.instance_methods.find{|meth| (name + '_graphkit').to_sym == meth}
		options[:graphkit_name] = name
		return send(method, options)
	end

	raise "Graph #{name} not found"
	
end

module GraphKits
		
	def apar2_by_mode_vs_time_graphkit(options={})
		options[:direction] = :mode
		apar2_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def apar2_by_kx_vs_time_graphkit(options={})
		options[:direction] = :kx
		apar2_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def apar2_by_ky_vs_time_graphkit(options={})
		options[:direction] = :ky
		apar2_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def apar2_by_kxy_or_mode_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "'apar2_by_ky_vs_time' or 'apar2_by_kx_vs_time' or 'apar2_by_mode_vs_time': Apar^2 over time for a given kx or ky, integrated over the other direction, or apar^2 vs time for a given kx and ky"
		when :options
			return  [:ky, :ky_index, :kx, :kx_index]
		else
			kxy = options[:direction]
	
			# i.e. apar2_by_ky_vs_time or apar2_by_kx_vs_time or apar2_by_mode_vs_time
			
			nt_options = options.dup # 'no time' options
			nt_options.delete(:t_index) if nt_options[:t_index]
			nt_options.delete(:frame_index) if nt_options[:frame_index]
			phiax = axiskit("apar2_by_#{kxy}_over_time", nt_options)	
			kit = GraphKit.autocreate({x: axiskit('t', options), y: phiax})
			kit.data[0].title = "Phi^2 total: #{kxy} = #{options[kxy]}"	
			if options[:t_index]
# 				p 'hello'
				array_element = options[:t_index_window] ? options[:t_index] - options[:t_index_window][0] : options[:t_index] - 1
# 				p phiax.data.size, array_element
# 				p options[:t_index], options[:t_index_window]
				time = DataKit.autocreate({x: {data: GSL::Vector.alloc([list(:t)[options[:t_index]]])}, y: {data: GSL::Vector.alloc([phiax.data[array_element]]) } })
				time.pointsize = 3.0
# 				p time
# 				kit.data[0].axes[:x].data = -kit.data[0].axes[:x].data
				kit.data.push time
				if options[:with_moving_efn] and kxy==:ky
					tmax = kit.data[0].axes[:x].data[-1]
# 					p 'tmax', tmax
					
					theta_max = @g_exb * tmax * options[:ky] * 2 * Math::PI / list(:kx)[2]
					kit.each_axiskit(:x) do |axiskit|
# 						p axiskit
						axiskit.data = axiskit.data / tmax * theta_max - theta_max
					end
				end
			end
			if options[:norm]
				xrange, yrange = kit.plot_area_size
				kit.each_axiskit(:y) do |axiskit|
					axiskit.data /= yrange[1] / (options[:height] or 1.0)
				end
			end
			kit.log_axis = 'y'
			#kit.data[0].title = "gs2:#@run_name"
			kit.data[0].with = "l" #"linespoints"
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def efnim_graphkit(options={})
		options[:imrc] = :im
		efn_graphkit(options)
	end


	def efnmag_graphkit(options={})
		options[:imrc] = :mag
		efn_graphkit(options)
	end

	def efn_graphkit(options={})
		case options[:command]
		when :help
			return  "Plot the eigenfunction along the extended domain. Options mag, norm, z can be specified by using a short hand in the name of the graph, eg. efnmagnormz, efnmag, efnnorm etc. If the range is set to 0, it plots the whole eigenfunction. Otherwise it plot a small bit of it. Only specify kx or kx_index if magnetic shear is 0."
		when :options
			return	[:mag, :norm, :z, :flip, :range, :kx_index, :ky_index, :kx, :ky, :strongest_non_zonal_mode]
		when :plot, nil
				eputs "Starting efn, this can take a while..."
				options[:imrc] ||= :real
				ep options
				options.convert_to_index(self, :ky)
				#ep 'converted options: ', options
				
				#decode naming scheme
				#mag = nil
				#options[:imrc] = :real
				#case name
				#when /im/
					#options[:imrc] = :im
				#when /mag/
					#options[:imrc] = :mag
					#options[:mag] = true
				#when /corr/
					#options[:imrc] = :corr
				#end 
				#options[:flip] = true if name =~ /flip/
				#options[:norm] = true if name =~ /norm/
				#options[:rev] = true if name =~ /rev/
				#options[:z] = true if name =~ /z/
				
				
				kit = GraphKit.autocreate({x: axiskit('theta_along_field_line', options), y: axiskit('phi_along_field_line', options)})
# 	                        ep kit
				kit.data[0].title =  "gs2:efn#{options[:imrc]}:#@run_name"
				kit.title = "#{options[:mag] ? "Magnitude of" : ""} Eigenfunction for ky=#{list(:ky)[options[:ky_index]]}, g_exb=#{@g_exb.to_f.to_s}, shat=#{@shat.to_s}"
				kit.file_name = "efn_for_#@run_name"
	# 			kit.pointsize = 1.0
				kit.modify(options)
				kit.title += sprintf(" time = %3.1f", list(:t)[options[:t_index]]) if options[:t_index]
				kit.data[0].with = "linespoints"
	# 			kit.data[0].axes[:x].data *= -1 #if options[:rev]
				#(eputs 'reversing'; gets)
				if (@s_hat_input||@shat).abs >= 1.0e-5
					range = options[:range] == 0 ? nil : (options[:range] or options[:z] ? 1 / (@s_hat_input||@shat) : 2 * Math::PI / (@s_hat_input||@shat))
					kit.xrange = [-range, range] if range
				end
				return kit
		end
	end

	alias :eigenfunction_graphkit :efn_graphkit
		
	def es_heat_flux_vs_ky_vs_kx_graphkit(options={})
		case options[:command]
		when :help
			return  "Graph of electrostatic contribution to heat flux at a given time vs kx and ky"
		when :options
			return [:with]
		else
			zaxis = axiskit('es_heat_flux_over_ky_over_kx', options)
			zaxis.data = zaxis.data.transpose
			kit = GraphKit.autocreate({y: axiskit('ky', options), x: axiskit('kx', options), z: zaxis})
			kit.title = "Heat flux"
			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.file_name = options[:graphkit_name]
			return kit
		end
	end

	def es_heat_flux_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "Heat flux vs time for each species."
		when :options
			return [:t_index_window, :species_index]
		else
			kit = GraphKit.autocreate({x: axiskit('t', options), y: axiskit('es_heat_flux', options)})
			kit.data[0].title = "#{species_type(options[:species_index])} hflux: #@run_name"
			kit.data[0].with = "l" #"lines"
			kit.file_name = options[:graphkit_name]
			kit
		end 
	end

	def es_mom_flux_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "Momentum flux vs time for each species."
		when :options
			return [:t_index_window, :species_index]
		else
			kit = GraphKit.autocreate({x: axiskit('t', options), y: axiskit('es_mom_flux_over_time', options)})
			kit.data[0].title = "#{species_type(options[:species_index])} momflux: #@run_name"
			kit.data[0].with = "l" #"lines"
			kit.file_name = options[:graphkit_name]
# 			kit.log_axis = 'y'
			return kit
		end
	end

	def transient_es_heat_flux_amplification_vs_kx_graphkit(options={})
		options[:kxy] = :kx
		transient_es_heat_flux_amplification_vs_kxy_graphkit(options)
	end

	def transient_es_heat_flux_amplification_vs_ky_graphkit(options={})
		options[:kxy] = :ky
		transient_es_heat_flux_amplification_vs_kxy_graphkit(options)
	end

	def transient_es_heat_flux_amplification_vs_kxy_graphkit(options={})
		case options[:command]
		when :help
			return "transient_es_heat_flux_amplification_vs_ky or transient_es_heat_flux_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_es_heat_flux_amplification_vs_kx_vs_ky/. "
		when :options
			return []
		else
			#raise "Growth Rates are not available in non-linear mode" if @nonlinear_mode == "on"
			kxy = options[:kxy]
			kit = GraphKit.autocreate({x: axiskit(kxy.to_s, options), y: axiskit("transient_es_heat_flux_amplification_over_#{kxy}", options)})
			kit.title  = "Transient Amplification of the ES Heat flux for species #{options[:species_index]} by #{kxy}"
			kit.data[0].with = "lp"
			kit.data[0].title = @run_name
			kit.file_name = options[:graphkit_name]
			kit
		end 
	end

	def transient_amplification_vs_kx_graphkit(options={})
		options[:kxy] = :kx
		transient_amplification_vs_kxy_graphkit(options)
	end

	def transient_amplification_vs_ky_graphkit(options={})
		options[:kxy] = :ky
		transient_amplification_vs_kxy_graphkit(options)
	end

	def transient_amplification_vs_kxy_graphkit(options={})
		case options[:command]
		when :help
			return "transient_amplification_vs_ky or transient_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_amplification_vs_kx_vs_ky/. "
		when :options
			return []
		else
			#raise "Growth Rates are not available in non-linear mode" if @nonlinear_mode == "on"
			kxy = options[:kxy]
			kit = GraphKit.autocreate({x: axiskit(kxy.to_s, options), y: axiskit("transient_amplification_over_#{kxy}", options)})
			kit.title  = "Transient Amplification by #{kxy}"
			kit.data[0].with = "lp"
			kit.data[0].title = @run_name
			kit.file_name = options[:graphkit_name]
			kit
		end 
	end


	def growth_rate_vs_kx_graphkit(options={})
		options[:kxy] = :kx
		growth_rate_vs_kxy_graphkit(options)
	end

	def growth_rate_vs_ky_graphkit(options={})
		options[:kxy] = :ky
		growth_rate_vs_kxy_graphkit(options)
	end

	def growth_rate_vs_kxy_graphkit(options={})
		case options[:command]
		when :help
			return "growth_rate_vs_ky or growth_rate_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /growth_rate_vs_kx_vs_ky/. "
		when :options
			return []
		else
			raise "Growth Rates are not available in non-linear mode" if @nonlinear_mode == "on"
			kxy = options[:kxy]
			kit = GraphKit.autocreate({x: axiskit(kxy.to_s, options), y: axiskit("growth_rate_over_#{kxy}", options)})
			kit.title  = "Growth Rates by #{kxy}"
			kit.data[0].with = "lp"
			kit.data[0].title = @run_name
			kit.file_name = options[:graphkit_name]
			kit
		end 
	end

	def growth_rate_vs_kx_vs_ky_graphkit(options={})
		case options[:command]
		when :help
			return "3D plot of growth rates vs ky and kx for phi^2"
		when :options
			return [:with]
		else
			zaxis = axiskit('growth_rate_over_ky_over_kx', options)
			zaxis.data = zaxis.data.transpose
			kit = GraphKit.autocreate({y: axiskit('ky', options), x: axiskit('kx', options),  z: zaxis})
			kit.title = "Growth Rate by kx and ky"
			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.file_name = options[:graphkit_name]
			kit
		end
	end

	def hflux_tot_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return "Graph of total heat flux vs time. No options"
		when :options
			return []
		else
			kit = GraphKit.autocreate({x: axiskit('t', options), y: axiskit('hflux_tot', options)})
			kit.data[0].title = "hflux tot:#@run_name"
			kit.data[0].with = "l" #"lines"
			kit.file_name = options[:graphkit_name]
# 			kit.log_axis = 'y'
			kit
		end
	end

	def kpar_spectrum_graphkit(options={})
		case options[:command]
		when :help
			return  "Graph of the k_parallel at a given kx and ky"
		when :options
			return  [:kx, :ky, :strongest_non_zonal_mode]
		else
			kit = GraphKit.autocreate({x: axiskit('kpar', options), y: axiskit('spectrum_over_kpar', options)})
			kit.data[0].title = "Spectrum at t = #{list(:t).values.max}"
			kit.file_name = options[:graphkit_name]
			kit.data[0].with = 'lp'
			kit
		end
	end

	def kx_spectrum_graphkit(options={})
		options[:kxy] = :kx
		kxy_spectrum_graphkit(options)
	end

	def ky_spectrum_graphkit(options={})
		options[:kxy] = :ky
		kxy_spectrum_graphkit(options)
	end

	def kxy_spectrum_graphkit(options={})
		case options[:command]
		when :help
			return "ky_spectrum or kx_spectrum: Graph of phi^2 vs kx or ky"
		when :options
			return  [:t, :t_index]
		else
			# ie ky_spectrum or kx_spectrum
			kxy = options[:kxy] 
			kit = GraphKit.autocreate({x: axiskit(kxy.to_s, options), y: axiskit("spectrum_over_#{kxy}", options)})
			kit.title  = "#{kxy} Spectrum"
			kit.file_name = options[:graphkit_name] + options[:t_index].to_s
			kit.data[0].with = 'lp'
			kit.ylabel = "Phi^2 #{kxy}^2"
			kit.pointsize = 2.0
			kit
		end
	end
	def lagrangian_kx_graphkit(options={})
		case options[:command]
		when :help
			return "A graph of the evolution of a single Lagrangian kx vs Eulerian kx and ky. Principally for debugging purposes"
		when :options
			return  []
		else
		  kyax = axiskit('ky', options)
			kyk = list(:ky).keys
			kx_list = list(:kx)
			begin
				kx_data = kyk.map do |ky_index| 
					options[:ky_index] = ky_index
					options[:kx_index] =1 
					ekx_index = eulerian_kx_index(options)
					kx_list[ekx_index]
				end
			rescue #ArgumentError
				kyk.pop
				retry
			end
			lky = list(:ky)
			kyax.data = kyk.map{|k| lky[k]}
			kit = GraphKit.autocreate(x: {data: kx_data, title: 'Eulerian kx (i.e. location on the GS2 grid)'}, y: kyax )
			kit.data[0].title = 'Lagrangian kx=0'
			kit.xrange = [kx_list.values.min, 0]
			kit.yrange = [0, lky.values.max]

			kit.title = 'Evolution of a Single Lagrangian kx vs ky'
			kit.file_name = 'lagrangian_kx_graphkit'
			return kit



		end
	end

	def phi_flux_tube_boundary_surface_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of cartesian coordinates, on one specified side of the flux tube (specified using the options :coordinate (:x, :y, :theta) and :side (:min, :max))"
		when :options
			return  [:Rgeo, :n0, :rho_star, :t_index, :nakx, :naky, :gs2_coordinate_factor, :xmax, :xmin, :ymax, :ymin, :thetamax, :thetamin, :ncopies, :side, :coordinate, :torphi_values]
		else
				phi = options[:phi] || phi_real_space_gsl_tensor(options)
			if options[:ncopies]
				ops = options.dup
				ops.delete(:ncopies)
				return (options[:ncopies].times.map do |n|
					ops[:ncopy] = n
					phi_flux_tube_boundary_surface_graphkit(ops)
				end).sum
			end
				side = options[:side]
				coord = options[:coordinate]
				shp = phi.shape
				#raise "  asdfal" unless shp
				opside = side == :max ? :min : :max
				ops = options.dup
				ops[coord+opside] = ops[coord+side] || 
					case side
					when :max
						case coord
						when :y
							shp[0]-1
						when :x
							shp[1]-1
						when :theta
							shp[2]-1
						end

						
					when :min
						0
					end
					opscut = ops.dup
					opscut[:ymin] = ((opscut[:ymin]||0) - (options[:ymin]||0))
					opscut[:ymax] = ((opscut[:ymax]||shp[0]-1) - (options[:ymin]||0))
					#ep 'opscut', opscut, (opscut[:xmax]||shp[1]-1), shp
					opscut[:xmin] = ((opscut[:xmin]||0) - (options[:xmin]||0))
					opscut[:xmax] = ((opscut[:xmax]||shp[1]-1) - (options[:xmin]||0))
					#ep 'opscut', opscut
					#ep 'opscut', opscut
					opscut[:thetamin] = ((opscut[:thetamin]||0) - (options[:thetamin]||0))
					opscut[:thetamax] = ((opscut[:thetamax]||shp[2]-1) - (options[:thetamin]||0))
				#p opscut;#gets
				if options[:cyl]
					coords = cylindrical_coordinates_gsl_tensor(ops) #.transpose(0,1,2,3)
				else
					coords = cartesian_coordinates_gsl_tensor(ops)
				end
				#ep ['coords', coords.shape, phi.shape]; gets
				newshape = coords.shape.slice(1..3)
				x = coords[0, false]; x.reshape!(*newshape)
				y = coords[1, false]; y.reshape!(*newshape)
				z = coords[2, false]; z.reshape!(*newshape)
				#ep shp , '...'
				range = [
					opscut[:ymin]||0..opscut[:ymax]||shp[0]-1, 
					opscut[:xmin]||0..opscut[:xmax]||shp[1]-1, 
					opscut[:thetamin]..opscut[:thetamax]

					#((ops[:thetamin]||0) - (options[:thetamin]||0))..((ops[:thetamax]||shp[2]-1) - (options[:thetamin]||0))
				]

				#ep ['range', range, 'phi.shape', phi.shape]
				phiside = phi[*range ]
				#ep ['coords', x.shape, phiside.shape]; gets
				kit = GraphKit.quick_create([x,y,z,phiside])

				# Create a map of this surface onto a continuous surface between two poloidal planes if ops[:torphi_values] is specified
				if torphi_values = ops[:torphi_values]
					raise "Can't take a poloidal cut at constant y or theta" if coord == :y #or coord == :theta
					raise "Can't take a poloidal cut with a limited y range (Remove :ymin and :ymax from ops)" if options[:ymin] or options[:ymax]

					torphi_values = torphi_values.sort
					cyls = cylindrical_coordinates_gsl_tensor(ops.absorb(extra_points: true)) 
					torphi_const0 = constant_torphi_surface_gsl_tensor(ops.absorb(torphi: torphi_values[0]))
					#ep torphi_const0, ops
					#torphi_const1 = constant_torphi_surface_gsl_tensor(ops.absorb(torphi: torphi_values[1]))
					raise "torphi_should be of rank 1: #{torphi_const0.shape}" unless torphi_const0.shape.include? 1
					
					#Get the number of points in the  y direction between the two poloidal planes
					deltorphi = cyls[2,-1,0,0] - cyls[2,0,0,0]
					i = i1 = istart = torphi_const0[0,0]
					displacement = torphi_values[0] - cyls[2,i,0,0]  # which copy of the flux tube are we in?
					#ep 'displacement', displacement
					shpc = cyls.shape
					ny = shpc[1]-1 # remove extra point
					while cyls[2,i%ny,0,0] + displacement < torphi_values[1]
					#ep ['displacement', displacement, 'torphi', torphi = cyls[2,i%ny,0,0] + displacement, 'i', i]
					 displacement=(cyls[2,i%ny+1,0,0]+displacement-cyls[2,0,0,0]) if (i+1)%ny == 0
						i+=1
					end
					ysize = i - i1
						#ep ['torphil', torphi = cyls[2,i%ny,0,0] + displacement, 'ysize', ysize]
					#exit


					
					#ysize = (torphi_const0[-1] - torphi_const1[0] ).abs + 1
					#case coord
					#when :x
					shpside = phiside.shape
					phicut = GSL::Tensor.float(ysize, shpside[1], shpside[2])
					xcut = GSL::Tensor.float(ysize, shpside[1], shpside[2])
					ycut = GSL::Tensor.float(ysize, shpside[1], shpside[2])
					zcut = GSL::Tensor.float(ysize, shpside[1], shpside[2])
					#ep shpside
					shpcut = phicut.shape
					for k in 0...shpcut[2] # 1 of these 2 loops 
					for j in 0...shpcut[1] # will have size 1

						istart = torphi_const0[j,k]
						displacement = torphi_values[0] - cyls[2,istart,j,k] # where in the torus does this copy of the flux tube start compared to where we want to be?

						#p torphi_const0[i,j], torphi_const1[i,j]
					for n in 0...shpcut[0] # index along our cut surface
						i = istart + n
						#icut = i%ysize
						#cyls = cylindrical_coordinates_gsl_tensor(ops.absorb(extra_points: true)) #.transpose(0,1,2,3)
					#ep ['displacement', displacement, n, i, ny, 'torphi', torphi = cyls[2,i%ny,j,k] +displacement]
						# We chop off the surface
						# between two y gridpoints
						# Hence we must interpolate
						# the field between those
						# gridpoints
						if n == 0
							torphi0 = torphi_values[0]
							s = Math.sin(torphi0)
							c = Math.cos(torphi0)
							d2 = cyls[2,(i+1)%ny,j,k] + displacement - torphi0 
							d1 = torphi0 - (cyls[2,i%ny,j,k] + displacement)
							dfac = d1 / (d1+d2)
						elsif n == ysize-1
							torphi1 = torphi_values[1]
							s = Math.sin(torphi1)
							c = Math.cos(torphi1)
							d2 = cyls[2,(i+1)%ny,j,k] + displacement - torphi1 
							d1 = torphi1 - (cyls[2,i%ny,j,k] + displacement)
							dfac = d1 / (d1+d2)
						else
							s = Math.sin(cyls[2,i%ny,j,k] +  displacement)
							c = Math.cos(cyls[2,i%ny,j,k] + displacement)
							dfac = 0
						end

						#ep [(phiside[i%ny,j,k] * (1-dfac) +  phiside[(i+1)%ny,j,k] * dfac).class, i,j,k, phiside.shape, ny  ]

						phicut[n,j,k] = phiside[i%ny,j,k] * (1-dfac) +  phiside[((i+1)%ny),j,k] * dfac
						
						rad = cyls[0,i%ny,j,k] 
						xcut[n,j,k] = rad * c
						ycut[n,j,k] = rad *s
						zcut[n,j,k] = cyls[1,i%ny,j,k] 
					  displacement=(cyls[2,i%ny+1,j,k]+displacement-cyls[2,0,j,k]) if (i+1)%ny == 0
						#displacement+=cyls[2,i,j,k] if (i+1)%ny == 0
					end
					#exit
					end
					end

					kit = GraphKit.quick_create([xcut,ycut,zcut,phicut])
				end






				kit.data[0].gp.with = "pm3d"
				return kit
		end

	end
	def phi_real_space_standard_representation_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of cartesian coordinates showing showing the standard way of representing the turbulence, with two poloidal cuts and the inner and outer radial surfaces. Multiple copies of the flux tube are used to fill the space."
		when :options
			return  [:Rgeo, :n0, :rho_star, :t_index, :nakx, :naky,  :xmax, :xmin, :thetamax, :thetamin,  :torphi_values]
		else
			phi = phi_real_space_gsl_tensor(options)
			options[:phi] = phi
			poloidal_planes = options[:torphi_values]
			#ep poloidal_planes; gets
			raise "Please set options[:torphi_values] to an array of two values" unless poloidal_planes.kind_of? Array and poloidal_planes.size==2
			options[:torphi] = poloidal_planes[0]
			kit1 = phi_real_space_poloidal_plane_graphkit(options)
			options[:torphi] = poloidal_planes[1]
			kit2 = phi_real_space_poloidal_plane_graphkit(options)
			options[:coordinate] = :x
			options[:side] = :min
			kit3 = phi_flux_tube_boundary_surface_graphkit(options)
			options[:side] = :max
			kit4 = phi_flux_tube_boundary_surface_graphkit(options)
			#kit= kit1+kit2
			#kit = kit3 +kit4
			kit = kit1+kit2+kit3+kit4
			if options[:thetamax] or options[:thetamin]
				options[:coordinate] = :theta
				options[:side] = :min
				kit5 = phi_flux_tube_boundary_surface_graphkit(options)
				options[:side] = :max
				kit6 = phi_flux_tube_boundary_surface_graphkit(options)
				kit += kit5+kit6
			end
			kit.gp.pm3d = "depthorder"
				kit.xlabel = 'X (m)'
				kit.ylabel = 'Y (m)'
				kit.zlabel = "'Z (m)' rotate by 90"
			kit.title = "3-D Potential"
			return kit
		end
	end
	def phi_real_space_poloidal_plane_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of cartesian coordinates showing a cut at one toroidal angle, with multiple periodic copies of the flux tube used to fill the whole circle.."
		when :options
			return  [:Rgeo, :n0, :rho_star, :t_index, :nakx, :naky,  :xmax, :xmin, :thetamax, :thetamin, :torphi]
		else
			#if options[:ncopies]
				#ops = options.dup
				#ops.delete(:ncopies)
				#return (options[:ncopies].times.map do |n|
					#ops[:ncopy] = n
					#phi_real_space_surface_graphkit(ops)
				#end).sum
			#end
			#zaxis = axiskit('phi0_over_x_over_y', options)
			#zaxis.data = zaxis.data.transpose
			#shape = zaxis.data.shape
			#carts = cartesian_coordinates_gsl_tensor(options)
			#torphiout = 2.6
			torphiout = options[:torphi]
			phi = options[:phi] || phi_real_space_gsl_tensor(options)
			torphi_const =  constant_torphi_surface_gsl_tensor(options)
			cyls = cylindrical_coordinates_gsl_tensor(options.absorb({extra_points: true}))
			#p torphi_const[0,true].to_a; 
			#p 'sh', cyls.shape[1], '','','',''; 
			#exit
			#carts = cartesian_coordinates_gsl_tensor(options)
			shp = phi.shape
			shpc = cyls.shape
			#ep 'shapes', shp, cyls.shape
			new_phi = GSL::Tensor.alloc(shp[1], shp[2])
			new_X = GSL::Tensor.alloc(shp[1], shp[2])
			new_Y = GSL::Tensor.alloc(shp[1], shp[2])
			new_Z = GSL::Tensor.alloc(shp[1], shp[2])
			#y = gsl_vector('y', options)
			#y = y.connect([2*y[-1] - y[-2]].to_gslv).dup
			#c = Math.cos(torphiout)
			#s = Math.sin(torphiout)
			#theta_vec = gsl_vector('theta', options)
			#x_vec = gsl_vector('x', options)
			#y_vec = gsl_vector('y', options)
			#phi.iterate do |i,j,k|
			#lastbracketed = nil
			#lastj = -1
			#phi.iterate_row_maj do |i,j,k|
			for k in 0...shp[2] #theta loop
			for j in 0...shp[1] #xloop
				#raise "Missed #{[j,k].inspect}, #{lastj}" unless lastj == j-1
				#ep [i,j,k], cyls.shape
				i = torphi_const[j,k]
				torphi1 = cyls[2,i,j,k]
				torphi2 = cyls[2,i+1,j,k]
				#ep cyls[2,true,j,k].to_a
				deltorphi = cyls[2,shpc[1]-1,j,k] - cyls[2,0,j,k]
				#raise "Periodicity not satisfied: #{(2*Math::PI/deltorphi+1.0e-8)%1.0}, #{(2*Math::PI/deltorphi+1.0e-5)}" unless ((2*Math::PI/deltorphi)+1.0e-8)%1.0 < 1.0e-5
				m1 = (torphi1 )%deltorphi 
				m2 = (torphi2 )%deltorphi 
				m3 = (torphiout)%deltorphi
				#bracketed = ((m1-m3).abs < 1.0e-4) || (
										#(m2-m3.abs) > 1.0e-4 &&
										#(m2 - m3) *
					    			#(m1 - m3) *
										#(m2 - m1) *
										#(torphi2 - torphi1) < 0)
				##p 'n0', (2*Math::PI/deltorphi).round
				#bracketed2 = (2*Math::PI/deltorphi + 1).round.times.inject(false) do |b,n|
					#epsn = 1.0e-4
					#eps2 = 1.0e-4
					#upp = torphiout + deltorphi * n
					#lwr = torphiout - deltorphi * n
					##measure = ((torphi1 < upp or (torphi1 - upp).abs < epsn) and  upp+epsn < torphi2) or ((torphi1 < lwr or (torphi1-lwr).abs < eps) and lwr + eps < torphi2)
					#a1 = a2 = a3 = b1 = b2 = b3 = 'q'
					##measure = ((
						###a1=((torphi1-upp).abs < (torphi2-upp).abs) and 
						##(a1=((torphi2-upp).abs > 1.0e-7)) and 
						##(a2=((torphi1 < upp or (torphi1 - upp).abs < epsn))) and  
						##(a3=(upp < torphi2))
					##) or (
						###b1=((torphi1-lwr).abs < (torphi2-lwr).abs) and 
						##(b1=((torphi2-lwr).abs > 1.0e7)) and 
						##(b2=(torphi1 < lwr or (torphi1-lwr).abs < epsn)) and 	
						##(b3 = (lwr  < torphi2))
					##))
						#a1=((torphi2-upp).abs > eps2)
						#a2=((torphi1 < upp or (torphi1 - upp).abs < epsn))
						#a3=(upp < torphi2)
						#b1=(torphi2-lwr).abs > eps2
						#b2=(torphi1 < lwr or (torphi1-lwr).abs < epsn)
						#b3 = (lwr  < torphi2)
					#measure = ((a1 and a2 and a3) or (b1 and b2 and b3))
					##ep 'measure', measure,  [torphi1, torphi2, upp, lwr , n, deltorphi, a1, a2, a3, b1, b2, b3, (torphi2-lwr).abs, (torphi2-lwr).abs > eps2] if measure #; gets if [j,k] == [5,8] # if measure
					#b or measure
				#end
				##bracketed = bracketed2
				#raise "Measures don't agree #{bracketed}, #{bracketed2}" unless bracketed2 == bracketed
#
#
				##d2 = torphi2 - torphiout
				##d1 = torphiout - torphi1
				#if bracketed
					#y1 = 
					d2 = m2 - m3
					d1 = m3 - m1
					if torphi2 > torphi1
						d1+=deltorphi if d1 < 0
					else
						d2-=deltorphi if d2 > 0
					end
					dfac = d1 / (d1+d2)

					#n = 0
					#loop do

					#ep [torphi1, torphi2, cyls[2,shpc[1]-1,j,k] , cyls[2,0,j,k], deltorphi, m1, m2, m3]
				#ep (mod2 - mod3)/(mod1 - mod3) < 0,"********" 
				#raise "Doubled up" if lastbracketed == [j,k]
				#raise "Missed: #{i},#{j}, #{lastbracketed.inspect} " unless lastbracketed == [j-1,k] or lastbracketed == [shp[1]-1, k-1] if lastbracketed
				#ep bracketed,"********" 
					#ep ['bracketed', i, j, k]
					#ep ['phi', phi[i,j,k]]
					#new_phi[j,k] = theta_vec[k]
					new_phi[j,k] = phi[i,j,k] * (1-dfac) +  phi[(i+1)%shp[0],j,k] * dfac
					#raise "Mismatched radii" unless
					#new_X[j,k] = cyls[0,i,j,k] * Math.cos(cyls[2,i,j,k])
					#new_X[j,k] =  x_vec[j]
					new_X[j,k] = cyls[0,i,j,k] * Math.cos(torphiout)
					#new_Y[j,k] = cyls[0,i,j,k] * Math.sin(cyls[2,i,j,k])
					#new_Y[j,k] = y_vec[i]
					new_Y[j,k] = cyls[0,i,j,k] * Math.sin(torphiout)
					#new_Z[j,k] = theta_vec[k]
					new_Z[j,k] = cyls[1,i,j,k]
					#lastbracketed = [j,k]
					#lastj = j
				#end
			end # xloop
			end # theta loop
			#exit
			kit =  GraphKit.quick_create([new_X, new_Y, new_Z, new_phi])
			kit.xlabel = 'X'
			kit.ylabel = 'Y'
			kit.data[0].gp.with = "pm3d"
			return kit








			#end

			sides = [:max, :min]
			kits = []
			#[].each_with_index do |coord,i|
				#sides.each_with_index do |side,j|
					[[:y, :min, 0], [:y, :max, Math::PI]].each do |coord, side, toroidalphi|
					#[0].each do |toroidalphi|
					#raise unless i.kind_of? Integer
					#return kit if coord + side == :xmax
					options[:phi] = phi
					options[:side] = side
					options[:coordinate] = coord
					options[:toroidal_projection] = coord == :x ? nil : toroidalphi
					next if coord == :x and toroidalphi == Math::PI
					kit = phi_flux_tube_boundary_surface_graphkit(options)
					kits.push kit
					end
				#end
			#end
			k = kits.sum
			#k = kits[5] + kits[4] + kits[3] + kits[2] + kits[1] + kits[0]
			#k = kits[5] + kits[4] + kits[1] + kits[0]
			k.gp.pm3d = "depthorder"
			#k.compress_datakits = true
			return k

		end
	end
	#def phi_flux_tube_poloidal_cut_graphkit(options={})
		######case options[:command]
		######when :help
			######return  "The potential as a function of cartesian coordinates, cut at two toroidal angles."
		######when :options
			######return  [:rgbformulae, :limit, :t_index]
		######end
		#poloidal_planes = options[:torphi_values]
		##ep poloidal_planes; gets
		#raise "Please set options[:torphi_values] to an array of two values" unless poloidal_planes.kind_of? Array and poloidal_planes.size==2
		#field = options[:phi] || phi_real_space_gsl_tensor(options)
		#options[:torphi] = poloidal_planes[0]
		#torphi_const1 = constant_torphi_surface_gsl_tensor(options.dup.absorb(no_copies: true))
		#options[:torphi] = poloidal_planes[1]
		#torphi_const2 = constant_torphi_surface_gsl_tensor(options.dup.absorb(no_copies: true))
		#x1 = SparseTensor.new(2) # First poloidal face
		#y1 = SparseTensor.new(2)  
		#z1 = SparseTensor.new(2)
		#field1 = SparseTensor.new(2)

		#x2 = SparseTensor.new(2)  # Second poloidal 
		#y2 = SparseTensor.new(2)  # face
		#z2 = SparseTensor.new(2)
		#field2 = SparseTensor.new(2)

		#xmaxx = SparseTensor.new(2) # Connecting 
		#ymaxx = SparseTensor.new(2) # surface of the
		#zmaxx = SparseTensor.new(2) # flux tube
		#fieldmaxx = SparseTensor.new(2) # max x
		#xminx = SparseTensor.new(2) # Connecting 
		#yminx = SparseTensor.new(2) # surface of the
		#zminx = SparseTensor.new(2) # flux tube
		#fieldminx = SparseTensor.new(2)
		#xmaxth = SparseTensor.new(2) # Connecting 
		#ymaxth = SparseTensor.new(2) # surface of the
		#zmaxth = SparseTensor.new(2) # flux tube
		#fieldmaxth = SparseTensor.new(2)
		#xminth = SparseTensor.new(2) # Connecting 
		#yminth = SparseTensor.new(2) # surface of the
		#zminth = SparseTensor.new(2) # flux tube
		#fieldminth = SparseTensor.new(2)

		#xsurf = [xmaxx, xminx, xmaxth, xminth]
		#ysurf = [ymaxx, yminx, ymaxth, yminth]
		#zsurf = [zmaxx, zminx, zmaxth, zminth]
		#fieldsurf = [fieldmaxx, fieldminx, fieldmaxth, fieldminth]

		#cyls = cylindrical_coordinates_gsl_tensor(options)
		#shp = cyls.shape
		#ysize = shp[1]

		## Find out how far each cut
		## poloidal face of the flux
		## tube extends in theta and x
		#k1min = GSL::Tensor.int(shp[2])
		#k1max = GSL::Tensor.int(shp[2])
		#k2min = GSL::Tensor.int(shp[2])
		#k2max = GSL::Tensor.int(shp[2])
		#k1min[true] = shp[3]; k1max[true] = 0
		#k2min[true] = shp[3]; k2max[true] = 0
		#j1min = GSL::Tensor.int(shp[3])
		#j1max = GSL::Tensor.int(shp[3])
		#j2min = GSL::Tensor.int(shp[3])
		#j2max = GSL::Tensor.int(shp[3])
		#j1min[true] = shp[2]; j1max[true] = 0
		#j2min[true] = shp[2]; j2max[true] = 0
		##ep j1max;gets
		#for j in 0...shp[2]
			#for k in 0...shp[3]
				#i1 = torphi_const1[j,k]
				#i2 = torphi_const2[j,k]
				#i1 = 0 if i1==-1
				#i2 = 0 if i2==-1
				#i1 = ysize-1 if i1==ysize
				#i2 = ysize-1 if i2==ysize
				#if i1%ysize==i1
					#k1min[j] = [k, k1min[j]].min
					#k1max[j] = [k, k1max[j]].max
					#j1min[k] = [j, j1min[k]].min
					#j1max[k] = [j, j1max[k]].max
				#end
				#if i2%ysize==i2
					#k2min[j] = [k, k2min[j]].min
					#k2max[j] = [k, k2max[j]].max
					#j2min[k] = [j, j2min[k]].min
					#j2max[k] = [j, j2max[k]].max
				#end
			#end
		#end
		#lineorder = SparseTensor.new(2)
		#surfaces = SparseTensor.new(2)
		##Now we generate the order of the line
		##scans between one cross section and 
		##the other (the lines have to be in order 
		## around the cross section)
		## We also specify which surface(s) the
		## line lies on
		#o = 0
		##min x surface
		#j=0
		#for k in ([k1min[0],k2min[0]].min)..([k1max[0], k2max[0]].max)
			#lineorder[j,k] = o
			#surfaces[j,k]||=[]
			#surfaces[j,k].push 1
			#o+=1
		#end
		## max theta surface
		#for j in 0...shp[2]
			#lineorder[j,[k1max[j],k2max[j]].max]=o
			#surfaces[j,[k1max[j],k2max[j]].max]||=[]
			#surfaces[j,[k1max[j],k2max[j]].max].push 2
			#o+=1
		#end
		#j = shp[2]-1
		## max x surface
		#ep ([k1max[j], k2max[j]].max)..([k1min[j],k2min[j]].min)
		#for k in (-[k1max[j], k2max[j]].max)...(-[k1min[j],k2min[j]].min)
			#ep 'k', k
			#lineorder[j,-k] = o
			#surfaces[j,-k]||=[]
			#surfaces[j,-k].push 0
			#o+=1
		#end
		## max theta surface
		#for j in -(shp[2]-1)..0
			#lineorder[-j,[k1min[j],k2min[j]].min]=o
			#surfaces[-j,[k1min[j],k2min[j]].min]||=[]
			#surfaces[-j,[k1min[j],k2min[j]].min].push 3
			#o+=1
		#end
		#ep lineorder, surfaces;gets

		##ep j1max.to_a;gets
		#for j in 0...shp[2]
			#for k in 0...shp[3]
				#i10 = i1 = torphi_const1[j,k]
				#i20 = i2 = torphi_const2[j,k]
				#ep ['i1', i1, 'i2', i2]

				## Include extra points just outside 
				## if necessary
				#i10 = 0 if i1==-1 #
				#i20 = 0 if i2==-1
				#i10 = ysize - 1 if i1 == ysize
				#i20 = ysize - 1 if i2 == ysize
				##if i1%ysize==i1 || i1==-1
				#if i10%ysize==i10
					## First assign all the points on the first cut poloidal face
					#r = cyls[0,i10,j,k]
					#z = cyls[1,i10,j,k]
					#tphi = poloidal_planes[0]
					##tphi = cyls[2,i1,j,k] 
					##ep ['tphi', tphi, 'i1', i1, cyls[2,i1,j,k]%(2*Math::PI), cyls[2,(i1+1)%ysize,j,k]%(2*Math::PI)]
					#s = Math.sin(tphi)
					#c = Math.cos(tphi)
					#x1p = x1[j,k] = r*c# X (not gs2 x!)
					#y1p = y1[j,k] = r*s# Y (not gs2 y!)
					#z1p = z1[j,k] = z# Z
					##ep [2,(i1+1)%ysize,j,k, cyls.shape]
					## I think this bit will be wrong
					## when i1 = -1 or ysize - 1
					## but it won't be a big error
					## --- need to fix
					##d2 = cyls[2,(i1+1)%ysize,j,k] - tphi
					##d1 = tphi - cyls[2,i1,j,k]
					##dfac = d1 / (d1+d2)
					##f1p = field1[j,k] = field[i1,j,k] * (1-dfac) +  field[i1%ysize,j,k] * dfac
					#f1p = (field1[j,k] = field[i1%ysize,j,k]  +  field[(i1+1)%ysize,j,k])/2.0
				#end
				#if i20%ysize==i20
					#r = cyls[0,i20,j,k]
					#z = cyls[1,i20,j,k]
					#tphi = poloidal_planes[1]
					#s = Math.sin(tphi)
					#c = Math.cos(tphi)
					#x2[j,k] = r*c# X (not gs2 x!)
					#y2[j,k] = r*s# Y (not gs2 y!)
					#z2[j,k] = z# Z
					##ep [2,(i2+1)%ysize,j,k, cyls.shape]
					##d2 = cyls[2,(i2+1)%ysize,j,k] - tphi
					##d1 = tphi - cyls[2,i2,j,k]
					##dfac = d1 / (d1+d2)
					#field2[j,k] = (field[i2%ysize,j,k]  +  field[(i2+1)%ysize,j,k])/2.0
					##field2[j,k] = field[i2,j,k] * (1-dfac) +  field[i2%ysize,j,k] * dfac
				#end
			#end
		#end
		#for j in 0...shp[2]
			#for k in 0...shp[3]
				#i10 = i1 = torphi_const1[j,k]
				#i20 = i2 = torphi_const2[j,k]
				##ep ['i1', i1, 'i2', i2]

				## Include extra points just outside 
				## if necessary
				#i10 = 0 if i1==-1 #
				#i20 = 0 if i2==-1
				#i10 = ysize - 1 if i1 == ysize
				#i20 = ysize - 1 if i2 == ysize
				#if lineorder[j,k]
					#n=0
					#sign21 = ((i2-i1)/(i2-i1).abs).round
					## Go from y on the first surface to
					## y on the next surface
					#ep ['jk',j,k, (i10+n*sign21)*sign21 < i20*sign21]
					#while (i10+n*sign21)*sign21 < i20*sign21
						#ka  = k
						#ia10 = i10
						#ia20 = i20
						#while (ia10+n*sign21)%ysize != ia10+n*sign21
							#ep [j,k,ka, 'ia10', ia10, ia10+n*sign21, lineorder[j,k], surfaces[j,k], 'n', n, 'ni', ia10 + n * sign21, 'ysize', ysize, 'max', k1max[j], k2max[j], 'min', k1min[j], k2min[j] ]
							## We have left the flux surface...
							## increase/decrease theta to get back in
							#if k > k1max[j] or k > k2max[j]  # thetamax
								#ka-=1
							#elsif k < k1min[j] or k < k2min[j]  # thetamin
								#ka+=1
							#else
								#raise "Got lost!"
							#end
							#ia10 = torphi_const1[j,ka]
							#ia10 = 0 if ia10==1
							#ia10 = ysize - 1 if ia10 == ysize
							#ia20 = torphi_const2[j,ka]
							#ia20 = 0 if ia20==1
							#ia20 = ysize - 1 if ia20 == ysize
						#end
						#ni = n*sign21+ia10
						## First end
						#if ni==ia10
							#xs = x1[j,ka]
							#ys = y1[j,ka]
							#zs = z1[j,ka]
							#fs  = field1[j,ka]
						## Second ed
						#elsif ni*sign21 == ia20*sign21 - 1
							#xs = x2[j,ka]
							#ys = y2[j,ka]
							#zs = z2[j,ka]
							#fs  = field2[j,ka]
						## Connecting line
						#else
							#tphi = cyls[2,ni,j,ka]
							#s = Math.sin(tphi)
							#c = Math.cos(tphi)
							#r = cyls[0,ni,j,ka]
							#xs = r*c
							#ys = r*s
							#zs = cyls[1,ni,j,ka]
							#fs = field[ni,j,ka] 
						#end


						#surfaces[j,k].each do |m|
							#o = lineorder[j,k]
							#xsurf[m][o,n] = xs
							#ysurf[m][o,n] = ys
							#zsurf[m][o,n] = zs
							#fieldsurf[m][o,n] = fs
						#end

						#n+=1
					#end
				#end


				#next

					## Fill in from the first poloidal plane to the second poloidal plane or the edge of the box, which ever is closer
					#if [j1min[k], j1max[k]].include? j or (false and [k1min[j],k1max[j]].include? k)
					#m = (j==j1min[k] ? 1 : 0)
					
					#xsurf[m][j,k,0] = r*c# X (not gs2 x!)
					#ysurf[m][j,k,0] = r*s# Y (not gs2 y!)
					#zsurf[m][j,k,0] = z# Z
					#fieldsurf[m][j,k,0] = field[i1,j,k] 
					#sign21 = ((i2-i1)/(i2-i1).abs).round
					##ep sign21; gets
					#n=1
					##while (i1+n*sign21)%ysize==i1+n*sign21  and (i1+n*sign21)*sign21 < i2*sign21
					#while (i10+n*sign21)%ysize==i10+n*sign21  and (i10+n*sign21)*sign21 < i20*sign21
						#ni = n*sign21+i10
						##ep [[2,ni,j,k], cyls.shape] 
						#tphi = cyls[2,ni,j,k]
						#s = Math.sin(tphi)
						#c = Math.cos(tphi)
						#r = cyls[0,ni,j,k]
						#z = cyls[1,ni,j,k]
						#xsurf[m][j,k,n] = r*c# X (not gs2 x!)
						#ysurf[m][j,k,n] = r*s# Y (not gs2 y!)
						#zsurf[m][j,k,n] = z# Z
						#fieldsurf[m][j,k,n] = field[ni,j,k]
						#n+=1
					#end
					#end
				##end
				##if i2%ysize==i2 || i2==-1
				#if i20%ysize==i20
					#r = cyls[0,i2,j,k]
					#z = cyls[1,i2,j,k]
					#tphi = poloidal_planes[1]
					#s = Math.sin(tphi)
					#c = Math.cos(tphi)
					#x2[j,k] = r*c# X (not gs2 x!)
					#y2[j,k] = r*s# Y (not gs2 y!)
					#z2[j,k] = z# Z
					##ep [2,(i2+1)%ysize,j,k, cyls.shape]
					#d2 = cyls[2,(i2+1)%ysize,j,k] - tphi
					#d1 = tphi - cyls[2,i2,j,k]
					#dfac = d1 / (d1+d2)
					#field2[j,k] = field[i2,j,k] * (1-dfac) +  field[i2%ysize,j,k] * dfac
				#end
				#next
					##if (i1%ysize==i1 || i1==-1) and ( [0,shp[2]-1].include? j or (false and [k1min[j],k1max[j]].include? k))
					#if i10%ysize==i10  and ( [j1min[k], j1max[k]].include? j or (false and [k1min[j],k1max[j]].include? k))
					#m = (j==j1min[k] ? 1 : 0)
						#xsurf[m][j,k,n] = x2[j,k]
						#ysurf[m][j,k,n] = y2[j,k]
						#zsurf[m][j,k,n] = z2[j,k]
						#fieldsurf[m][j,k,n] = field2[j,k]
					#elsif false
						#if [0,shp[2]-1].include? j or [k2min[j],k2max[j]].include? k
						#sign12 = ((i1-i2)/(i1-i2).abs).round
						#xouter[j,k,0] = x2[j,k]
						#youter[j,k,0] = y2[j,k]
						#zouter[j,k,0] = z2[j,k]
						#fieldouter[j,k,0] = field2[j,k]
						##ep sign12; gets

						## Fill in from the first poloidal plane to the second poloidal plane or the edge of the box, which ever is closer
						#n=1
						#while (i2+n*sign12)%ysize==i2+n*sign12  and (i2+n*sign12)*sign12 < i1*sign12
							#ni = n*sign12+i2
							##ep [[1,ni,j,k], cyls.shape] 
							#tphi = cyls[2,ni,j,k]
							#s = Math.sin(tphi)
							#c = Math.cos(tphi)
							#r = cyls[0,ni,j,k]
							#z = cyls[1,ni,j,k]
							#xouter[j,k,n] = r*c# X (not gs2 x!)
							#youter[j,k,n] = r*s# Y (not gs2 y!)
							#zouter[j,k,n] = z# Z
							#fieldouter[j,k,n] = field[ni,j,k]
							#n+=1
						#end
						#end

					##end
				#end
			#end
		#end
		#ep xsurf[3]
		#ep xsurf[2]
		#kit1 = GraphKit.quick_create([x1,y1,z1,field1])
		#kit2 = GraphKit.quick_create([x2,y2,z2,field2])
		#kits = 4.times.map{|i| GraphKit.quick_create([xsurf[i],ysurf[i],zsurf[i],fieldsurf[i]])}
		#return kit1 + kit2 + kits.sum
	#end





	def phi_magnifying_glass_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of cartesian coordinates, plus a close up section."
		when :options
			return  [:Rgeo, :n0, :rho_star, :t_index, :nakx, :naky, :xmax, :xmin, :thetamax, :thetamin, :magnify]
		end

		ops = options.dup
		ops.delete(:thetamin)
		ops.delete(:thetamax)
		flux_tube = phi_real_space_surface_graphkit(ops)
		options[:thetamin] #= options[:theta_magnify]
		options[:thetamax] #= options[:theta_magnify] + 4
		#options[:torphi_values] = [:options
			#correct_3d_options(options)
		#factors = geometric_factors_gsl_tensor(options)
			#xfac = 1.0 / options[:rho_star_actual]
			#yfac = options[:rhoc_actual] / options[:q_actual] / options[:rho_star_actual]	
			#y = gsl_vector('y', options)
			#x = gsl_vector('x', options)
			#kmin = 0 #options[:thetamin]
			#kmax = -1 #options[:thetamax]

			#ep 'factors.shape', factors.shape
			#torphi1 = y[-1]/yfac - factors[2,kmin] - x[-1]/xfac * factors[5,kmin]
			#torphi2 = factors[2,kmax]
			cyls = cylindrical_coordinates_gsl_tensor(options)
			#torphi1 = cyls[2,0,true,0].max
			torphi1 = cyls[2,0,true,-1].max
			torphi2 = cyls[2,-1,true,-1].min
			torphi3 = cyls[2,0,true,-1].max
			torphi4 = cyls[2,-1,true,0].min
			delta = 0.5
			if torphi1 > torphi3
				tv = [torphi1, [torphi1 + delta, torphi4].min]
			else
				tv = [torphi3, [torphi3 + delta, torphi2].min]
			end

			ep 'torphiVaues', options[:torphi_values] = tv; 


		#ep options; gets
		section = phi_real_space_standard_representation_graphkit(options)
		#dR1 = (section.data.map{|dk| 
			#(dk.x.data.max**2 + dk.y.data.max**2)**0.5
		#}.max +
		 #section.data.map{|dk|
			#(dk.x.data.min**2 + dk.y.data.min**2)**0.5
		#}.min)/2 
		#z1max = section.data.map{|dk| dk.z.data.max}.max 
		#z1min = section.data.map{|dk| dk.z.data.min}.min )/2

		magnify = options[:magnify]

		#p section; gets
		blacked_out = section.dup
		fmin = section.data.map{|dk| dk.f.data.min}.min 
		if magnify > 2
			blacked_out.data.each{|dk|
				dk.x.data *= 1.05
				dk.y.data *= 1.05
				dk.z.data *= 1.05
				dk.f.data.fill!(fmin)
				#dk.gp.with = 'p lc rgb "#000000"  ps 3.0 '
			}
		end
		section.data.each do |dk|
			magnify = options[:magnify]
			#ep dk.x.data.shape
			dk.x.data *= magnify
			dk.y.data *= magnify
			dk.z.data *= magnify
		end
		dR = (section.data.map{|dk| 
			(dk.x.data.max**2 + dk.y.data.max**2)**0.5
		}.max +
		 section.data.map{|dk|
			(dk.x.data.min**2 + dk.y.data.min**2)**0.5
		}.min)/2 

		#xmag = (section.data.map{|dk| dk.x.data.max}.max +
		 #section.data.map{|dk| dk.x.data.min}.min)/2 
		#ymag = (section.data.map{|dk| dk.y.data.max}.max +
		 #section.data.map{|dk| dk.y.data.min}.min)/2 

		#ep 'dR', dR; gets
		#dy = (section.data.map{|dk| dk.y.data.max}.max +
		 #section.data.map{|dk| dk.y.data.min}.min)/2 
		dz = (section.data.map{|dk| dk.z.data.max}.max +
		 section.data.map{|dk| dk.z.data.min}.min )/2
		section.data.each do |dk|
			#dk.x.data -= options[:Rgeo]*(magnify* (1.0-Math.exp(-(magnify-1)**2)))*c/1.2
			#dk.x.data -= (options[:Rgeo]-options[:rhoc_actual])*([magnify-2, 0].max)*c* (1.0-Math.exp(-(magnify-1)**2))
			ep 'rhoc', options[:rhoc_actual]
			#dR = options[:Rgeo]*(0.75*magnify - magnify * (1.0 - options[:rhoc_actual]/options[:Rgeo]))* (1.0-Math.exp(-(magnify-1)**2))
			c = Math.cos(tv[0])	
			s = Math.sin(tv[0])	
			#dR = options[:Rgeo]*(1.5 - magnify * (1.0 - options[:rhoc_actual]))* (1.0-Math.exp(-(magnify-1)**2))

			rout = 2.0
			dk.x.data -= (dR - options[:Rgeo] * rout) * c* (1.0-Math.exp(-(magnify-1)**2))
			#ep dk.x.data.shape; gets
			#dk.y.data -= (options[:Rgeo]-options[:rhoc_actual])*([magnify-2, 0].max)*s* (1.0-Math.exp(-(magnify-1)**2))
			dk.y.data -= (dR - options[:Rgeo] * rout) * s* (1.0-Math.exp(-(magnify-1)**2))
			dk.z.data-=dz* (1.0-Math.exp(-(magnify-1)**2))
		end
		#return section
		kit =    flux_tube + section
		kit.data.each{|dk| dk.gp.with = "pm3d"}
	  kit	+= blacked_out
		kit.data.each{|dk| dk.gp.with = "pm3d"}
		kit.xlabel = 'X'
		ep phideg =  tv[0]/Math::PI * 180
		kit.gp.view = [",#{(180-phideg)%360},2.0", "equal xyz"]
		kit.gp.hidden3d = ""
			kit.gp.pm3d = "depthorder"
				kit.xlabel = 'X (m)'
				kit.ylabel = 'Y (m)'
				kit.zlabel = "'Z (m)' rotate by 90"
			kit.title = "3-D Potential"
	
		kit
	end


	def phi_real_space_surface_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of cartesian coordinates, plotted on the six outer surfaces of constant x, y and theta."
		when :options
			return  [:Rgeo, :n0, :rho_star, :t_index, :nakx, :naky, :gs2_coordinate_factor, :xmax, :xmin, :ymax, :ymin, :thetamax, :thetamin, :ncopies]
		else
			if options[:ncopies]
				ops = options.dup
				ops.delete(:ncopies)
				return (options[:ncopies].times.map do |n|
					ops[:ncopy] = n
					phi_real_space_surface_graphkit(ops)
				end).sum
			end
			#zaxis = axiskit('phi0_over_x_over_y', options)
			#zaxis.data = zaxis.data.transpose
			#shape = zaxis.data.shape
			#carts = cartesian_coordinates_gsl_tensor(options)
			phi = options[:phi] || phi_real_space_gsl_tensor(options)
			sides = [:max, :min]
			kits = []
			[:y, :x, :theta].each_with_index do |coord,i|
				sides.each_with_index do |side,j|
					raise unless i.kind_of? Integer
					#return kit if coord + side == :xmax
					options[:phi] = phi
					options[:side] = side
					options[:coordinate] = coord
					kit = phi_flux_tube_boundary_surface_graphkit(options)
					kits.push kit
				end
			end
			k = kits.reverse.sum
			#k = kits[5] + kits[4] + kits[3] + kits[2] + kits[1] + kits[0]
			#k = kits[5] + kits[4] + kits[1] + kits[0]
			k.gp.pm3d = "depthorder"
			kit = k
			kit.data.each{|dk| dk.title = "notitle"}
			#k.compress_datakits = true
			if options[:gs2_coordinate_factor] and options[:gs2_coordinate_factor] == 1.0
				kit.xlabel = 'x (rho_i)'
				kit.ylabel = 'y (rho_i)'
				kit.zlabel = "'theta' rotate by 90"
			else
				kit.xlabel = 'X (m)'
				kit.ylabel = 'Y (m)'
				kit.zlabel = "'Z (m)' rotate by 90"
			end
			kit.title = "3-D Potential"
			return k

		end
	end
	def phi_real_space_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of cartesian coordinates"
		when :options
			return  [:rgbformulae, :limit, :t_index]
		else
			if options[:ncopies]
				ops = options.dup
				ops.delete(:ncopies)
				return (options[:ncopies].times.map do |n|
					ops[:ncopy] = n
					phi_real_space_graphkit(ops)
				end).sum
			end
			#zaxis = axiskit('phi0_over_x_over_y', options)
			#zaxis.data = zaxis.data.transpose
			#shape = zaxis.data.shape
			#carts = cartesian_coordinates_gsl_tensor(options)
			phi = phi_real_space_gsl_tensor(options)
			if options[:cyl]
				carts = cylindrical_coordinates_gsl_tensor(options) #.transpose(0,1,2,3)
			else
				carts = cartesian_coordinates_gsl_tensor(options)
			end

				newshape = carts.shape.slice(1..3)

			x = carts[0, false]; x.reshape!(*newshape)
			y = carts[1, false]; y.reshape!(*newshape)
			z = carts[2, false]; z.reshape!(*newshape)

			unless options[:cyl]
			#x = x.transpose(2,1,0)
			#y = y.transpose(2,1,0)
			#z = z.transpose(2,1,0)
			#phi = phi.transpose(2,1,0)
			end


			kit = GraphKit.autocreate({
				x: GraphKit::AxisKit.autocreate({data: x, title: "X", units: "m"}), 
				y: GraphKit::AxisKit.autocreate({data: y, title: "Y", units: "m"}), 
				z: GraphKit::AxisKit.autocreate({data: z, title: "Z", units: "m"}),
				f: GraphKit::AxisKit.autocreate({data: phi, title: "Phi"})
			})
# 			kit.xrange = [0,shape[0]]
# 			kit.yrange = [0,shape[1]]
			kit.cbrange = options[:limit] if options[:limit]
			kit.title = "Phi"
			#kit.palette = "rgbformulae #{(options[:rgbformulae] or "-3,3,0")}"
			#kit.view = "map"
			kit.data[0].with = "pm3d"
			#kit.gp.pm3d = "interpolate 2,2"
			#kit.title = "Phi at the outboard midplane"
# 			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def phi_gs2_space_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential as a function of y, x and theta"
		when :options
			return  [:rgbformulae, :limit, :t_index]
		else
			#zaxis = axiskit('phi0_over_x_over_y', options)
			#zaxis.data = zaxis.data.transpose
			#shape = zaxis.data.shape
			kit = GraphKit.autocreate({
				x: GraphKit::AxisKit.autocreate({data: gsl_vector('y',options), title: "y", units: "rho_#{species_letter}"}), 
				y: AxisKit.autocreate({data: gsl_vector('x',options), title: "x", units: "rho_#{species_letter}"}), 
				z: GraphKit::AxisKit.autocreate({data: gsl_vector('theta',options), title: "theta"}),
				f: GraphKit::AxisKit.autocreate({data: phi_real_space_gsl_tensor(options), title: "Phi"})
			})
# 			kit.xrange = [0,shape[0]]
# 			kit.yrange = [0,shape[1]]
			kit.cbrange = options[:limit] if options[:limit]
			kit.title = "Phi"
			#kit.palette = "rgbformulae #{(options[:rgbformulae] or "-3,3,0")}"
			#kit.view = "map"
			kit.data[0].with = "pm3d"
			#kit.title = "Phi at the outboard midplane"
# 			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def phi0_vs_x_vs_y_graphkit(options={})
		case options[:command]
		when :help
			return  "The potential at the outboard midplane"
		when :options
			return  [:rgbformulae, :limit]
		else
			zaxis = axiskit('phi0_over_x_over_y', options)
			zaxis.data = zaxis.data.transpose
			shape = zaxis.data.shape
			kit = GraphKit.autocreate({x: GraphKit::AxisKit.autocreate({data: (GSL::Vector.alloc((0...shape[0]).to_a)/shape[0] -  0.5) * (@x0 or @y0 * @jwist / 2 / Math::PI / @shat), title: "x", units: "rho_#{species_letter}"}), y: AxisKit.autocreate({data: (GSL::Vector.alloc((0...shape[1]).to_a)/shape[1] - 0.5) * @y0, title: "y", units: "rho_#{species_letter}"}), z: zaxis})
# 			kit.xrange = [0,shape[0]]
# 			kit.yrange = [0,shape[1]]
			kit.zrange = options[:limit] if options[:limit]
			kit.title = "Spectrum"
			kit.palette = "rgbformulae #{(options[:rgbformulae] or "-3,3,0")}"
			kit.view = "map"
			kit.style = "data pm3d"
			kit.title = "Phi at the outboard midplane"
# 			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.file_name = options[:graphkit_name]
			kit
		end
	end

	def growth_rate_by_mode_vs_time_graphkit(options={})
		options[:direction] = :mode
		growth_rate_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def growth_rate_by_kx_vs_time_graphkit(options={})
		options[:direction] = :kx
		growth_rate_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def growth_rate_by_ky_vs_time_graphkit(options={})
		options[:direction] = :ky
		growth_rate_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def growth_rate_by_kxy_or_mode_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "'growth_rate_by_ky_vs_time' or 'growth_rate_by_kx_vs_time': Growth rate vs time for a given kx or ky, integrated over the other direction"
		when :options
			return  [:ky, :ky_index, :kx, :kx_index]
		else
			kxy = options[:direction]
			phiax = axiskit("growth_rate_by_#{kxy}_over_time", options)	
      x = axiskit('t', options)
			x.data = x.data.subvector(0, x.data.size-1)
			kit = GraphKit.autocreate({x: x , y: phiax})
			kit.data[0].title = "Growth Rate: #{kxy} = #{options[kxy]}"	
			kit.data[0].title = "gs2:#@run_name"
			kit.data[0].with = "l" #"linespoints"
			kit.file_name = options[:graphkit_name]
			kit
		end
	end

	def es_heat_by_mode_vs_time_graphkit(options={})
		options[:direction] = :mode
		es_heat_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def es_heat_by_kx_vs_time_graphkit(options={})
		options[:direction] = :kx
		es_heat_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def es_heat_by_ky_vs_time_graphkit(options={})
		options[:direction] = :ky
		es_heat_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def es_heat_by_kxy_or_mode_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "'es_heat_by_ky_vs_time' or 'es_heat_by_kx_vs_time': Electrostatic Heat Flux vs Time for a given kx or ky, integrated over the other direction"
		when :options
			return  [:ky, :ky_index, :kx, :kx_index, :species_index]
		else
			kxy = options[:direction]
			nt_options = options.dup # 'no time' options
			#nt_options.delete(:t_index) if nt_options[:t_index]
			#nt_options.delete(:frame_index) if nt_options[:frame_index]
			phiax = axiskit("es_heat_by_#{kxy}_over_time", nt_options)	
			kit = GraphKit.autocreate({x: axiskit('t', options), y: phiax})
			kit.data[0].title = "ES Heat Flux: #@run_name #{kxy} = #{options[kxy]}"	
			kit.log_axis = 'y'
			#kit.data[0].title = "gs2:#@run_name"
			kit.data[0].with = "l" #"linespoints"
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def phi2_by_mode_vs_time_graphkit(options={})
		options[:direction] = :mode
		phi2_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def phi2_by_kx_vs_time_graphkit(options={})
		options[:direction] = :kx
		phi2_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def phi2_by_ky_vs_time_graphkit(options={})
		options[:direction] = :ky
		phi2_by_kxy_or_mode_vs_time_graphkit(options)
	end

	def phi2_by_kxy_or_mode_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "'phi2_by_ky_vs_time' or 'phi2_by_kx_vs_time': Phi^2 over time for a given kx or ky, integrated over the other direction"
		when :options
			return  [:ky, :ky_index, :kx, :kx_index]
		else
			kxy = options[:direction]
	
			# i.e. phi2_by_ky_vs_time or phi2_by_kx_vs_time or phi2_by_mode_vs_time
			
			nt_options = options.dup # 'no time' options
			nt_options.delete(:t_index) if nt_options[:t_index]
			nt_options.delete(:frame_index) if nt_options[:frame_index]
			phiax = axiskit("phi2_by_#{kxy}_over_time", nt_options)	
			kit = GraphKit.autocreate({x: axiskit('t', options), y: phiax})
			kit.data[0].title = "Phi^2 total: #{kxy} = #{options[kxy]}"	
			if options[:t_index]
# 				p 'hello'
				array_element = options[:t_index_window] ? options[:t_index] - options[:t_index_window][0] : options[:t_index] - 1
# 				p phiax.data.size, array_element
# 				p options[:t_index], options[:t_index_window]
				time = DataKit.autocreate({x: {data: GSL::Vector.alloc([list(:t)[options[:t_index]]])}, y: {data: GSL::Vector.alloc([phiax.data[array_element]]) } })
				time.pointsize = 3.0
# 				p time
# 				kit.data[0].axes[:x].data = -kit.data[0].axes[:x].data
				kit.data.push time
				if options[:with_moving_efn] and kxy==:ky
					tmax = kit.data[0].axes[:x].data[-1]
# 					p 'tmax', tmax
					
					theta_max = @g_exb * tmax * options[:ky] * 2 * Math::PI / list(:kx)[2]
					kit.each_axiskit(:x) do |axiskit|
# 						p axiskit
						axiskit.data = axiskit.data / tmax * theta_max - theta_max
					end
				end
			end
			if options[:norm]
				xrange, yrange = kit.plot_area_size
				kit.each_axiskit(:y) do |axiskit|
					axiskit.data /= yrange[1] / (options[:height] or 1.0)
				end
			end
			kit.log_axis = 'y'
			#kit.data[0].title = "gs2:#@run_name"
			kit.data[0].with = "l" #"linespoints"
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def apar2_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "Graph of apar^2 vs time integrated over all space. No options"
		when :options
			return []
		else
			kit = GraphKit.autocreate({x: axiskit('t'), y: axiskit('apar2_over_time', options)})
			kit.data[0].with = "l" #"linespoints"
			kit.data[0].title = "Apar^2 total:#@run_name"	
			kit.file_name = options[:graphkit_name]
			kit.log_axis = 'y'
			kit
		end
	end
	def phi2tot_vs_time_graphkit(options={})
		case options[:command]
		when :help
			return  "Graph of phi^2 vs time integrated over all space. No options"
		when :options
			return []
		else
			kit = GraphKit.autocreate({x: axiskit('t'), y: axiskit('phi2tot_over_time', options)})
			kit.data[0].with = "l" #"linespoints"
			kit.data[0].title = "Phi^2 total:#@run_name"	
			kit.file_name = options[:graphkit_name]
			kit.log_axis = 'y'
			kit
		end
	end
	def spectrum_graphkit(options={})
		case options[:command]
		when :help
			return "Graph of phi^2 at a given time vs kx and ky"
		when :options
			return  [:with]
		else
# 	                        p @name_match 
			#options[:times_kx4] = true if @name_match[:kxsq]
			zaxis = axiskit('spectrum_over_ky_over_kx', options)
			zaxis.data = zaxis.data.transpose
			kit = GraphKit.autocreate({y: axiskit('ky', options), x: axiskit('kx', options), z: zaxis})
			kit.title = "#{options[:log] ?  "Log ": ""}Spectrum (phi^2#{options[:times_kx2] ? " * kx^2" : ""}#{options[:times_kx4] ? " * kx^4" : ""})"
			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def spectrum_vs_kpar_vs_ky_graphkit(options={})
		case options[:command]
		when :help
			return "Graph of phi^2 * ky^2 at a given time vs kpar and ky"
		when :options
			return  [:with, :log, :no_zonal, :no_kpar0]
		else
# 	                        p @name_match 
			#options[:times_kx4] = true if @name_match[:kxsq]
			zaxis = axiskit('spectrum_over_ky_over_kpar', options)
			zaxis.data = zaxis.data.transpose
			kit = GraphKit.autocreate({y: axiskit('ky', options), x: axiskit('kpar', options.modify({ky_index: 1, kx_index: 1})), z: zaxis})
			kit.title = "#{options[:log] ?  "Log ": ""}Spectrum (phi^2 ky^2)"
			kit.data[0].with = (options[:with] or 'pm3d palette')
			kit.gp.view = "map" if options[:map]
			kit.file_name = options[:graphkit_name]
			kit
		end
	end
	def vspace_diagnostics_graphkit(options={})
		case options[:command]
		when :help
			return "Plots vspace diagnostics. All lines shouldn't stray much above 0.1 - otherwise large amounts of the distribution function is in the higher k velocity space and velocity space is probably unresolved. (NB This graph is here temporarily (ha ha) until I add the vspace diagnostics to the NetCDF file (or the apocalypse, whichever is sooner) EGH)"
		when :options
			return []
		else
			raise "Velocity space diagnostics not found" unless FileTest.exist? "#@directory/#@run_name.lpc" or FileTest.exist? "#@directory/#@run_name.vres"
			lpc = GSL::Vector.filescan("#@directory/#@run_name.lpc") rescue [[0], [0], [0]]
			vres = GSL::Vector.filescan("#@directory/#@run_name.vres") rescue [[0], [0], [0]]
			xaxis = AxisKit.autocreate({data: lpc[0], title: "Time"})
			data = [[lpc[0], lpc[1], "Pitch Angle Harmonics (lpc)"], [lpc[0], lpc[2], "Energy Harmonics (lpc)"], [vres[0], vres[1], "Pitch Angle Harmonics (vres)"], [vres[0], vres[2], "Energy Harmonics (vres)"]]
			kits = data.inject([]) do |arr, (x, vector, title)|
				arr + [GraphKit.autocreate({x: AxisKit.autocreate({data: x, title: "Time"}), y: AxisKit.autocreate({data: vector, title: title})})]
			end
			kit = kits.sum
# 				exit
			kit.title = "Velocity Space Diagnostics"
			kit.ylabel = "Fraction of Dist Func Contained"
			kit.file_name = options[:graphkit_name]
# 			kit.log_axis = 'y'
			kit
		end
	end
	def zonal_spectrum_graphkit(options={})
		case options[:command]
		when :help
			return  "zonal_spectrum: Graph of kx^4 phi^2 vs kx for ky=0"
		when :options
			return  [:t, :t_index]
		else
			options[:times_kx4] = true
			kit = GraphKit.autocreate({x: axiskit('kx', options), y: axiskit("zonal_spectrum", options)})
			kit.title  = "Zonal Spectrum"
			kit.file_name = options[:graphkit_name] + options[:t_index].to_s
			kit.data[0].with = 'lp'
			kit.pointsize = 2.0
			kit
		end
	end


end

include GraphKits

end
end
