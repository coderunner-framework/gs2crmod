##########################################
# Calculations for GS2 Code Runner Module
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
class Gs2 



def calculate_time_averaged_fluxes
	eputs 'Calculating time averaged fluxes'
	calculate_saturation_time_index unless @saturation_time_index
	return unless FileTest.exist?(netcdf_filename)
	@hflux_tot_stav = saturated_time_average('hflux_tot_over_time', {})
	@hflux_tot_stav_error = saturated_time_average_error('hflux_tot_over_time', {})
	@hflux_tot_stav_std_dev = saturated_time_average_std_dev('hflux_tot_over_time', {})
	@phi2_tot_stav = saturated_time_average('phi2tot_over_time', {})
	#@par_mom_flux_stav = saturated_time_average('par_mom_flux_over_time', {}) rescue nil
	#@perp_mom_flux_stav = saturated_time_average('perp_mom_flux_over_time', {}) rescue nil
	@es_mom_flux_stav = {}
	@es_heat_flux_stav = {}
	@es_mom_flux_stav_error = {}
	@es_mom_flux_stav_std_dev = {}
	@es_heat_flux_stav_error = {}
	@es_heat_flux_stav_std_dev = {}

	@nspec.times do |i|
		species_index = i + 1
		@es_mom_flux_stav[species_index]  = saturated_time_average('es_mom_flux_over_time', {species_index: species_index})
		@es_heat_flux_stav[species_index]  = saturated_time_average('es_heat_flux_over_time', {species_index: species_index})
		@es_mom_flux_stav_error[species_index]  = saturated_time_average_error('es_mom_flux_over_time', {species_index: species_index})
		@es_mom_flux_stav_std_dev[species_index]  = saturated_time_average_std_dev('es_heat_flux_over_time', {species_index: species_index})
		@es_heat_flux_stav_error[species_index]  = saturated_time_average_error('es_heat_flux_over_time', {species_index: species_index})
		@es_heat_flux_stav_std_dev[species_index]  = saturated_time_average_std_dev('es_heat_flux_over_time', {species_index: species_index})
	end
# 	ep @es_mom_flux_stav, @es_heat_flux_stav
end

alias :ctaf :calculate_time_averaged_fluxes

def saturated_time_average(name, options)
# 	calculate_saturation_time_index unless @saturation_time_index
# 	p 'sat', @saturation_time_index, 'max', list(:t).keys.max
	raise "saturation_time_index not calculated for #@run_name" unless @saturation_time_index
	options[:t_index_window] = [@saturation_time_index, list(:t).keys.max - 1]
	#ep gsl_vector(name, {}).size
	#ep name, options
	begin
		vec = gsl_vector(name, options)
	rescue GSL::ERROR::EINVAL
		# IF the vector doesn't have enough values for each timestep (due to run aborting early?), this error will be thrown.
		options[:t_index_window] = [@saturation_time_index, gsl_vector(name, {}).size]
		retry
	rescue NoMethodError
		eputs "Warning: could not calculate #{name} saturated time average"
		return nil
	end
	
		                                                               
		tvec = gsl_vector('t', options)

		                                                               
	dt = tvec.subvector(1, tvec.size - 1) - tvec.subvector(0, tvec.size - 1)
	trapezium = (vec.subvector(1, tvec.size - 1) + vec.subvector(0, tvec.size - 1)) / 2.0
	return trapezium.mul(dt).sum / dt.sum
end

def saturated_time_average_error(name, options)
# 	calculate_saturation_time_index unless @saturation_time_index
	options[:t_index_window] = [@saturation_time_index, list(:t).keys.max]
	begin
		vec = gsl_vector(name, options)
		tavg = GSL::Vector.alloc(vec.size)
		vec.size.times.each{|i| tavg[i] = vec.subvector(i+1).mean}
	rescue NoMethodError
		eputs "Warning: could not calculate #{name} saturated_time_average_error"
		return nil
	end
# 	tavg = 0.0; i = 0
	
# 	tavg_vec = vec.collect{|val| tavg += val; tavg = tavg / (i+=1); tavg}
# 	ind = GSL::Vector.indgen(vec.size)
# 	i = 0
# 	begin 
# 		fit = GSL::Fit::linear(ind.subvector(i, ind.size - i) , vec.subvector(i, ind.size - i))
# # 		p fit[1].abs - 100.0 * fit[4].abs
# 		i += 1
# 		(eputs "Not Saturated"; break) if i > vec.size * 0.9
# 	end while (fit[1].abs - Math.sqrt(fit[4].abs)) > 0 
# 	p fit
# 	fit_vec = ind * fit[1] + fit[0]
# # 	p tavg.size
# 	# GraphKit.autocreate({x: {data: gsl_vector(name, {})}})
# 	(GraphKit.autocreate({x: {data: tavg}}) + GraphKit.autocreate({x: {data: vec}}) + GraphKit.autocreate({x: {data: fit_vec}})).gnuplot
	return tavg.sd
end

def saturated_time_average_std_dev(name, options)
# 	calculate_saturation_time_index unless @saturation_time_index
	options[:t_index_window] = [@saturation_time_index, list(:t).keys.max]
	begin
		vec = gsl_vector(name, options)
	rescue NoMethodError
		eputs "Warning: could not calculate #{name} saturated_time_average_std_dev"
		return nil
	end
	return vec.sd
end


# I.e. the time at which the primary modes are saturated and the fluxes settle around a long term average. 

def calculate_saturation_time_index(show_graph = false)
	
	eprint "Checking for saturation..."

	#hflux = gsl_vector('hflux_tot_over_time', {})
	hflux = gsl_vector('phi2tot_over_time', {})
	
 	#eputs 'got hflux'
	#ep 'hflux', hflux
	
	#Check if it's decayed to 0
	if hflux[-1] < 1.0e-10
		for i in 1..hflux.size
# 			raise "negative heat flux: #{hflux[-i]} " if hflux[-i] < 0
			(break) unless hflux[- i] < 1.0e-10
		end
		if i > hflux.size * 1.0/10.0 #i.e if was 0 for more than a tenth of the time
			@saturated = true
			@saturation_time_index = hflux.size - i + 1
			eputs "saturation time = #{list(:t)[@saturation_time_index]}"
			GraphKit.quick_create([gsl_vector('t',{}), hflux]).gnuplot(log_axis: 'y') if show_graph
			return
		end
	end
		
	# Get initial estimate for saturation time
	for i in 0...hflux.size
		rem = hflux.subvector(i, hflux.size - i)
		break if (hflux[i] - rem.mean).abs < rem.sd / 2.0
		break if i > 3.0/4.0*hflux.size
	end
	
	@saturation_time_index = [i + 1, hflux.size - 2].min
	
# 	fit = GSL::Fit::linear(GSL::Vector.indgen(rem.size), rem)
# 	
# 	slope, covar11 = fit[1], fit[4]
# 	range = [slope + Math.sqrt(covar11), slope - Math.sqrt(covar11)]
# 	
# 	unless range.min < 0 and range.max > 0
# 		eputs "Warning: This run (#{id}) has probably not reached a saturated state: the estimated slope of the heat flux is in this range: #{range.inspect}"
# 		@saturated = false
# 	end
# 	
# 	ep fit
	
# 	eputs "Saturation time estimate', @saturation_time_index = i + 1
# 	t_vec[@saturation_time_index - 1]
	max_t_index = list(:t).keys.max
	max_t = list(:t).values.max
	min_t = list(:t).values.min
	#hflux = gsl_vector('hflux_tot_over_time', {:t_index_window => [@saturation_time_index, max_t_index]})
	hflux = gsl_vector('phi2tot_over_time', {:t_index_window => [@saturation_time_index, max_t_index]})
	t_vec = gsl_vector('t', {:t_index_window => [@saturation_time_index, max_t_index]})
# 	p t_vec[0]
	i = 0
	t_arr = []; conf_arr = []
	loop do
		eprint '.'
		
# 		GraphKit.autocreate(x: {data: t_vec}, y: {data: hflux}).gnuplot
		
		lomb = GSL::SpectralAnalysis::Lomb.alloc(t_vec.subvector(i, t_vec.size - i),  hflux.subvector(i, hflux.size - i))
		fs, periodogram = lomb.calculate_periodogram(1.0, 4.0, [0]) #(1.0) #0.1 * hflux.size / ( hflux.size - i))
# 		lomb.graphkit.gnuplot
		
# 		eputs 'Confidence that lowest frequency is not noise is: '
		# pnoise is the probability of the strength of the lowest frequency signal in the heat flux given a hypothesis of gaussian noise. If it is high there is a low likelihood that there is a signal at the lowest frequency: ie. within that window the heat flux has reached a stationary state
		pnoise = lomb.pnull(periodogram[0])
		t_arr.push t_vec[i]; conf_arr.push pnoise
		
		(@saturated = true; break) if pnoise > 0.9
		step = (hflux.size / 25.0).to_i
		step = 1 if step==0
		i += step
		#(@saturated = false; i ; break) if (i >= t_vec.size or t_vec[i] > (max_t - min_t) * 2.0 / 3.0 + min_t )
		(@saturated = false; break) if (i >= t_vec.size or t_vec[i] > (max_t - min_t) * 2.0 / 3.0 + min_t )
		@saturation_time_index += step	
#		ep '---i,t,size',i, t_vec[i], t_vec.size
	end
	(kit = GraphKit.autocreate({x: {data: t_vec}, y: {data: hflux / hflux.max}}, {x: {data: t_arr}, y: {data: conf_arr}}); kit.data[1].with = 'lp'; kit.gnuplot) if show_graph #(log_axis: 'y')
# 	puts 
	if @saturated
# 		p i
		eputs "saturation time = #{list(:t)[@saturation_time_index]}"
	else
		eputs "run not saturated"
	end
		
	return
	exit
	# Get regularly spaced t vector
	
# 	
# 	t_delta_vec = GSL::Vector.alloc(t_vec.size - 1)
# 	t_delta_vec.size.times.each{|i| t_delta_vec[i] = t_vec[i+1] - t_vec[i]}
# 	
# 	ep t_delta_vec.max, t_delta_vec.min
# 	
# 	even_t = GSL::Vector.linspace(t_vec.min, t_vec.max, ((t_vec.max - t_vec.min) / t_delta_vec.max).round )
# 	
# # 	even_t = []
# # 	tm = t = t_vec[t_delta_vec.max_index]
# 	
# # 	loop do
# # 		even_t.push t
# 		
# # 	
# 	ep even_t.size, t_vec.size
# 	
# 	min_delt = t_delta_vec.min
# 	p even_t.any?{|el| bool = (not t_vec.any?{|ele| (ele - el).abs < 1.0e-1 * min_delt}); ep el if bool; bool}
# 	
# 	ep t_vec.dup.delete_if{|el| not (el - 71.3).abs < 0.5}
# 	
# 	exit
	
	
	
	
	return
	
	# Calculate a series of time averaged segments
	pieces = hflux.pieces(20) # split into 20 pieces
	avgs = GSL::Vector.alloc(pieces.map{|vec| vec.sum/vec.size})
	# Calculate their variance
	mean = (avgs.sum/avgs.size)
	sig = Math.sqrt((avgs.square - mean**2).sum/avgs.size)
	# Discount any at the start which are more than one standard deviation away from the average - they are from the linear growth phase
	t_index = 1
	kept_avgs = avgs.dup
	for i in 0...pieces.size
		if (avgs[i] - mean).abs > sig
			kept_avgs.delete_at(i)
			t_index += pieces[i].size
		else
			break
		end
	end
	eputs "Warning: probably not saturated" if [kept_avgs, kept_avgs.reverse].include? kept_avgs.sort
	ep kept_avgs
	@saturation_time_index = t_index
# 	p t_index, list(:t)[t_index]
end

alias :csti :calculate_saturation_time_index

#Actually, this doesn't calculate the frequencies but reads them from run_name.out. Requires write_line to be .true.
#
def calculate_frequencies
		@real_frequencies = FloatHash.new
		gs2_out = FileUtils.tail(@run_name + ".out", list(:ky).size*list(:kx).size)
# 		a  = gs2_out.split("\n")
		final_timestep_list = gs2_out #a.slice((a.size-@ky_list.size*@kx_list.size-1)..a.size-1).join("\n")
 		log(final_timestep_list.slice(-2..-1))
# 		eputs final_timestep_list
		f = LongRegexen::FLOAT.verbatim
		logi(f)
		@frequency_at_ky_at_kx||= FloatHash.new
		ky_values = []
		regex = Regexp.new( "^.*aky=\\s*(?<aky>#{f})\s*akx=\\s*(?<akx>#{f}).*omav=\\s*(?<re>#{f})\\s*(?<gr>#{f})")
		final_timestep_list.scan(regex) do
			aky = eval($~[:aky])
			akx = eval($~[:akx])
			@frequency_at_ky_at_kx[aky] = FloatHash.new unless ky_values.include? aky
			ky_values.push aky
			@frequency_at_ky_at_kx[aky][akx] = eval($~[:re])
		end
end
def calculate_growth_rates_and_frequencies
        return if @grid_option == "single" and @aky == 0.0 # no meaningful results
	Dir.chdir(@directory) do
		logf(:calculate_growth_rates_and_frequencies)
		logd

		calculate_frequencies
		
# 		get_list_of(:ky, :kx)
		@growth_rates= FloatHash.new
			#raise CRFatal.new("Unknown value of ky read from output file: #{data[:aky].to_f}. Not in list:\n#{list(:ky).values.inspect}") 
# 		pp @ky_list
		
		# With zero magnetic shear, calculate growth rates for both kx and ky
		#if @shat and @shat.abs < 1.0e-5 and @nx and @nx > 1 
			to_calc = [:kx, :ky]
			@growth_rate_at_kx ||= FloatHash.new
		#else
			#to_calc = [:ky]
		#end
		
		@growth_rate_at_ky ||= FloatHash.new
 		eputs
#		p @growth_rate_at_kx; exit
		to_calc.each do |kxy|
			growth_rates = send(:growth_rate_at_ + kxy)
		list(kxy).values.sort.each do |value|
			
			#p growth_rates.keys, value, growth_rates[value.to_f-0.0],
			#growth_rates.class, growth_rates.keys.include?(value); exit
	
			next if growth_rates.keys.include? value

			
			Terminal.erewind(1)
			#ep growth_rates.keys
			eputs sprintf("Calculating growth rate for #{kxy} = % 1.5e#{Terminal::CLEAR_LINE}", value) 
			

					# Mode has 0 growth rate at ky==0
			(growth_rates[value] = 0.0; next) if value == 0.0 and kxy == :ky 
			if @g_exb_start_timestep
				t_index_window = [1, [(g_exb_start_timestep-1)/@nwrite, list(:t).keys.max].min]
				#ep "t_index_window", t_index_window
			else
				t_index_window = nil
			end
			if list(kxy).size == 1
				phi2_vec = gsl_vector("phi2tot_over_time", t_index_window: t_index_window)
			else
				phi2_vec = gsl_vector("phi2_by_#{kxy}_over_time", kxy=>value, :t_index_window=> t_index_window)
			end
			(growth_rates[value] = 0.0; next) if phi2_vec.min <= 0.0
			growth_rates[value] = calculate_growth_rate(phi2_vec)
			(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rate\n----------\n\n"; growth_rates[value] = -1; next) if growth_rates[value] == "NaN"
		end
		end
		
 		write_results
		
# 		ep "growth_rate_at_ky", @growth_rate_at_ky
		if ENV['GS2_CALCULATE_ALL']
		trap(0){eputs "Calculation of spectrum did not complete: run 'cgrf' (i.e. calculate_growth_rates_and_frequencies) for this run. E.g. from the command line \n $ coderunner rc 'cgrf' -j #{@id}"; exit}
		@growth_rate_at_ky_at_kx ||= FloatHash.new
		list(:ky).values.sort.each do |kyv|
			# MJL 2013-11-07: The line below originally used ||= instead of =. I'm not sure why, since ||= does not seem to work.
			@growth_rate_at_ky_at_kx[kyv] = FloatHash.new
			list(:kx).values.sort.each do |kxv|	
				# MJL 2013-11-07: I'm not sure why this next line was originally included. It seemed to cause almost all k's to be skipped.
				#next if @growth_rate_at_ky_at_kx[kyv].keys.include? kxv
				Terminal.erewind(1)
				eputs sprintf("Calculating growth rate for kx = % 1.5e and ky = % 1.5e#{Terminal::CLEAR_LINE}", kxv, kyv) 
				(@growth_rate_at_ky_at_kx[kyv][kxv] = 0.0; next) if kyv == 0.0 # Mode has 0 growth rate at ky==0
				phi2_vec = gsl_vector("phi2_by_mode_over_time", {:kx=>kxv, :ky=>kyv})
				(@growth_rate_at_ky_at_kx[kyv][kxv] = 0.0; next) if phi2_vec.min <= 0.0
				@growth_rate_at_ky_at_kx[kyv][kxv] = calculate_growth_rate(phi2_vec)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rates\n----------\n\n"; @growth_rate_at_ky_at_kx[kyv][kxv] = -1; next) if @growth_rate_at_ky_at_kx[kyv][kxv] == "NaN" 
			end
			write_results
		end
		trap(0){}
		end
		@growth_rates = @growth_rate_at_ky
		@max_growth_rate = @growth_rates.values.max
		@fastest_growing_mode = @growth_rates.key(@max_growth_rate)
    @freq_of_max_growth_rate = @real_frequencies[@fastest_growing_mode] rescue nil
		ep @max_growth_rate, @growth_rates
		@decaying = (@max_growth_rate < 0) if @max_growth_rate
		@ky = @aky if @aky
		if @grid_option == "single"
# 			ep @aky, @growth_rates
			@gamma_r = @growth_rates[@aky.to_f]
			@gamma_i = @real_frequencies[@aky.to_f]
		end
# 		ep @gamma_r
		
		
# 		eputs @growth_rates; gets
	end
end

alias :cgrf :calculate_growth_rates_and_frequencies

def calculate_growth_rate(vector, options={})
	raise "This vector should be positive definite" if vector.min < 0.0
	offset = 0
	length = vector.length
  while vector[offset] == 0.0
    offset+=1
    return 0.0 if offset == vector.length
  end
	growth_rate = GSL::Fit::linear(gsl_vector(:t).subvector(offset, length-offset), 0.5*GSL::Sf::log(vector.subvector(offset, length - offset)))[1]
	divisor = 1
	while (growth_rate.to_s == "NaN")
			#This corrects the growth rate if phi has grown all the way to NaN during the simulation
		divisor *= 2
		length = (vector.size.to_f / divisor.to_f).floor
# 				p length
		return "NaN" if length <= offset + 1
		growth_rate = GSL::Fit::linear(gsl_vector(:t).subvector(offset, length-offset), 0.5*GSL::Sf::log(vector.subvector(offset, length-offset)))[1]
	end	
	growth_rate
end

# Not needed for GS2 built after 16/06/2010

def corrected_mom_flux_stav
	par_mom_flux_stav - perp_mom_flux_stav
end

def calculate_transient_amplifications
  return if @grid_option == "single" and @aky == 0.0 # no meaningful results
	Dir.chdir(@directory) do
		# With zero magnetic shear, calculate amplifications for both kx and ky
		if @shat and @shat.abs < 1.0e-5 and @nx > 1 
			to_calc = [:kx, :ky]
			@transient_amplification_at_kx ||= FloatHash.new
		else
			to_calc = [:ky]
		end
		
		@transient_amplification_at_ky ||= FloatHash.new
 		eputs
		to_calc.each do |kxy|
			transient_amplifications = send(:transient_amplification_at_ + kxy)
			list(kxy).values.sort.each do |value|
			
				#p transient_amplifications.keys, value, transient_amplifications[value.to_f-0.0],
				#transient_amplifications.class, transient_amplifications.keys.include?(value); exit
		
				next if transient_amplifications.keys.include? value

				
				Terminal.erewind(1)
				#ep transient_amplifications.keys
				eputs sprintf("Calculating transient amplification for #{kxy} = % 1.5e#{Terminal::CLEAR_LINE}", value) 
				

						# Mode has 0 growth rate at ky==0
				(transient_amplifications[value] = 0.0; next) if value == 0.0 and kxy == :ky 
				phi2_vec = gsl_vector("phi2_by_#{kxy}_over_time", {kxy=>value})
				#(transient_amplifications[value] = 0.0; next) if phi2_vec.min <= 0.0
				transient_amplifications[value] = calculate_transient_amplification(phi2_vec)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rate\n----------\n\n"; transient_amplifications[value] = -1; next) if transient_amplifications[value].to_s == "NaN"
			end
		end
		
 		write_results
		
# 		ep "transient_amplification_at_ky", @transient_amplification_at_ky
		if ENV['GS2_CALCULATE_ALL']
		trap(0){eputs "Calculation of spectrum did not complete: run 'cgrf' (i.e. calculate_transient_amplifications_and_frequencies) for this run. E.g. from the command line \n $ coderunner rc 'cgrf' -j #{@id}"; exit}
		@transient_amplification_at_ky_at_kx ||= FloatHash.new
		list(:ky).values.sort.each do |kyv|
			@transient_amplification_at_ky_at_kx[kyv] ||= FloatHash.new
			#p @transient_amplification_at_ky_at_kx[kyv]
			list(:kx).values.sort.each do |kxv|	
				next if @transient_amplification_at_ky_at_kx[kyv].keys.include? kxv
				Terminal.erewind(1)
				eputs sprintf("Calculating growth rate for kx = % 1.5e and ky = % 1.5e#{Terminal::CLEAR_LINE}", kxv, kyv) 
				(@transient_amplification_at_ky_at_kx[kyv][kxv] = 0.0; next) if kyv == 0.0 # Mode has 0 growth rate at ky==0
				phi2_vec = gsl_vector("phi2_by_mode_over_time", {:kx=>kxv, :ky=>kyv})
				#(@transient_amplification_at_ky_at_kx[kyv][kxv] = 0.0; next) if phi2_vec.min <= 0.0
				@transient_amplification_at_ky_at_kx[kyv][kxv] = calculate_transient_amplification(phi2_vec)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rates\n----------\n\n"; @transient_amplification_at_ky_at_kx[kyv][kxv] = -1; next) if @transient_amplification_at_ky_at_kx[kyv][kxv].to_s == "NaN" 
			end
			write_results
		end
		trap(0){}
		end
		@transient_amplifications = @transient_amplification_at_ky
		@max_transient_amplification = @transient_amplifications.values.max
		@most_amplified_mode = @transient_amplifications.key(@max_transient_amplification)
		#@freq_of_max_transient_amplification = @real_frequencies[@fastest_growing_mode]
		#ep @max_transient_amplification, @transient_amplifications
		#@decaying = (@max_transient_amplification < 0) if @max_transient_amplification
		@ky = @aky if @aky
		#if @grid_option == "single"
## 			ep @aky, @transient_amplifications
			#@gamma_r = @transient_amplifications[@aky.to_f]
			#@gamma_i = @real_frequencies[@aky.to_f]
		#end
# 		ep @gamma_r
		
		
# 		eputs @transient_amplifications; gets
	end
end

alias :cta :calculate_transient_amplifications


def calculate_transient_es_heat_flux_amplifications
  return if @grid_option == "single" and @aky == 0.0 # no meaningful results

	@transient_es_heat_flux_amplification_at_species_at_kx = []
	@transient_es_heat_flux_amplification_at_species_at_ky = []
  @transient_es_heat_flux_amplification_at_species_at_ky_at_kx = []
	for species_index in 1..nspec

	Dir.chdir(@directory) do
		# With zero magnetic shear, calculate amplifications for both kx and ky
		if @shat and @shat.abs < 1.0e-5 and @nx > 1 and !@ikx_init and false
			to_calc = [:kx, :ky]
			@transient_es_heat_flux_amplification_at_species_at_kx[species_index-1] ||= FloatHash.new
		else
			to_calc = [:ky]
		end
		
		@transient_es_heat_flux_amplification_at_species_at_ky[species_index-1] ||= FloatHash.new
 		eputs
		to_calc.each do |kxy|
			transient_es_heat_flux_amplifications = send(:transient_es_heat_flux_amplification_at_species_at_ + kxy)[species_index-1]
			list(kxy).values.sort.each do |value|
			
				#p transient_es_heat_flux_amplifications.keys, value, transient_es_heat_flux_amplifications[value.to_f-0.0],
				#transient_es_heat_flux_amplifications.class, transient_es_heat_flux_amplifications.keys.include?(value); exit
		
				next if transient_es_heat_flux_amplifications.keys.include? value

				
				Terminal.erewind(1)
				#ep transient_es_heat_flux_amplifications.keys
				eputs sprintf("Calculating transient amplification for #{kxy} = % 1.5e#{Terminal::CLEAR_LINE}", value) 
				

						# Mode has 0 growth rate at ky==0
				(transient_es_heat_flux_amplifications[value] = 0.0; next) if value == 0.0 and kxy == :ky 
				phi2_vec = gsl_vector("es_heat_by_#{kxy}_over_time", {kxy=>value, species_index: species_index})
				#(transient_es_heat_flux_amplifications[value] = 0.0; next) if phi2_vec.min <= 0.0
				transient_es_heat_flux_amplifications[value] = calculate_transient_amplification(phi2_vec)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rate\n----------\n\n"; transient_es_heat_flux_amplifications[value] = -1; next) if transient_es_heat_flux_amplifications[value].to_s == "NaN"
			end
		end
		
 		write_results
		
# 		ep "transient_es_heat_flux_amplification_at_species_at_ky", @transient_es_heat_flux_amplification_at_species_at_ky
		if ENV['GS2_CALCULATE_ALL']
		trap(0){eputs "Calculation of spectrum did not complete: run 'ctehfa' (i.e. calculate_transient_es_heat_flux_amplifications) for this run. E.g. from the command line \n $ coderunner rc 'ctehfa' -j #{@id}"; exit}
		@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1] ||= FloatHash.new
		list(:ky).values.sort.each do |kyv|
			@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv] ||= FloatHash.new
			#p @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[kyv]
			list(:kx).values.sort.each do |kxv|	
				next if @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv].keys.include? kxv
				Terminal.erewind(1)
				eputs sprintf("Calculating growth rate for kx = % 1.5e and ky = % 1.5e#{Terminal::CLEAR_LINE}", kxv, kyv) 
				(@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv] = 0.0; next) if kyv == 0.0 # Mode has 0 growth rate at ky==0
				phi2_vec = gsl_vector("phi2_by_mode_over_time", {:kx=>kxv, :ky=>kyv})
				#(@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[kyv][kxv] = 0.0; next) if phi2_vec.min <= 0.0
				@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv] = calculate_transient_es_heat_flux_amplification(phi2_vec)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rates\n----------\n\n"; @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv] = -1; next) if @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv].to_s == "NaN" 
			end
			write_results
		end
		trap(0){}
		end
		#@max_transient_es_heat_flux_amplification = @transient_es_heat_flux_amplifications.values.max
		#@most_amplified_mode = @transient_es_heat_flux_amplifications.key(@max_transient_es_heat_flux_amplification)
		#@ky = @aky if @aky
	end
	end # for species_index in 1..nspec
end

def calculate_transient_es_heat_flux_amplifications
  return if @grid_option == "single" and @aky == 0.0 # no meaningful results
	@transient_es_heat_flux_amplification_at_species_at_kx = []
	@transient_es_heat_flux_amplification_at_species_at_ky = []
  @transient_es_heat_flux_amplification_at_species_at_ky_at_kx = []
	for species_index in 1..nspec

	Dir.chdir(@directory) do
		#trap(0){eputs "Calculation of transient amplification did not complete: run 'ctehfa' (i.e. calculate_transient_es_heat_flux_amplifications) for this run. E.g. from the command line \n $ coderunner rc 'ctehfa' -j #{@id}"; exit}
		es_heat_flux_narray = netcdf_file.var('es_heat_by_k').get

		ep es_heat_flux_narray.shape, '', '', ''
		t_size = es_heat_flux_narray.shape[-1]
		temp_es_heat = GSL::Vector.alloc(t_size)
		@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1] ||= FloatHash.new
		ky_index=0
		list(:ky).values.each do |kyv|
			ky_index +=1
			kx_index = 0
			@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv] ||= FloatHash.new
			#p @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[kyv]
			list(:kx).values.each do |kxv|	
				kx_index +=1
				next if @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv].keys.include? kxv
				Terminal.erewind(1)
				eputs sprintf("Calculating transient amplification for kx = % 1.5e and ky = % 1.5e#{Terminal::CLEAR_LINE}", kxv, kyv) 
				(@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv] = 0.0; next) if kyv == 0.0 # Mode has 0 growth rate at ky==0
				#phi2_vec = gsl_vector("phi2_by_mode_over_time", {:kx=>kxv, :ky=>kyv})

				#ep({kx: kxv})
				#hash = ({kx: kxv})
				#ep(hash.convert_to_index(self, :kx))
				#kx_index = {kx: kxv}.convert_to_index(self, :kx)[:kx_index]
				#ky_index = {ky: kyv}.convert_to_index(self, :ky)[:ky_index]
				#p t_size
				for i in 0...t_size
					#p i
					temp_es_heat[i] = es_heat_flux_narray[kx_index-1, ky_index -1, species_index-1, i]
					#temp_es_heat[i] = es_heat_flux_narray[i, species_index-1, ky_index-1, kx_index-1]
				end
				#(@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[kyv][kxv] = 0.0; next) if phi2_vec.min <= 0.0
				@transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv] = calculate_transient_amplification(temp_es_heat)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rates\n----------\n\n"; @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv] = -1; next) if @transient_es_heat_flux_amplification_at_species_at_ky_at_kx[species_index-1][kyv][kxv].to_s == "NaN" 
			end
			write_results
		end
		#trap(0){}
	# With zero magnetic shear, calculate amplifications for both kx and ky

	  #return
		if @shat and @shat.abs < 1.0e-5 and @nx > 1 and !@ikx_init and false
			to_calc = [:kx, :ky]
			@transient_es_heat_flux_amplification_at_species_at_kx[species_index-1] ||= FloatHash.new
		else
			to_calc = [:ky]
		end
		
		@transient_es_heat_flux_amplification_at_species_at_ky[species_index-1] ||= FloatHash.new
 		eputs
		to_calc.each do |kxy|
			transient_es_heat_flux_amplifications = send(:transient_es_heat_flux_amplification_at_species_at_ + kxy)[species_index-1]
			list(kxy).values.sort.each do |value|
			
				#p transient_es_heat_flux_amplifications.keys, value, transient_es_heat_flux_amplifications[value.to_f-0.0],
				#transient_es_heat_flux_amplifications.class, transient_es_heat_flux_amplifications.keys.include?(value); exit
		
				next if transient_es_heat_flux_amplifications.keys.include? value

				
				Terminal.erewind(1)
				#ep transient_es_heat_flux_amplifications.keys
				eputs sprintf("Calculating transient amplification for #{kxy} = % 1.5e#{Terminal::CLEAR_LINE}", value) 
				

						# Mode has 0 growth rate at ky==0
				(transient_es_heat_flux_amplifications[value] = 0.0; next) if value == 0.0 and kxy == :ky 
				phi2_vec = gsl_vector("es_heat_by_#{kxy}_over_time", {kxy=>value, species_index: species_index})
				#(transient_es_heat_flux_amplifications[value] = 0.0; next) if phi2_vec.min <= 0.0
				transient_es_heat_flux_amplifications[value] = calculate_transient_amplification(phi2_vec)
				(eputs "\n\n----------\nIn #@run_name:\n\nphi2_by_#{kxy}_over_time is all NaN; unable to calculate growth rate\n----------\n\n"; transient_es_heat_flux_amplifications[value] = -1; next) if transient_es_heat_flux_amplifications[value].to_s == "NaN"
			end
		end
		
 		write_results
		
# 		ep "transient_es_heat_flux_amplification_at_species_at_ky", @transient_es_heat_flux_amplification_at_species_at_ky

		#@max_transient_es_heat_flux_amplification = @transient_es_heat_flux_amplifications.values.max
		#@most_amplified_mode = @transient_es_heat_flux_amplifications.key(@max_transient_es_heat_flux_amplification)
		#@ky = @aky if @aky
	end # Dir.chdir
	end # for species_index in 1..nspec
end

alias :ctehfa :calculate_transient_es_heat_flux_amplifications
alias :ctehfa :calculate_transient_es_heat_flux_amplifications


def calculate_transient_amplification(vector, options={})
  if @g_exb_start_timestep 
    return GSL::Sf::log(vector[(@g_exb_start_timestep/@nwrite).to_i...-1].max / vector[(@g_exb_start_timestep/@nwrite).to_i])/2
  else
    eputs "Warning: set g_exb_start_timestep to calculate transient amplifications."
    return 0
  end
end

def ctan
	list(:ky).each do |(ky_index, ky)|
		eputs "ky: #{ky}"
		phi_vec = gsl_vector("phi2_by_ky_over_time", ky_index: ky_index)
		t_element	= 0
		old = phi_vec[0]

	  loop do 
			t_element+=1
			#print t_element, ',', phi_vec.size
			new = phi_vec[t_element]
			break if new > old or t_element == phi_vec.size - 1
			old = new
		end
		
		if t_element == phi_vec.size - 1
			@transient_amplification_at_ky[ky] = -1
			eputs "No Min"
			next
		end
		first_min = t_element

		eputs "ky: #{ky}, first_min: #{first_min}"
	  loop do 
			t_element+=1
			#print t_element, ',', phi_vec.size
			new = phi_vec[t_element]
			break if new < old or t_element == phi_vec.size - 1
		end
		if t_element == phi_vec.size - 1
			@transient_amplification_at_ky[ky] = -1
			next
		end
		@transient_amplification_at_ky[ky] = phi_vec.subvector(t_element, phi_vec.size - t_element).max
	end
end

def max_trans_phi
	phivec = gsl_vector('phi2tot_over_time')
	offset = 30
	phivec.subvector(20, phivec.size - 20).max
end

def max_es_heat_amp(species_index)
	@transient_es_heat_flux_amplification_at_species_at_ky[species_index-1].values.max
end

def calculate_spectral_checks
	kx = gsl_vector('kx')
	ky = gsl_vector('ky')
	ky_spec = gsl_vector('spectrum_over_ky')
	kx_spec = gsl_vector('spectrum_over_kx')
	kpar_spec = gsl_vector('spectrum_over_kpar', ky_index: ky_spec.max_index + 1, kx_index: 1)
	
	@spectrum_check = []
	[kx_spec, ky_spec, kpar_spec].each do |spec|
		begin
			ends_max = [spec[0], spec[-1]].max + (10.0**(-9))
			p ends_max 		
			p spec.max
			check = (Math.log(spec.max/ends_max)/Math.log(10)).round
		rescue
			check= -10
		end
		@spectrum_check.push check
	end

  #Calculate peak kx, ky spectrum values and associated phi2 values
  @ky_spectrum_peak_idx = ky_spec.max_index 
  @ky_spectrum_peak_ky = ky[@ky_spectrum_peak_idx] 

  #Also want to know the phi2 at the energy containing scales and for ZFs
  #Pick phi2 at the final time step.
  phi_vec = gsl_vector('phi2_by_ky_over_time', ky_index:@ky_spectrum_peak_idx)
  @ky_spectrum_peak_phi2 = phi_vec[-1] 
  phi_vec = gsl_vector('phi2_by_ky_over_time', ky_index:1)
  @phi2_zonal = phi_vec[-1] 
end

def calculate_vspace_checks
	@vspace_check = ['lpc_pitch_angle', 'vres_pitch_angle', 'lpc_energy',  'vres_energy'].map do |name|
		saturated_time_average(name, {}) 
	end
		
end

alias :cvc :calculate_vspace_checks

def spec_chec(min, *dirns)
	return @spectrum_check.zip([0, 1, 2]).inject(true) do |bool, (check,dirn)|
		unless dirns.include? dirn
			bool and true
		else
			unless check >= min
				false
			else
				bool and true
			end
		end
	end
end

def sc(min)
	return @spectrum_check.min >= min
end
alias :csc :calculate_spectral_checks

end
end

