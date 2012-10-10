##########################################
# = Code Runner GS2 Module              
##########################################
#
# Authors: Edmund Highcock
# Copyright: 2009 Edmund Highcock
#
# This is free software released under the GPL v3 
#
# This module allows easy running of the plasma turbulence simulation code gs2 using Code Runner, by automatically organising, naming and submitting runs, and analysing the run data.
#
# See Code Runner documentation, or documentation for individual methods.
#
# Notes
#
# index variables, e.g. kx_index, ky_index etc always refer to the 1-based Fortran index, to keep correspondance with the gs2 indices. Element variables, e.g. kx_element, always refer to the 0-based C/ruby index
#
# raw NumRu::NetCDF grids are in Fortran row-major order. This means that when you access grids using the NetCDF function NetCDF#get, you must specify the indices in fortran order (but 0-based!). The NetCDF#get function then returns a C-like NArray with the indices in the opposite order. You can convert this to a Ruby Array using the method NArray#to_a (the indices will still be in the same order).


#



class CodeRunner
	
	
	
	# This is a customised subclass of CodeRunner::Run which allows CodeRunner to submit and analyse simulations from the gyrokinetic flux tube code GS2, which is principally used for simulating plasmas in magnetic confinement fusion.
	#
	# It performs two distinct roles: submitting simulations and analysing the data.
	#
	# = Submitting Simulations
	#
	# This principally involves generating the input file, which is a very nontrivial task. In order to do this, it maintains a complete record of every possible input parameter for GS2, as well as what datatype that parameter is, and sometimes what values it is allowed to take. This allows that not only to generate the input file, but to check that the input file makes sense. However, although generating the input file works beautifully, the set of sanity checks that it makes is not exhaustive: intelligent use is still required!
	# 
	# In tandem with this, it maintains a whole set of tools for manipulating its database of input parameters. This includes updating their allowed values and also editing and accessing help for every input parameter.
	#
	# = Analysing Simulations
	#
	# The amount of analysis possible on GS2 data is enormous, and CodeRunner hasn't got close to getting there. What it can do is:
	# 
	# * Check if the run is complete by comparing the number of completed timesteps against nstep
	# * Calculate growth rates for linear runs.
	# * Check if non-linear runs have saturated and calculate fluxes for those runs.
	# * Automatically plot a huge variety of different graphs, ranging from simple plots of heat flux versus time to three-dimensional plots of the spectrum and potential.

class Gs2 < Run::FortranNamelist

GS2_CRMOD_VERSION = Version.new(Gem.loaded_specs['gs2crmod'].version.to_s)


def agk?
	false
end

MODULE_FOLDER = File.dirname(File.expand_path(__FILE__))

# Include the other files
@code_module_folder = folder = File.dirname(File.expand_path(__FILE__)) # i.e. the directory this file is in
setup_namelists(folder)
require folder + '/graphs.rb'
require folder + '/gsl_data.rb'
require folder + '/gsl_data_3d.rb'
require folder + '/check_convergence.rb'
require folder + '/calculations.rb'
require folder + '/ingen.rb'
require folder + '/properties.rb'
require folder + '/test_gs2.rb'

NaN = GSL::NAN
# GSL::Neg


def code_run_environment
	case CodeRunner::SYS
	when /helios/
		<<EOF
module load intel
module load bullxmpi
module load netcdf
EOF
	else
		""
	end
end

eval(%[
], GLOBAL_BINDING)


################################################
# Quantities that are calculated or determined by CodeRunner
# after the simulation has ended, i.e. quantities
# that are not available from the GS2 output files.
################################################


@results = [
	:converged,
	:decaying,
	:growth_rates,
	:real_frequencies,
	:growth_rates_by_ky, # deprecated
	:growth_rates_by_kx, # deprecated
	:growth_rate_at_ky,
	:growth_rate_at_kx,
	:growth_rate_at_ky_at_kx,
	:real_frequencies_by_ky,
	:max_growth_rate,
	:fastest_growing_mode,
	:freq_of_max_growth_rate,
	:ky,
	:gamma_r,
	:gamma_i,
	:run_time,
	:hflux_tot_stav,
	:phi2_tot_stav,
	:saturation_time_index,
	:es_heat_flux_stav,
	:es_mom_flux_stav,
	:hflux_tot_stav_error,
	:es_heat_flux_stav_error,
	:es_mom_flux_stav_error,
	:saturated,
	:shot_time,
	:spectrum_check,
	:par_mom_flux_stav,
	:perp_mom_flux_stav,
	:transient_amplification_at_kx,
	:transient_amplification_at_ky,
	:transient_amplification_at_ky_at_kx,
	:transient_es_heat_flux_amplification_at_species_at_kx,
	:transient_es_heat_flux_amplification_at_species_at_ky,
	:transient_es_heat_flux_amplification_at_species_at_ky_at_kx,
	:vspace_check
]


###############################################
# Other useful information about the run
###############################################

@gs2_run_info = [:time, :percent_of_total_time, :checked_converged, :is_a_restart, :restart_id, :restart_run_name, :completed_timesteps]

@run_info = @gs2_run_info.dup

##############################################################
# For backwards compatibility with CodeRunner version 0.5.0
##############################################################

@run_info_0_5_0 = {	
		time: :to_f,
		percent_of_total_time: :to_f,
		checked_converged: :to_b
	}

@results_0_5_0 = {
		converged: :to_b,
		decaying: :to_b,
		:growth_rates => :to_h,
		:real_frequencies => :to_h,
#		:ky_list => :to_h,
#		:kx_list => :to_h,
		:growth_rates_by_ky => :to_s,
		:real_frequencies_by_ky => :to_s,
		:max_growth_rate => :to_f,
		:fastest_growing_mode => :to_f,
		:freq_of_max_growth_rate => :to_f,
		:ky => :to_f,
		:gamma_r => :to_f,
		:gamma_i => :to_f,
		:run_time => :to_f
#		:theta_list => :to_h
	}

###############################################################

@uses_mpi = true

@modlet_required = false

@use_graphs = false
Phi = Struct.new("Phi", :phi, :ri, :theta_index, :kx_index, :ky_index)

@naming_pars = []

# def self.finish_setting_up_class
# 	@@variables += [
# end

# This method, as its name suggests, is called whenever CodeRunner is asked to analyse a run directory.this happens if the run status is not :Complete, or if the user has specified recalc_all(-A on the command line) or reprocess_all (-a on the command line).
#
# the structure of this function is very simple: first it calls get_status to determine the directory status, i.e. :Complete, :Incomplete, :NotStarted or :Failed, then it gets the time, which is the GS2 time at the end of the run, and it also gets the run_time, which is the wall clock time of the run. Finally,if non-linear mode is switched off, it calls calculate_growth_rates_and_frequencies, and if the non-linear mode is switched on, it calls calculate_time_averaged_fluxes.

def process_directory_code_specific
	run_namelist_backwards_compatibility

	unless @status == :Queueing
		get_status
	end
	
	eputs "Run #@status: #@run_name" if [:Complete,:Failed].include? @status 

	try_to_get_error_file
	@sys = @@successful_trial_system

	return if @status == :NotStarted or @status == :Failed or @status == :Queueing
	begin
		percent_complete = get_completed_timesteps/@nstep
		@percent_of_total_time = percent_complete
	rescue
		get_time
		@percent_of_total_time = @time / (@delt*@nstep) * 100.0	 rescue 0.0
	end
	return if @status == :Incomplete

	get_run_time

	calculate_results

end
def calculate_results
	return if ENV['CODE_RUNNER_NO_ANALYSIS'] =~ /true/

	
	eputs "Analysing run"	
	
	if @nonlinear_mode == "off"
	
		calculate_growth_rates_and_frequencies
	elsif @nonlinear_mode == "on"
		calculate_saturation_time_index
		calculate_time_averaged_fluxes
		calculate_spectral_checks
		begin 
			calculate_vspace_checks
		rescue
		end
	end

	@growth_rates ||={}
	@real_frequencies ||={}
end
	

# Try to read the runtime in minutes from the GS2 standard out.

def get_run_time
	logf(:get_run_time)
	output = @output_file || try_to_get_output_file
	return nil unless output
	begin
		Regexp.new("total from timer is:\\s*#{LongRegexen::NUMBER}", Regexp::IGNORECASE).match FileUtils.tail(output, 300) 
		logi $~
		@run_time = $~[:number].to_f
	rescue
		@run_time = nil
	end
end

# Output useful information from the NetCDF file. If no names are provided, output a list of all variables in the NetCDF file. <tt>names</tt> can either be a symbol or an array of symbols, in which case information will be output for the variables with those names. If values are provided, for example :dims,:get, :ndims, this information is retrieved from the file for every variable named.
#	ncdump
#	ncdump(:hflux)
#	ncdump([:hflux, :phi])
#	ncdump([:hflux, :phi], :dims)


def ncdump(names=nil, values=nil, extension = '.out.nc')
	names = [names] unless !names or names.class == Array
	names.map!{|name| name.to_s} if names
	pp NumRu::NetCDF.open(@run_name + extension).vars(names).to_a.sort{|var1, var2| var1.name <=> var2.name}.map{|var| values ? [var.name, var.send(values)] : var.name.to_sym}
end


#

def generate_phantom_runs
	@phantom_runs = []
	logf(:generate_phantom_runs)
	return if @grid_option == "single" and @scan_type == "none"
	begin 
		list(:ky) # This will fail unless the run has output the netcdf file
	rescue
		return
	end
	return unless @status == :Complete #and @converged
	log(@run_name)
	if @grid_option == "box" and @nonlinear_mode == "off"
		@ky = nil
# 		raise CRFatal.new("no @ky_list") unless @ky_list
# 		log list(:ky)
		list(:ky).each do |id, ky|
			phantom_run = create_phantom #self.dup
			phantom_run.ky = ky
			phantom_run.gamma_r = @growth_rates[ky]
			phantom_run.gamma_i = @real_frequencies[ky]
			log @runner.phantom_ids
# 			log('@runner.class', @runner.class)
# 			@runner.add_phantom_run(phantom_run)
		end
	elsif @scan_type and @scan_type != "none" 
		t = gsl_vector('t')
		scan_vals = gsl_vector('scan_parameter_value')
		current = scan_vals[0]
		start = 0
		for i in 0...t.size
			if scan_vals[i] != current
				phantom = create_phantom
				phantom.scan_index_window = [start+1, i] #remember indexes are elements + 1
				#ep 'scan_index_window', phantom.scan_index_window
				phantom.scan_parameter_value = current
				phantom.growth_rate_at_ky = nil
				phantom.growth_rate_at_kx = nil
				phantom.growth_rate_at_ky_at_kx = nil
				phantom.calculate_results
				current = scan_vals[i]
				start = i
			end
		end
	end
end



def get_time
	begin 
		lt = list(:t)
		return lt.values.max if lt.size>0
	rescue
	end
	time = nil
# 	eputs 	File.readlines(@run_name +".out").slice(-4..-1).reverse.join( "\n"); gets
	raise CRFatal.new("Couldn't find outfile #{@run_name}.out") unless FileTest.exist?(@run_name + ".out") 
	tail = FileUtils.tail("#@run_name.out", 4)
	#File.readlines(@run_name +".out").slice(-4..-1).reverse.join( "\n")
	tail.sub(LongRegexen::FLOAT) do
# 		eputs $~.inspect
		time =   $~[:float].to_f
	end  #if FileTest.exist? (@run_name +".out")
	#raise CRFatal.new("couldn't get the time from #{tail}") unless time
	@time = time
end

def get_completed_timesteps	
 	#raise CRFatal.new("Couldn't find outfile #{@run_name}.out") unless FileTest.exist?(@run_name + ".out")
	#p 'try to get completed_timesteps', Dir.pwd, 'nwrite', @nwrite, 'delt', @delt
	@completed_timesteps = (list(:t).size - 1) * (@nwrite || 1)
	#p 'tried to get completed_timesteps'
	#rescue
	#`grep time= #@run_name.out`.split.size
# 	File.read("#@run_name.out").scan(/^\s+time\s*=\s+/).size * @nwrite
end

def incomplete
	return (not 100 == percent_complete)
end

def parameter_transition(run)
end
# @@executable_location = nil
# def executable_location
# 	return "~/gs2_newterm" #(@@executable_location || ($gs2_new_term ? "~/gs2_newterm" : "~/gs2"))
# end
# 
# def executable_name
# 	"gs2"
# end

@code_long = "GS2 Gyrokinetic Flux Tube Code"

@excluded_sub_folders =[]

attr_accessor :theta_list, :ky_list, :ky_graphs, :eigenfunctions, :ky_list, :t_list
attr_accessor :scan_index_window, :scan_parameter_value

class << self
	aliold(:check_and_update)
	def check_and_update
		old_check_and_update
		@readout_list = (@variables + @results - [:growth_rates_by_ky, :growth_rates, :real_frequencies, :real_frequencies_by_ky, :ky_list, :kx_list, :theta_list, :t_list])
	end
end

def data_string
	logf(:data_string)
	return "" unless @converged unless @grid_option == 'single'
	logi(@ky, @growth_rates, @real_frequencies)
# 	log(:@@readout_list, @@readout_list)
	return rcp.readout_list.inject(""){|str,(var,type_co)| str+"#{(send(var) || "0")}\t"} + "\n" 

# 	@ky ? (@@variables + @@results - ).inject(""){|str,(var,type_co)| str+"#{(send(var) || "0")}\t"} + sprintf("%e\t%e\t%e\n", @ky, @growth_rates[@ky], @real_frequencies[@ky]) : (@@variables + @@results).inject(""){|str,(var,type_co)| str+"#{(send(var) || "0")}\t"} + sprintf("%e\t%e\t%e\n",  @fastest_growing_mode, @max_growth_rate, @freq_of_max_growth_rate)
end

def percent_complete
	@completed_timesteps ? @completed_timesteps.to_f / @nstep.to_f * 100.0 : @percent_of_total_time
end

def print_out_line
	logf(:print_out_line)
	name = @run_name
	name += " (res: #@restart_id)" if @restart_id
	name += " real_id: #@real_id" if @real_id
	beginning = sprintf("%2d:%d %-60s %1s:%2.1f(%s) %3s%1s %1s",  @id, @job_no, name, @status.to_s[0,1],  @run_time.to_f / 60.0, @nprocs.to_s, percent_complete, "%", @converged.to_s)
	 if @ky
		beginning += sprintf("%3s %4s %4s", @ky, @growth_rates[@ky], @real_frequencies[@ky])
	 elsif @nonlinear_mode == "off"
	    beginning += sprintf("%3s %4s %4s", 
	     @fastest_growing_mode, @max_growth_rate, 
	    @freq_of_max_growth_rate)
	 elsif @nonlinear_mode == "on"
# 		 p @hflux_tot_stav
		 beginning += "       sat:#{saturated.to_s[0]}" 
		 beginning += sprintf(" hflux:%1.2e", @hflux_tot_stav) if  @hflux_tot_stav 
		 beginning += sprintf("+/-%1.2e", @hflux_tot_stav_error) if  @hflux_tot_stav_error
		 beginning += sprintf(" momflux:%1.2e", @es_mom_flux_stav.values.sum) if @es_mom_flux_stav and @es_mom_flux_stav.values[0]
		 beginning += '  SC:' + @spectrum_check.map{|c| c.to_s}.join(',') if @spectrum_check 
		 beginning += '  VC:' + @vspace_check.map{|c| sprintf("%d", ((c*10.0).to_i rescue -1))}.join(',') if @vspace_check 
	 end
	 beginning += "  ---#{@comment}" if @comment
	 beginning
	
end    


def get_list_of(*args)
	#args can be any list of e.g. :ky, :kx, :theta, :t ... 
	logf(:get_list_of)
	refresh = args[-1] == true ? true : false
	args.pop if args[-1] == true
	logd
	Dir.chdir(@directory) do
		args.each do |var|
# 			eputs "Loading #{var}"
			list_name = var + :_list
			log list_name
			
# 			self.class.send(:attr_accessor, list_name)
 			next if (cache[list_name] and [:Failed, :Complete].include? status and not refresh)
			
			cache[list_name] = {}
			netcdf_file.var(var.to_s).get.to_a.each_with_index do |value, element|
# 				print '.'
				cache[list_name][element+1]=value
			end
#			eputs send(var+:_list)
		end
	end
	logfc :get_list_of
	return cache[args[0] + :_list] if args.size == 1
end

alias :list :get_list_of

def visually_check_growth_rate(ky=nil)
	logf :visually_check_growth_rate
	phi_vec = gsl_vector(:phi2_by_ky_over_time, {ky: ky})
	t_vec = gsl_vector(:t)
	constant, growth_rate = GSL::Fit::linear(t_vec, 0.5*GSL::Sf::log(phi_vec)).slice(0..1)
	eputs growth_rate

	graph = @@phi2tot_vs_time_template.graph(["#{constant} * exp (2 * #{growth_rate} * x)"], [[[t_vec, phi_vec], "u 1:2 title 'phi2tot #{@run_name}' w p"]], {"set_show_commands" => "\nset log y\n", "point_size"=>'1.0'})
# 	eputs graph.inline_data.inspect
	graph.show
	gets
	graph.kill

end


def show_graph
	thegraph = special_graph('phi2tot_vs_time_all_kys')
	thegraph.title += " for g_exb = #{@g_exb.to_f.to_s}"
	thegraph.show
	sleep 1.5
# 	@decaying = Feedback.get_boolean("Is the graph decaying?")
	thegraph.kill
end

# @@phi2tot_vs_time_template = {title: "Phi^2 Total vs Time", xlabel: " Time ", ylabel: "Phi^2 Total"})




def restart(new_run)
	#new_run = self.dup
	(rcp.variables).each{|v| new_run.set(v, send(v)) if send(v)}
	SUBMIT_OPTIONS.each{|v| new_run.set(v, self.send(v)) unless new_run.send(v)}
	#(rcp.results + rcp.gs2_run_info).each{|result| new_run.set(result, nil)}
	new_run.is_a_restart = true
	new_run.ginit_option = "many"
	new_run.delt_option = "default"
	#if Dir.entries(@directory).include? "nc"
		#old_restart_run_name =  (@restart_run_name or Dir.entries(@directory + '/nc').grep(/\.nc/)[0].sub(/\.nc\.\d+$/, ''))
		#new_run.restart_file = File.expand_path("#@directory/nc/#{old_restart_run_name}.nc")
	#else
		#new_run.restart_file = File.expand_path("#@directory/#@run_name.nc")
	#end
	new_run.restart_id = @id
	new_run.restart_run_name = @run_name
	@runner.nprocs = @nprocs if @runner.nprocs == "1" # 1 is the default so this means the user probably didn't specify nprocs 
	raise "Restart must be on the same number of processors as the previous run: new is #{new_run.nprocs.inspect} and old is #{@nprocs.inspect}" if !new_run.nprocs or new_run.nprocs != @nprocs
# 	@runner.parameters.each{|var, value| new_run.set(var,value)} if @runner.parameters
#   ep @runner.parameters
  new_run.run_name = nil
	new_run.naming_pars = @naming_pars
	new_run.update_submission_parameters(new_run.parameter_hash.inspect, false) if new_run.parameter_hash 
	new_run.naming_pars.delete(:restart_id)
	new_run.generate_run_name
	#@runner.submit(new_run)
	new_run
end	

# Return a list of restart file paths (relative to the run directory). 

def list_of_restart_files
	Dir.chdir(@directory) do
		files = Dir.entries.grep(/^\.\d+$/)
		files = Dir.entries.grep(/\.nc(?:\.\d|_ene)/) if files.size == 0
		if files.size == 0
			(Dir.entries.find_all{|dir| FileTest.directory? dir} - ['.', '..']).each do |dir|
				files = Dir.entries(dir).grep(/\.nc(?:\.\d|_ene)/).map{|file| dir + "/" + file}
				break if files.size == 0
			end
		end #if files.size == 0
		return files
	end # Dir.chdir(@directory) do
end

alias :lorf :list_of_restart_files

# Put restart files in the conventional location, i.e. nc/run_name.proc

def standardize_restart_files
	Dir.chdir(@directory) do
		FileUtils.makedirs('nc')
		list_of_restart_files.each do |file|
			proc_id = file.scan(/\.\d+$|_ene$/)[0]
			#p 'proc_id', proc_id
			FileUtils.mv(file, "nc/#@run_name.nc#{proc_id}")
		end
	end
end

# Delete all the restart files (irreversible!)
#

def delete_restart_files(options={})
	return unless Feedback.get_boolean("Deleting restart files. This action cannot be reversed. Do you wish to continue?") unless options[:no_confirm]
	list_of_restart_files.each{|file| FileUtils.rm file}
end


	


def species_letter
	species_type(1).downcase[0,1]
end

def species_type(index)	
	if rcp.variables.include? :type_1 
		type = send(:type_ + index.to_sym)
	else
		types = rcp.variables.find_all{|var| var.to_s =~ /^type/}.map{|var| send(var)}
		type = types[index.to_i - 1]
	end
	type
end


# Returns true if this run  has not been restarted, false if it has. This allows one to get data from the final run of  a series of restarts.

def no_restarts
	raise NoRunnerError unless @runner
	!(@runner.runs.find{|run| run.restart_id == @id})
end


def restart_chain
	if @restart_id
		return @runner.run_list[@restart_id].restart_chain
	end
	chain = []
	currid = @id
	loop do
		chain.push currid
		break unless (restrt = @runner.runs.find{|run| run.restart_id == currid})
		currid = restrt.id
	end
	return chain
end




	

def get_status
# 	eputs 'Checking Status'
	logf(:get_status)

	Dir.chdir(@directory) do
		if @running
			if FileTest.exist?(@run_name + ".out") and FileUtils.tail(@run_name + ".out", 5).split(/\n/).size > 4 and FileUtils.tail(@run_name + ".out", 200) =~ /t\=/
				@status = :Incomplete
			else
				@status = :NotStarted
			end
			
		else	
			if FileTest.exist?(@run_name + ".out") and FileUtils.tail(@run_name + ".out", 5).split(/\n/).size > 4
				#eputs "HERE", @scan_type
				if  @nonlinear_mode == "off" and FileUtils.tail(@run_name + ".out",200) =~ /omega converged/
					eputs 'Omega converged...'
					@status = :Complete
			  elsif @scan_type and @scan_type != "none" and FileUtils.tail(@run_name + ".par_scan",200) =~ /scan\s+is\s+complete/i
					eputs 'Scan complete...'
					@status = :Complete
				elsif @nonlinear_mode == "on" or !@omegatol or @omegatol < 0.0 or (@exit_when_converged and @exit_when_converged.fortran_false?)
				   	eputs 'No omegatol'
					if FileTest.exist?(@run_name + ".out.nc")
						get_completed_timesteps
					else		
						eputs "Warning: no netcdf file #@run_name.out.nc"
						@status = :Failed
						return
					end
						#ep "completed_timesteps", @completed_timesteps
					eputs "#{percent_complete}% of Timesteps Complete"
					if percent_complete == 100.0
						@status = :Complete
					elsif percent_complete > 5 and FileUtils.tail(output_file, 200) =~ /total from timer is/
						@status = :Complete
					else 
						@status = :Failed
					end		
				else
					@status = :Failed
				end
			else 
				@status=:Failed
			end
		end
	end
end


	def self.modify_job_script(runner, runs, script)
		if CODE_OPTIONS[:gs2] and CODE_OPTIONS[:gs2][:list]
			if (list_size = CODE_OPTIONS[:gs2][:list]).kind_of? Integer
				raise "The total number of runs must be a multiple of the list size!" unless runs.size % list_size == 0
				pieces = runs.pieces(runs.size/list_size)
			else
				pieces = [runs]
			end
			script = ""
			pieces.each do |runs|
				#ep 'there is a list'
				FileUtils.makedirs('job_lists')
				jid = "#{runs[0].id}-#{runs[-1].id}"
				list_file = "job_lists/gs2_list_#{jid}.list"
				File.open(list_file,'w') do |file|
					file.puts runs.size
					file.puts runs.map{|r| "#{r.relative_directory}/#{r.run_name}"}.join("\n")
				end
				raise "runs must all have the same nprocs" unless runs.map{|r| r.nprocs}.uniq.size == 1 
				runs.each do |r| 
					# Make sure the restart file name includes the relative directory for
					# list runs
					reldir = r.relative_directory 
					rdir = r.restart_dir
					#puts rdir[0...reldir.size] == reldir, rdir[0...reldir.size], reldir
					#raise ""
					if rdir
						r.restart_dir = reldir + '/' + rdir if not rdir[0...reldir.size] == reldir
					else
						r.restart_dir = reldir
					end
					Dir.chdir(r.directory){r.write_input_file}
				end
				np = runs[0].nprocs.split('x').map{|n| n.to_i}
				np[0] *= runs.size
				nprocs = np.map{|n| n.to_s}.join('x')
				@runner.nprocs = nprocs
				ls = ListSubmitter.new(@runner, nprocs, list_file, jid)
				script << ls.run_command 
			end
		end
		return script
	end	

	class ListSubmitter
		include CodeRunner::SYSTEM_MODULE
		@uses_mpi = true
		attr_reader :executable_location, :executable_name, :parameter_string
		attr_reader :job_identifier
		def initialize(runner, nprocs, list_file, jid)
			@executable_location = runner.executable_location
			@executable_name = runner.executable_name
			@parameter_string = list_file
			@job_identifier = jid
			@nprocs = nprocs
		end
		def rcp 
			self.class.rcp
		end
		def self.rcp
			@rcp ||= CodeRunner::Run::RunClassPropertyFetcher.new(self)
		end

	end #class ListSubmitter

def recheck
	logf(:recheck)
	Dir.chdir(@directory) do
		logi('@runner.object_id', @runner.object_id)
		log('@runner.class',  @runner.class)
		runner = @runner
		instance_variables.each{|var| instance_variable_set(var, nil) unless var == :@runner}
		begin File.delete(".code_runner_run_data") rescue Errno::ENOENT end
		begin File.delete("code_runner_results.rb") rescue Errno::ENOENT end
		logi(:@checked_converged, @checked_converged)
		logi('@runner.object_id after reset', @runner.object_id)		
		log('@runner.class',  @runner.class)
		process_directory
	end
end


def generate_input_file
	raise CRFatal("No Input Module File Given or Module Corrupted") unless methods.include? (:input_file_text)
	run_namelist_backwards_compatibility
	if @restart_id
		@runner.run_list[@restart_id].restart(self)
		eputs 'Copying Restart files', ''
		FileUtils.makedirs('nc')
		#old_dir = File.dirname(@restart_file)
		@restart_file = "#@restart_run_name.nc" #+ File.basename(@restart_file) #.sub(/\.nc/, '')
		@restart_dir = "nc"
		#files = Dir.entries(old_dir).grep(/\.nc(?:\.\d|_ene)/)
		#files = Dir.entries(old_dir).grep(/^\.\d+$/) if files.size == 0
		files = @runner.run_list[@restart_id].list_of_restart_files.map do |file|
			@runner.run_list[@restart_id].directory + "/" + file
		end
		files.each_with_index do |file , index|
			eputs "\033[2A" # Terminal jargon - go back one line
			eputs "#{index+1} out of #{files.size}"
			num = file.scan(/(?:\.\d+|_ene)$/)[0]
			#FileUtils.cp("#{old_dir}/#{file}", "nc/#@restart_file#{num}")
			FileUtils.cp(file, "nc/#@restart_file#{num}")
		end
	elsif @save_for_restart.fortran_true?
		@restart_dir = "nc"
		#if CODE_OPTIONS[:gs2] and CODE_OPTIONS[:gs2][:list]
			#FileUtils.makedirs "#{@runner.root_folder}/#@restart_dir"
		#else
			FileUtils.makedirs @restart_dir
		#end
		@restart_file = "#@run_name.nc"

	end
	
	# Let Gs2 know how much wall clock time is available. avail_cpu_time is a GS2 input parameter.
	@avail_cpu_time = @wall_mins * 60 if @wall_mins

	#  Automatically set the number of  nodes to be the maximum possible without parallelising over x, if the user has left the number of nodes unspecified.
	
	set_nprocs


	######### Check for errors and inconsistencies 
	ingen
	#########

	write_input_file
end

def write_input_file
	File.open(@run_name + ".in", 'w'){|file| file.puts input_file_text}
end

def set_nprocs

	if (nprocs_in = @nprocs) =~ /^x/
		max = max_nprocs_no_x
		nodes = 0
		@nprocs = "#{nodes}#{nprocs_in}"	 
		loop do
			nodes += 1
			@nprocs = "#{nodes}#{nprocs_in}"	 
			if actual_number_of_processors > max 
				nodes -= 1
				@nprocs = "#{nodes}#{nprocs_in}"	 
				break
			end
		end
	end
end

def actual_number_of_processors
  raise "Please specify the processor layout using the -n or (n:) option" unless @nprocs
	@nprocs.split('x').map{|n| n.to_i}.inject(1){|ntot, n| ntot*n}
end


def parameter_string
		return "#{@run_name}.in"
end


def self.list_code_commands
	puts (methods - superclass.methods).sort
end

def self.add_variable_to_namelist(namelist, var, value)
	var = :stir_ + var if namelist == :stir
	super(namelist, var, value)
end

def input_file_header
		run_namelist_backwards_compatibility
	<<EOF
!==============================================================================
!  		GS2 INPUT FILE automatically generated by CodeRunner 
!==============================================================================
!
!  GS2 is a gyrokinetic flux tube initial value turbulence code 
!  which can be used for fusion or astrophysical plasmas.
!  
!  	See http://gyrokinetics.sourceforge.net
!
!  CodeRunner is a framework for the automated running and analysis 
!  of large simulations. 
!
!  	See http://coderunner.sourceforge.net
!      by CodeRunner version #{CodeRunner::CODE_RUNNER_VERSION.to_s}
!
!==============================================================================

EOF
end

def self.defaults_file_header
	<<EOF1
######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created #{Time.now.to_s}                                           #
#                                                                    #
######################################################################

@defaults_file_description = ""
EOF1
end


# Customize this method from Run::FortranNamelist by saying when diagnostics are not switched on.

def namelist_text(namelist, enum = nil)
	hash = rcp.namelists[namelist]
	text = ""
	ext = enum ? "_#{enum}" : ""
	text << "!#{'='*30}\n!#{hash[:description]} #{enum} \n!#{'='*30}\n" if hash[:description]
	text << "&#{namelist}#{ext}\n"
	hash[:variables].each do |var, var_hash|
		code_var = (var_hash[:code_name] or var)
		cr_var = var+ext.to_sym 
# 		ep cr_var, namelist
		if send(cr_var) and (not var_hash[:should_include] or  eval(var_hash[:should_include]))
# 				var_hash[:tests].each{|tst| eval(tst).test(send(cr_var), cr_var)}
			if String::FORTRAN_BOOLS.include? send(cr_var) # var is a Fortran Bool, not really a string
				output = send(cr_var).to_s
			elsif (v = send(cr_var)).kind_of? Complex
				output = "(#{v.real}, #{v.imag})"
			else
				output = send(cr_var).inspect
			end
			text << " #{code_var} = #{output} #{var_hash[:description] ? "! #{var_hash[:description]}": ""}\n"
		elsif namelist == :gs2_diagnostics_knobs or namelist == :diagnostics
			text << "  ! #{code_var} not specified --- #{var_hash[:description]}\n"
		end
	end
# # 	end
	text << "/\n\n"
	text
end

# def self.add_code_var
# 	rcp.namelists.each do |namelist, hash|
# 		hash[:variables].each do |var, var_hash|
# 			p var
# 			var_hash[:code_name] = var_hash[:gs2_name] if var_hash[:gs2_name]
# 		end
# 	end
# 	save_namelists
# end
		
	
def update_physics_parameters_from_miller_input_file(file)
	hash = self.class.parse_input_file(file)
	hash[:parameters].each do |var, val|
		set(var,val)
	end
	hash[:theta_grid_parameters].each do |var, val|
		next if  [:ntheta, :nperiod].include? var
		set(var, val)
	end
	hash[:dist_fn_knobs].each do |var, val|
		next unless [:g_exb].include? var
		set(var, val)
	end
	hash[:theta_grid_eik_knobs].each do |var, val|
		next unless [:s_hat_input, :beta_prime_input].include? var
		set(var, val)
	end
	
	hash[:species_parameters_2].each do |var, val|
		#next unless [:s_hat_input, :beta_prime_input].include? var
		set((var.to_s + '_2').to_sym, val)
	end
	hash[:species_parameters_1].each do |var, val|
		#next unless [:s_hat_input, :beta_prime_input].include? var
		set((var.to_s + '_1').to_sym, val)
	end
end



def renew_info_file
	Dir.chdir(@directory){make_info_file("#@run_name.in")}
end
	
# This method overrides a method defined in heuristic_run_methods.rb in the CodeRunner source. It is called when CodeRunner cannot find any of its own files in the folder being analysed. It takes a GS2 input file and generates a CodeRunner info file. This means that GS2 runs which were not run using CodeRunner can nonetheless be analysed by it. In order for it to be called the -H flag must be specified on the command line.
	
def run_heuristic_analysis
	ep 'run_heuristic_analysis', Dir.pwd
	infiles = Dir.entries.grep(/^[^\.].*\.in$/)
	ep infiles
	raise CRMild.new('No input file') unless infiles[0]
	raise CRError.new("More than one input file in this directory: \n\t#{infiles}") if infiles.size > 1
	input_file = infiles[0]
	ep 'asdf'
	@nprocs ||= "1"
	@executable ||= "/dev/null"
	make_info_file(input_file, false)
end
	
@source_code_subfolders = ['utils', 'geo']

attr_accessor :iphi00, :saturation_time #Necessary for back. comp. due to an old bug

folder = File.dirname(File.expand_path(__FILE__)) # i.e. the directory this file is in

SPECIES_DEPENDENT_NAMELISTS = eval(File.read(folder + '/species_dependent_namelists.rb'), binding, folder + '/species_dependent_namelists.rb')
# 
SPECIES_DEPENDENT_VARIABLES_WITH_HELP = SPECIES_DEPENDENT_NAMELISTS.values.inject({}) do |hash, namelist_hash|
	namelist_hash[:variables].each do |var, var_hash|
			hash[var] = var_hash[:help]
	end
	hash
end

SPECIES_DEPENDENT_VARIABLES = SPECIES_DEPENDENT_VARIABLES_WITH_HELP.keys
SPECIES_DEPENDENT_VARIABLES.each{|var| attr_accessor var} # for backwards compatibility

['i', 'e'].each do |n| 
	SPECIES_DEPENDENT_VARIABLES_WITH_HELP.each do |name, help|
		attr_accessor name + "_#{n}".to_sym #for backwards compatibility 
 	end
end

old_vars = %w[
	:TiTe
	:Rmaj
	:R_geo
	:invLp_input
	:D_hypervisc
	:D_hyperres
	:D_hyper
	:C_par
	:C_perp
].map{|n| n.to_s.sub(/^:/, '').to_sym}

old_vars.each do |var|
	alias_method(var, var.to_s.downcase.to_sym)
	alias_method("#{var}=".to_sym, "#{var.downcase}=".to_sym)
end




def run_namelist_backwards_compatibility
	SPECIES_DEPENDENT_VARIABLES.each do |var|
		set(var + "_1".to_sym, (send(var + "_1".to_sym) or send(var + "_i".to_sym) or send(var)))
		set(var + "_2".to_sym, (send(var + "_2".to_sym) or send(var + "_e".to_sym)))
	end
end


def stop
	`touch #@directory/#@run_name.stop`
end

  def vim_output
		system "vim -Ro #{output_file} #{error_file} #@directory/#@run_name.error #@directory/#@run_name.out "
	end
	alias :vo :vim_output
  def vim_stdout
		system "vim -Ro #{output_file} "
	end
	alias :vo1 :vim_stdout
	def plot_efit_file
		Dir.chdir(@directory) do
			text = File.read(@eqfile)
			text_lines = text.split("\n")
			first_line = text_lines[0].split(/\s+/)
			second_line = text_lines[1].split(/\s+/)
			nr = first_line[-2].to_i
			nz = first_line[-1].to_i
			rwidth = second_line[1].to_f
			zwidth = second_line[2].to_f
			rmag = second_line[3].to_f
			nlines = (nr.to_f/5.0).ceil
			nlines_psi = ((nr*nz).to_f/5.0).ceil
			start = 5
			f = text_lines[start...(start+=nlines)].join(" ").split(nil).map{|s| s.to_f}.to_gslv
			pres = text_lines[(start)...(start += nlines)].join(" ").split(nil).map{|s| s.to_f}.to_gslv 
			dumy = text_lines[(start)...(start += nlines)].join(" ").split(nil).map{|s| s.to_f}.to_gslv
			ffprime = text_lines[(start)...(start+= nlines)].join(" ").split(nil).map{|s| s.to_f}.to_gslv
			psi = text_lines[(start)...(start += nlines_psi)].join(" ")
			q = text_lines[(start)...(start += nlines)].join(" ").split(nil).map{|s| s.to_f}.to_gslv
			nbound = text_lines[start...start+=1].join(" ").to_i
			rz = text_lines[(start)...(start += nbound*2)].join(" ").split(/\s+/)
			rz.shift
			rbound, zbound, dummy = rz.inject([[], [], true]){|arr,val| arr[2] ? [arr[0].push(val), arr[1], false] : [arr[0], arr[1].push(val), true]}
			#rbound.shift

			psi = psi.split(/\s+/)
			psi.shift
			psi.map!{|v| v.to_f}
			psi_arr = SparseTensor.new(2)
			k = 0
			for i in 0...nz
				for j in 0...nr
					psi_arr[j,i] = psi[k]
					k+=1
				end
			end
			kit = GraphKit.quick_create([((0...nr).to_a.to_gslv - nr/2 - 1 )/(nr-1)*rwidth+rmag, ((0...nz).to_a.to_gslv-nz/2 + 1)/(nz-1) * zwidth, psi_arr], [rbound, zbound, rbound.map{|r| 0}])
			kit.gp.contour = ""
			kit.gp.view = "map"
			#kit.gp.nosurface = ""
			kit.gp.cntrparam = "levels 20"
			kit.data[0].gp.with = 'l'
			kit.data[1].gp.with = 'l lw 2 nocontours'
			kit.gnuplot
			
			kit2 = GraphKit.quick_create([pres/pres.max],[f/f.max],[q/q.max])
			kit2.data[0].title = 'Pressure/Max Pressure'
			kit2.data[1].title = 'Poloidal current function/Max poloidal current function'
			kit2.data[2].title = 'Safety factor/Max Safety factor'
			kit2.gnuplot
	    		


			#p ['f', f, 'p', pres, 'ffprime', ffprime, 'nlines', nlines, 'psi', psi, 'q', q, 'nbound', nbound, 'rbound', rbound, 'zbound', zbound]


		end
	end

	end # class GS2
	# For backwards compatibility

Gs2BoxNtRun = Gs2CycloneRun = Gs2BoxCollisionalRun = Gs2Jet42982Run = Gs2ArtunRun = Gs2LinskerRun = Gs2BarnesLinskerRun = Gs2BoxMovieRun = Gs2Run = Gs2
end # class CodeRunner

# ep CodeRunner::Gs2CycloneRun.ancestors

begin     
	require "numru/netcdf"
rescue LoadError
    eputs "Error: No NetCDF: data analysis for gs2 not possible"
end

class Float
      def <=>(other) # necessary because of netcdf quirks
				
      	  d = (self - other)
	  if d.abs / (self.abs + 1) < 1e-10
	     return 0
	  else 
	     return (d / d.abs).to_i
	  end
      end
      def ==(other)
				return false unless other.kind_of? Numeric
      	  return (self - other).abs < 1e-14
      end
end 
    
class Hash
	
#	puts self

	def convert_to_index(run, *names)
			if self[:strongest_non_zonal_mode]
	         ky_element, kx_element =  run.gsl_matrix('spectrum_over_ky_over_kx', no_zonal: true).max_index
	         p self[:kx_index] = kx_element + 1
	         p self[:ky_index] = ky_element + 1
					 self[:strongest_non_zonal_mode] = false
	    end


# 		ep run
		names.each do |name|
			if name == :kx
				if lkx = self[:lagrangian_kx]
					self[:lagrangian_kx_index] = list(:kx).key(lkx)
				end
				if lkxi = self[:lagrangian_kx_index] ||= self[:lkx_index]
					self[:kx_index] = run.eulerian_kx_index(kx_index: lkxi, ky_index: self[:ky_index], t_index: self[:t_index])
				end
			end

			#ep 'name', name
			self[:ky_index] = 1 if name == :ky and run.grid_option == "single" 
			self[:kx_index] = 1 if name == :kx and run.grid_option == "single" 
# 			ep run.list(name)
			self[name + :_index] ||= run.list(name).key(self[name]) || (raise ("#{name} not specified"))
		end
		
	end
	def setup_time_window
		self[:t_index_window] ||= [self[:t_index],self[:t_index]] if self[:t_index]
		self[:begin_element], self[:end_element] = (self[:t_index_window] ? self[:t_index_window].map{|ind| ind -1} : [0, -1])
	end

end
