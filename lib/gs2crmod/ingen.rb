class CodeRunner
class Gs2 

	MAX_NAME_SIZE = 310


def warning(message)
	eputs "Warning: " + message; sleep 0.3
end

class InputFileError < StandardError
end

def error(message)
	raise InputFileError.new("Error: " + message)
end

def test_failed(namelist, var, gs2_var, tst)
	return  <<EOF

---------------------------
	Test Failed
---------------------------

Namelist: #{namelist}	
Variable: #{var}
GS2 Name: #{gs2_var}
Value: #{send(var)}
Test: #{tst[:test]}
Explanation: #{tst[:explanation]}

---------------------------
EOF

end


def namelist_test_failed(namelist, tst)
	return  <<EOF

---------------------------
	Test Failed
---------------------------

Namelist: #{namelist}	
Test: #{tst[:test]}
Explanation: #{tst[:explanation]}

---------------------------
EOF

end

# Checks input parameters for inconsistencies and prints a report. 

def run_namelist_tests(namelist, hash, enum = nil)
	ext = enum ? "_#{enum}" : ""
	hash[:must_pass].each do |tst|
		error(namelist_test_failed(namelist, tst)) unless instance_eval(tst[:test])
	end if hash[:must_pass]
	hash[:should_pass].each do |tst|
		warning(namelist_test_failed(namelist, tst)) unless instance_eval(tst[:test])
	end if hash[:should_pass]
	hash[:variables].each do |var, var_hash|
		gs2_var = (var_hash[:gs2_name] or var)
		cr_var = var+ext.to_sym 
		if send(cr_var) and (not var_hash[:should_include] or  eval(var_hash[:should_include]))
			var_hash[:must_pass].each do |tst|
				error(test_failed(namelist, cr_var, gs2_var, tst)) unless send(cr_var).instance_eval(tst[:test])
			end if var_hash[:must_pass]
			var_hash[:should_pass].each do |tst|
				warning(test_failed(namelist, cr_var, gs2_var, tst)) unless send(cr_var).instance_eval(tst[:test])
			end if var_hash[:should_pass]
			if (var_hash[:allowed_values] or var_hash[:text_options])
				tst = {test: "#{(var_hash[:allowed_values] or var_hash[:text_options]).inspect}.include? self", explanation: "The variable must have one of these values"}
				error(test_failed(namelist, cr_var, gs2_var, tst)) unless send(cr_var).instance_eval(tst[:test])
			end

		end
	end
end

	
# Eventually, this will be a full port of the tool of the same name in the GS2 folder. At the moment it runs a limited set of tests for common errors in the input parameters (including type checking).

def ingen
	
	# Sections
	
	# Namelist Tests
	# Grids
	# Parallelisation
	# Initialisation
  # Diagnostics	
	# Misc
	
	# Namelist Tests
	
	rcp.namelists.each do |namelist, hash|
		next if hash[:should_include].kind_of? String and not eval(hash[:should_include])
		if en = hash[:enumerator]
			#ep 'en', en, namelist
			next unless send(en[:name])
			send(en[:name]).times do |i|
				run_namelist_tests(namelist, hash, i+1)
			end
		else
			run_namelist_tests(namelist, hash)
		end
	end

	
	# Grid Errors
	
	# naky
	warning("Setting naky when non-linear mode is on is not recommended.") if @naky and @nonlinear_mode == "on"
	
	warning("You have set both ny and naky; naky will override ny.") if @ny and @naky

	error("Boundary options should not be periodic with finite magnetic shear") if @boundary_option == "periodic" and ((@s_hat_input and @s_hat_input.abs > 1.0e-6) or (@shat and @shat.abs > 1.0e-6))

	error("abs(shat) should not be less that 1.0e-6") if @shat and @shat.abs < 1.0e-6
	error("abs(s_hat_input) should not be less that 1.0e-6") if @s_hat_input and @s_hat_input.abs < 1.0e-6
	
	# delt 
	
	error("Please specify delt") unless @delt
	error("delt <= 0") if @delt <= 0.0
	warning("Nonlinear run with delt_minimum unspecified.") if @nonlinear_mode=="on" and not @delt_minimum

	error("delt (#@delt) < delt_minimum") if @delt and @delt_minimum and @delt < @delt_minimum

	# negrid
	#
	warning('negrid < 8 is not a good idea!') if @negrid and @negrid < 8

	# Parallelisation Errors
	
		
	# Check whether we are parallelising over x
	warning("Parallelising over x: suggest total number of processors should be: #{max_nprocs_no_x}") if actual_number_of_processors > max_nprocs_no_x and not @grid_option == "single"

	
	# Initialisation Errors
	
	# Check if restart folder exists
	if @restart_file and  @restart_file =~ /^(?<folder>[^\/]+)\//
		folder = $~[:folder]
		warning("Folder #{folder}, specified in restart_file, not present. NetCDF save may fail") unless FileTest.exist?(folder)
	end

	error("Setting @restart_file as an empty string will result in hidden restart files.") if @restart_file == ""

	error("ginit_option is 'many' but is_a_restart is false") if @ginit_option == "many" and not @is_a_restart

	#Diagnostic errors
	#
	#Check whether useful diagnostics have been omitted.

	not_set = [:write_verr, :save_for_restart, :write_nl_flux, :write_final_fields, :write_final_moments].find_all do  |diagnostic|
		not (send(diagnostic) and send(diagnostic).fortran_true?)
	end

	if not_set.size > 0
		str = not_set.inject("") do |str, diagnostic|
			str + "\n\t#{diagnostic} --- " + rcp.namelists[diagnostics_namelist][:variables][diagnostic][:description] rescue str
			end
		warning("The following useful diagnostics were not set:" + str) if str.length > 0
	end

	warning("You are running in nonlinear mode but have not switched the nonlinear flux diagnostic.") if not (@write_nl_flux and @write_nl_flux.fortran_true?) and @nonlinear_mode == "on" 

	#{
		#write_verr: "Velocity space diagnostics will not be output for this run"
	#}.each do |var, warn|
		#warning(v"#{var} not set or .false. --- " + warn) unless send(var) and send(var).fortran_true?
	#end
	
	warning("You will write out diagnostics less than 50 times") if @nstep/@nwrite < 50
	
	#Miscellaneous errors.
	
	error("The run name for this run is too long. Please move some of the variable settings to the local defaults file.") if @relative_directory.size + @run_name.size > MAX_NAME_SIZE

	warning("You are submitting a nonlinear run with no dissipation.") if @nonlinear_mode == "on" and @hyper_option=="none" and @collision_model=="none"

	warning("You have no spacial implicitness: (bakdif) for one of your species. Be prepared for numerical instabilities!") if (1..@nspec).to_a.find{|i| bd = send("bakdif_#{i}") and bd == 0}

	warning("The system will abort with rapid timestep changes...") if !@abort_rapid_time_step_change or @abort_rapid_time_step_change.fortran_true?



end

#  A hash which gives the actual numbers of gridpoints indexed by their corresponding letters in the layout string.

def gridpoints
 	gridpoints = {'l' => @ngauss, 'e' => @negrid, 's' => @nspec}
	if @grid_option == "single"
		gridpoints.absorb({'x'=>1, 'y'=>1})
	else
		gridpoints.absorb({'x' => (2.0 * (@nx - 1.0) / 3.0  + 1.0).floor,  'y' => (@naky or ((@ny - 1.0) / 3.0  + 1.0).floor)})
	end
	return gridpoints
end

def cumulative_gridpoints
	c = 1
	error("Please specify layout") unless @layout
  @layout.split(//).reverse.inject({}){|hash, let| c*=gridpoints[let]; hash[let] = c; hash}
end
# 	ep parallelisation
def	max_nprocs_no_x 
	parallelisation = cumulative_gridpoints
	parallelisation[parallelisation.keys[parallelisation.keys.index('x') - 1]]
end


	def diagnostics_namelist
		:gs2_diagnostics_knobs
	end
end
end

