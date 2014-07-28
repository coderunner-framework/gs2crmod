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
    #gs2_var = (var_hash[:gs2_name] or var)
    cr_var = var+ext.to_sym
    value = send(cr_var)
    if value.kind_of? Array
      value.each{|v| test_variable(namelist, var, var_hash, ext, v)}
    else
      test_variable(namelist, var, var_hash, ext, value)
    end
  end
end

def test_variable(namelist, var, var_hash, ext, value)
    gs2_var = (var_hash[:gs2_name] or var)
    cr_var = var+ext.to_sym
    if value and (not var_hash[:should_include] or  eval(var_hash[:should_include]))
      var_hash[:must_pass].each do |tst|
        error(test_failed(namelist, cr_var, gs2_var, tst)) unless value.instance_eval(tst[:test])
      end if var_hash[:must_pass]
      var_hash[:should_pass].each do |tst|
        warning(test_failed(namelist, cr_var, gs2_var, tst)) unless value.instance_eval(tst[:test])
      end if var_hash[:should_pass]
      if (var_hash[:allowed_values] or var_hash[:text_options])
        tst = {test: "#{(var_hash[:allowed_values] or var_hash[:text_options]).inspect}.include? self", explanation: "The variable must have one of these values"}
        error(test_failed(namelist, cr_var, gs2_var, tst)) unless value.instance_eval(tst[:test])
      end

    end
end


# Eventually, this will be a full port of the ingen tool in the GS2 folder. At the moment it runs a limited set of tests for common errors in the input parameters (including type checking).

def check_parameters

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

  ###############
  # Grid Errors #
  ###############

  # naky
  warning("Setting naky when non-linear mode is on is not recommended.") if @naky and @nonlinear_mode == "on"

  warning("You have set both ny and naky; naky will override ny.") if @ny and @naky

  error("abs(shat) should not be less that 1.0e-6") if @shat and @shat.abs < 1.0e-6 and not agk?
  error("abs(s_hat_input) should not be less that 1.0e-6") if @s_hat_input and @s_hat_input.abs < 1.0e-6 and not agk?

  # delt

  error("Please specify delt") unless @delt
  error("delt <= 0") if @delt <= 0.0
  warning("Nonlinear run with delt_minimum unspecified.") if @nonlinear_mode=="on" and not @delt_minimum

  error("delt (#@delt) < delt_minimum") if @delt and @delt_minimum and @delt < @delt_minimum

  # negrid
  warning('negrid < 8 is not a good idea!') if @negrid and @negrid < 8

    # nakx
  warning("You have set both nx and ntheta0; ntheta0 will override nx.") if @nx and @ntheta0

  warning("Do you have a reason for setting equal_arc = true (default)? If not set false.") if @equilibrium_option=="eik" and (!@equal_arc or @equal_arc.fortran_true?)

  warning("Recommend nperiod > 1 for linear runs.") if @nonlinear_mode == "off" and (!@nperiod or @nperiod == 1)
  warning("Recommend nperiod = 1 for nonlinear runs.") if @nonlinear_mode == "on" and (@nperiod > 1)

  warning("Consider using field_option = local and associated optimizations.") if @field_option and @field_option == "implicit"

  #################################
  # Parallelisation/Layout Errors #
  #################################

  # Best linear run layout is lexys
  warning("The best layout for linear runs is usually lexys.") if @nonlinear_mode=="off" and not @layout=="lexys"

  # Best nonlinear run layout is xyles
        warning("The best layout for nonlinear runs is usually xyles.") if @nonlinear_mode=="on" and not @layout=="xyles"

  # Check whether we are parallelising over x
  warning("Parallelising over x: suggest total number of processors should be: #{max_nprocs_no_x}") if actual_number_of_processors > max_nprocs_no_x and not @grid_option == "single"

  #########################
  # Initialisation Errors #
  #########################

  # Check if restart folder exists
  if @restart_file and  @restart_file =~ /^(?<folder>[^\/]+)\//
    folder = $~[:folder]
    warning("Folder #{folder}, specified in restart_file, not present. NetCDF save may fail") unless FileTest.exist?(folder)
  end

  error("Setting @restart_file as an empty string will result in hidden restart files.") if @restart_file == ""

  error("ginit_option is 'many' but is_a_restart is false") if @ginit_option == "many" and not @is_a_restart

  error("chop_side should not be used (remove test if default changes from T to F)") if !@chop_side or @chop_side.fortran_true?

  #####################
  # Diagnostic errors #
  #####################

  #Check whether useful diagnostics have been omitted.

  not_set = [:write_verr, :save_for_restart, :write_nl_flux, :write_final_fields, :write_final_moments].find_all do  |diagnostic|
    not (send(diagnostic) and send(diagnostic).fortran_true?)
  end

  if not_set.size > 0
    str = not_set.inject("") do |s, diagnostic|
      s + "\n\t#{diagnostic} --- " + rcp.namelists[diagnostics_namelist][:variables][diagnostic][:description] rescue s
    end
    warning("The following useful diagnostics were not set:" + str) if str.length > 0
  end

  warning("You are running in nonlinear mode but have not switched the nonlinear flux diagnostic.") if not (@write_nl_flux and @write_nl_flux.fortran_true?) and @nonlinear_mode == "on"

  #{
    #write_verr: "Velocity space diagnostics will not be output for this run"
  #}.each do |var, warn|
    #warning(v"#{var} not set or .false. --- " + warn) unless send(var) and send(var).fortran_true?
  #end

  error("Please specify nwrite") unless @nwrite
  error("Please specify nstep") unless @nstep


  warning("You will write out diagnostics less than 50 times") if @nstep/@nwrite < 50

  ########################
  # Miscellaneous errors #
  ########################

  error("The run name for this run is too long. Please move some of the variable settings to the local defaults file.") if @relative_directory.size + @run_name.size > MAX_NAME_SIZE

  warning("You are submitting a nonlinear run with no dissipation.") if @nonlinear_mode == "on" and @hyper_option=="none" and @collision_model=="none"

  warning("You have no spacial implicitness: (bakdif) for one of your species. Be prepared for numerical instabilities!") if (1..@nspec).to_a.find{|i| bd = send("bakdif_#{i}") and bd == 0}

  warning("The system will abort with rapid timestep changes...") if !@abort_rapid_time_step_change or @abort_rapid_time_step_change.fortran_true?

  warning("local_field_solve is an old variable that should not really be used.") if @local_field_solve and  @local_field_solve.fortran_true?

  #############################
  # Boundary Condition Errors #
  #############################

  error("Boundary options should be linked with finite magnetic shear.") if (!@boundary_option or @boundary_option != "linked") and ((@s_hat_input and @s_hat_input.abs > 1.0e-6) or (@shat and @shat.abs > 1.0e-6))

  error("Set nonad_zero = true.") if @nonad_zero and not @nonad_zero.fortran_true?


  ###################
  # Spectrogk tests #
  ###################
  #
  if spectrogk?
    if @force_5d and @force_5d.fortran_true?
      warning("Must specify interpolation method with phi_method.") if not (@phi_method)
    end
  end

  ################
  # Damping Rate #
  ################

  error("Linear runs with hyperviscosity are NOT recommended!") if @nonlinear_mode=="off" and (@hyper_option and @hyper_option=="visc_only") and (@d_hypervisc and @d_hypervisc!=0)

  warning("Amplitude dependent part of hyperviscosity being ignored since const_amp = true") if (@hyper_option and @hyper_option=="visc_only") and (@const_amp and @const_amp.fortran_true?)

  ###################
  # Geometry Errors #
  ###################

  error("You must set bishop = 4 for Miller(local) geometry. Remember also that s_hat_input will override shat") if (@bishop!=4 and (@local_eq and @local_eq.fortran_true?))

  error("Shift should be > 0 for s-alpha equilibrium.") if @equilibrium_option=="s-alpha" and (@shift and @shift < 0)
  error("Shift should be < 0 for Miller equilibrium.") if @equilibrium_option=="eik" and @local_eq.fortran_true? and (@shift and @shift > 0)

  error("irho must be 2 for Miller equilibrium.") if @equilibrium_option=="eik" and @local_eq.fortran_true? and (@irho and @irho!=2)

  warning("Note that shat != s_hat_input") if @shat and @s_hat_input and @shat!=@s_hat_input

  ##################
  # Species Errors #
  ##################

  error("Must set z = -1 for electron species.") if (@type_2 and @z_2 and @type_2=='electron' and @z_2 != -1)


  #################
  # Optimisations #
  #################

  if CODE_OPTIONS[:gs2] and CODE_OPTIONS[:gs2][:show_opt]
    eputs("Optimisation Summary:")
    optimisation_flags.each do |flag|
      eputs("-------------------------  #{flag}: #{send(flag)}\n* #{rcp.variables_with_help[flag].gsub(/\n/, "\n\t").sub(/\A([^.]*.).*\Z/m, '\1')}") 
    end
    #not_set = [:operator, :save_for_restart, :write_nl_flux, :write_final_fields, :write_final_moments].find_all do  |diagnostic|
      #not (send(diagnostic) and send(diagnostic).fortran_true?)
    #end

    #if not_set.size > 0
      #str = not_set.inject("") do |s, diagnostic|
        #s + "\n\t#{diagnostic} --- " + rcp.namelists[diagnostics_namelist][:variables][diagnostic][:description] rescue s
      #end
      #warning("The following useful diagnostics were not set:" + str) if str.length > 0
    #end
  end
  
 


end

def optimisation_flags
  [
    :opt_redist_persist,
    :opt_redist_persist_overlap,
    :opt_redist_nbk,
    :opt_redist_init,
    :intmom_sub,
    :intspec_sub,
    #:local_field_solve,
    :do_smart_update,
    :field_subgath,
    :field_option,
    :field_local_allreduce,
    :field_local_allreduce_sub,
    :minnrow,
    :opt_init_bc,
    :opt_source
  ]
end

#  A hash which gives the actual numbers of gridpoints indexed by their corresponding letters in the layout string.

def gridpoints
  gridpoints = {'l' => @ngauss, 'e' => @negrid, 's' => @nspec}
  if @grid_option == "single"
    gridpoints.absorb({'x'=>1, 'y'=>1})
  else
    gridpoints.absorb({'x' => (@ntheta0 or (2.0 * (@nx - 1.0) / 3.0  + 1.0).floor),  'y' => (@naky or ((@ny - 1.0) / 3.0  + 1.0).floor)})
  end
  return gridpoints
end

def cumulative_gridpoints
  c = 1
  error("Please specify layout") unless @layout
  @layout.split(//).reverse.inject({}){|hash, let| c*=gridpoints[let]; hash[let] = c; hash}
end
#   ep parallelisation
def max_nprocs_no_x
  parallelisation = cumulative_gridpoints
  parallelisation[parallelisation.keys[parallelisation.keys.index('x') - 1]]
end


  def diagnostics_namelist
    :gs2_diagnostics_knobs
  end

  # Run the ingen tool on the input file
  def ingen
    Dir.chdir(@directory) do
      ing = File.dirname(File.expand_path(@executable)) + '/ingen'
      success = system "#{ing} #@run_name.in"
      warning("Could not run ingen... make sure that ingen is in the same folder as @executable and can be run on the login nodes if you want this to work") unless success
    end
  end
end
end


