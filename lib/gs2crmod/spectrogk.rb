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
	
	
	#CodeRunner.setup_run_class('g2') # Astrogk Inherits from GS2
	
	# This is a customised subclass of CodeRunner::Run which allows CodeRunner to submit and analyse simulations from the gyrokinetic flux tube code AstroGK, which is principally used for simulating plasmas in magnetic confinement fusion.
	#
	# It performs two distinct roles: submitting simulations and analysing the data.
	#
	# = Submitting Simulations
	#
	# This principally involves generating the input file, which is a very nontrivial task. In order to do this, it maintains a complete record of every possible input parameter for AstroGK, as well as what datatype that parameter is, and sometimes what values it is allowed to take. This allows that not only to generate the input file, but to check that the input file makes sense. However, although generating the input file works beautifully, the set of sanity checks that it makes is not exhaustive: intelligent use is still required!
	# 
	# In tandem with this, it maintains a whole set of tools for manipulating its database of input parameters. This includes updating their allowed values and also editing and accessing help for every input parameter.
	#
	# = Analysing Simulations
	#
	# The amount of analysis possible on AstroGK data is enormous, and CodeRunner hasn't got close to getting there. What it can do is:
	# 
	# * Check if the run is complete by comparing the number of completed timesteps against nstep
	# * Calculate growth rates for linear runs.
	# * Check if non-linear runs have saturated and calculate fluxes for those runs.
	# * Automatically plot a huge variety of different graphs, ranging from simple plots of heat flux versus time to three-dimensional plots of the spectrum and potential.


CodeRunner.setup_run_class('gs2', modlet: 'astrogk') # SpectroGK Inherits from AstroGK

class Gs2::Spectrogk  < Gs2::Astrogk




# Include the other files
@code_module_folder = folder = File.dirname(File.expand_path(__FILE__)) + '/spectrogk' # i.e. the directory this file is in
setup_namelists(folder)
#require folder + '/graphs.rb'
## require folder + '/namelist_tools.rb'
## require folder + '/input_file_tools.rb'
#require folder + '/gsl_data.rb'
#require folder + '/check_convergence.rb'
#require folder + '/calculations.rb'
#require folder + '/ingen.rb'
#require folder + '/properties.rb'
#require folder + '/test_gs2.rb'

NaN = GSL::NAN
# GSL::Neg





################################################
# Quantities that are calculated or determined by CodeRunner
# after the simulation has ended, i.e. quantities
# that are not available from the AstroGK output files.
################################################
# Get GS2 results and add any 
# astrogk ones (currently empty)

@results = rcp.results + [] 

###############################################
# Other useful information about the run
###############################################

@astrogk_run_info = rcp.gs2_run_info + []

@run_info = @astrogk_run_info.dup



@naming_pars = rcp.naming_pars + []

@code_long = "SpectroGK Astrophysical Full Spectral Gyrokinetic Code"

@excluded_sub_folders =[]

def input_file_header
	<<EOF
!==============================================================================
!  		SpectroGK INPUT FILE automatically generated by CodeRunner 
!==============================================================================
!
!  SpectroGK is a gyrokinetic flux tube initial value turbulence code 
!  which can be used for fusion or astrophysical plasmas. It is explicit
!  and treats velocity spectrally. 
!  
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
# Automatically generated defaults file for SpectroGK CodeRunner module#
#                                                                    #
# This defaults file specifies a set of defaults for SpectroGK which are#
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created #{Time.now.to_s}                                           #
#                                                                    #
######################################################################

@defaults_file_description = ""
EOF1
end

@source_code_subfolders = ['utils']

folder = File.dirname(File.expand_path(__FILE__)) # i.e. the directory this file is in

	end  # class SpectroGK

end # class CodeRunner


