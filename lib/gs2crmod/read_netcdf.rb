# This module reads data from the new diagnostics output file 
# <run_name>.cdf.
#
# It is intended to replace a lot of the function of gsl_data.rb
# which reads the old netcdf file. In particular, it defines a new
# generic reader function which can read any variable in the new
# netcdf file using a standard set of index constraints
#

class CodeRunner::Gs2
module ReadNetcdf

def new_netcdf_file
	if (open = @runner.run_list.keys.find_all{|id|  @runner.run_list[id].cache[:new_netcdf_file]}).size > 200
		open = open.sort_by{|id| @runner.run_list[id].cache[:new_netcdf_file_otime]}
		@runner.run_list[open[0]].new_ncclose
	end

	if cache[:new_netcdf_file] and not [:Complete, :Failed].include? @status
		new_ncclose
	end
	cache[:new_netcdf_file_otime] = Time.now.to_i
	cache[:new_netcdf_file] ||= NumRu::NetCDF.open(new_netcdf_filename)
	cache[:new_netcdf_file].sync
	cache[:new_netcdf_file]
end

def new_netcdf_filename
	@directory + '/' +  @run_name + '.cdf'
end
def new_ncclose
	cache[:new_netcdf_file].close
	cache.delete(:new_netcdf_file)
end


end #module ReadNetcdf
include ReadNetcdf

class NetcdfSmartReader
	def initialize(file)
		@file = file
	end
	def dimensions(varname)
		@file.var(varname).dims
	end
	def read_variable(varname, options)
		#start = get_start(dims, options)
		dims = dimensions(varname)
		narray = @file.var(varname).get('start'=>starts(dims, options), 'end'=>ends(dims, options))
		shape = narray.shape
		shape.delete_if{|i| i==1}
		#p 'shape', shape; STDIN.gets
		narray.reshape!(*shape)
		narray

	end
	def starts(dims, options)
		dims.map{|d| dim_start(d.name, options)}
	end
	def dim_start(name, options)
		sym = name.to_sym
		if i=options[sym + :_index]
			return i-1
		elsif i=options[sym + :_element]
			return i
		elsif i=options[sym + :min]
			return i
		else
			return 0
		end
	end
	def ends(dims, options)
		dims.map{|d| dim_end(d.name, options)}
	end
	def dim_end(name, options)
		sym = name.to_sym
		if i=options[sym + :_index]
			return i-1
		elsif i=options[sym + :_element]
			return i
		elsif i=options[sym + :max]
			return i
		else
			return -1
		end
	end
end # class NetcdfSmartReader

def netcdf_smart_reader
	NetcdfSmartReader.new(new_netcdf_file)
end

end 
