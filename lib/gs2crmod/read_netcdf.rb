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
		#p 'varname', varname
		@file.var(varname).dims
	end
	def read_variable(varname, options)
		#start = get_start(dims, options)
		if varname == 's' # dummy response for species index
			return NArray.float(@nspec||1)
		end
		dims = dimensions(varname)
		narray = @file.var(varname).get('start'=>starts(dims, options), 'end'=>ends(dims, options))
		if options[:modify_variable]
			hsh = dims.inject({}){|hsh, dim| 
				opts = options.dup
				opts[:modify_variable] = nil
				dimval = read_variable(dimension_variable_name(dim.name), opts)
				hsh[dim.name] = dimval
				hsh
			}
			narray = options[:modify_variable].call(varname, narray, hsh)
		end
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

	def axiskit(variable, options)
		GraphKit::AxisKit.autocreate(data: read_variable(variable, options), units: @file.var(variable).att('Units').get, title: @file.var(variable).att('Description').get.sub(/(,|summed|average).*$/, '').sub(/[vV]alues of (the )?/, '').sub(/ coordinate/, ''))
	end
	def dimension_variable_name(n)
		case n
		when 'X'
			'kx'
		when 'Y'
			'ky'
		when 'z'
			'kz'
		when 'e'
			'energy'
		when 'l'
			'lambda'
		when 't'
			n
		when 'm'
			'hermite'
		when 'p'
			'hankel'
		when 's'
			's'
		else
			raise "Unknown dimension #{n}"
		end
	end
	def check_no_r(non_flat_dims)
		raise "Please specify the r index for real or imaginary" if non_flat_dims.include? @file.dim('r')
	end
	def graphkit(variable, options)
		non_flat_dims=dimensions(variable).find_all{|dim| dim_start(dim.name, options) != dim_end(dim.name, options) and dim.length != 1}
		check_no_r(non_flat_dims)
		axiskits = non_flat_dims.map{|dim| dimvar = dimension_variable_name(dim.name); axiskit(dimvar, options)} + [axiskit(variable, options)]
		hash = {}
		axes = [:x, :y, :z, :f]
		axiskits.each_with_index{|ax, i| hash[axes[i]] = ax}
		kit = GraphKit.autocreate(hash)
		opts = options.dup
		opts.delete(:modify_variable)
		opts.delete(:graphkit_name)
		#kit.data[0].title += " with options: " + opts.to_s
		kit.data[0].title += " " + opts.to_s.gsub(/_(index|element)/, '')
		if kit.zlabel
			kit.zlabel = "'#{kit.zlabel}' rotate by 90"
			#kit.zlabel = nil
		end
		kit
	end
end # class NetcdfSmartReader

class OldNetcdfSmartReader < NetcdfSmartReader
	def dimension_variable_name(n)
	  if (dimnames = @file.dims.map{|dim| dim.name}).include? n
			#p ['dimnames', dimnames, n]
			n
		else
			raise "Unknown dimension #{r}: dimensions are: #{dimnames}"
		end
	end
	def check_no_r(non_flat_dims)
		raise "Please specify the ri index for real or imaginary" if non_flat_dims.include? @file.dim('ri')
	end
	def axiskit(variable, options)
		GraphKit::AxisKit.autocreate(data: read_variable(variable, options),  title: variable)
	end
end

def netcdf_smart_reader
	NetcdfSmartReader.new(new_netcdf_file)
end

def smart_graphkit(options)
	case options[:command]
	when :help
		"A smart graphkit is a direct plot of a given variable from the new netcdf file. The name of the graphkit is the name of the variable prefixed by 'cdf_'. To plot, for example, the heat flux vs time, you would give the graph name cdf_heat_flux_tot. You can use index specifiers in the the options; for example, to plot the potential as a function of kx and ky for a given time index, you would use the graph name cdf_phi2_by_mode, and the options {t_index: n}. To plot the potential as function of kx for a given ky and time would use the options {t_index, n, Y_index: m}. For each dimension you can specify the index, or a minium and/or maximum."
	when :options
		[:X_index, :Y_index, :t_index, :e_index, :l_index, :s_index, :Xmax, :Xmin, :X_element]
	else
		netcdf_smart_reader.graphkit(options[:graphkit_name].sub(/^cdf_/, ''), options)
	end
end
def old_smart_graphkit(options)
	case options[:command]
	when :help
		"An old smart graphkit is a direct plot of a given variable from the old netcdf file. The name of the graphkit is the name of the variable prefixed by 'nc_'. To plot, for example, the heat flux vs time, you would give the graph name nc_hflux_tot. You can use index specifiers in the the options; for example, to plot the potential as a function of kx and ky for a given time index, you would use the graph name nc_phi2_by_mode, and the options {t_index: n}. To plot the potential as function of kx for a given ky and time would use the options {t_index, n, ky_index: m}. For each dimension you can specify the index, or a minium and/or maximum."
	when :options
		[:kx_index, :ky_index, :t_index, :e_index, :l_index, :s_index, :kxmax, :kxmin, :kx_element]
	else
	 vars = OldNetcdfSmartReader.new(netcdf_file).graphkit(options[:graphkit_name].sub(/^nc_/, ''), options)
	end
end

def hyperviscosity_graphkit(options)
	raise "This only works for spectrogk"  unless spectrogk?
	options[:modify_variable] = Proc.new do |varname, narray, dimhash|
		#dimnames = dimhash.keys
		shape = narray.shape
		if  varname == "gnew2_ta"
			#p dimhash
			#p dimhash['Y']
			ky = dimhash['Y'].to_a.to_gslv
			kx = dimhash['X'].to_a.to_gslv.to_box_order
			shape = narray.shape
			for ig in 0...shape[0]
				for it in 0...shape[1]
					for ik in 0...shape[2]
						for il in 0...shape[3]
							for ie in 0...shape[4]
								for is in 0...shape[5]
									narray[ig,it,ik,il,ie,is]*=(ky[ik]**2.0 + kx[it]**2.0)**(2*@nexp)*@d_hypervisc
								end
							end
						end
					end
				end
			end
		end
		narray
	end
	options[:graphkit_name] = 'cdf_gnew2_ta'
	kit = smart_graphkit(options)
end
def hypercoll_graphkit(options)
	raise "This only works for spectrogk"  unless spectrogk?
	options[:modify_variable] = Proc.new do |varname, narray, dimhash|
		#dimnames = dimhash.keys
		p varname, dimhash
		if  varname == "gnew2_ta"
			shape = narray.shape
			m = dimhash['m']
			mmax = new_netcdf_file.var('hermite').get.to_a.size - 1
			p 'shape',shape
			for ig in 0...shape[0]
				for it in 0...shape[1]
					for ik in 0...shape[2]
						for il in 0...shape[3]
							for ie in 0...shape[4]
								for is in 0...shape[5]
									narray[ig,it,ik,il,ie,is]*=send(:nu_h_ + (is+1).to_sym)*(m[il]/mmax)**send(:nexp_h_ + (is+1).to_sym)
								end
							end
						end
					end
				end
			end
		end
		narray
	end
	options[:graphkit_name] = 'cdf_gnew2_ta'
	kit = smart_graphkit(options)
end
def lenardbern_graphkit(options)
	raise "This only works for spectrogk"  unless spectrogk?
	options[:modify_variable] = Proc.new do |varname, narray, dimhash|
		#dimnames = dimhash.keys
		if  varname == "gnew2_ta"
			m = dimhash['m']
			shape = narray.shape
			for ig in 0...shape[0]
				for it in 0...shape[1]
					for ik in 0...shape[2]
						for il in 0...shape[3]
							for ie in 0...shape[4]
								for is in 0...shape[5]
									narray[ig,it,ik,il,ie,is]*=send(:nu_ + (is+1).to_sym)*m[il]
								end
							end
						end
					end
				end
			end
		end
		narray
	end
	options[:graphkit_name] = 'cdf_gnew2_ta'
	kit = smart_graphkit(options)
end


end 
