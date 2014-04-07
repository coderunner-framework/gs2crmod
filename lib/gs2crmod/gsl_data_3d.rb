require 'narray'
class NArray
	 def _dump *ignored
		 Marshal.dump :typecode => typecode, :shape => shape, :data => to_s
	 end
	 def self._load buf
		 h = Marshal.load buf
		 typecode = h[:typecode]
		 shape = h[:shape]
		 data = h[:data]
		 to_na data, typecode, *shape
	 end
	 def inspect
		 #ep "called inspect"
		 "#{self.class}.to_narray(#{self.to_a.inspect})"
	 end
end

	class GSL::Tensor 
		class << self
			def method_missing(meth, *args)
				#ep 'calling... ', meth, args
				ans = new(NArray.send(meth, *args.reverse))
				#ep 'got', ans
				ans
			end
		end
		attr_reader :narray
		def self.alloc(*args)
			new(NArray.float(*args.reverse))
		end
		def initialize(narray)
			@narray = narray
		end
		def inspect
			"GSL::Tensor.new(#{@narray.inspect})"
		end
		def [](*args)
			#if args.inject(true){|b,i| b and i.kind_of? Integer}
				#@narray[*args.reverse]
			#else
				#self.class.new(@narray[*args.reverse])
			#end
		 	case 	ans =	@narray[*args.reverse]
			when Numeric
				ans
			else
				self.class.new(@narray[*args.reverse])
			end
				
		end
		def []=(*args, value)
		#def []=(*args)
			#ep 'args', args, value
			@narray[*args.reverse] = value
		end
		def shape
			@narray.shape.reverse
		end
		def reshape!(*args)
			#ep 'rags', args
			@narray.reshape!(*args.reverse)
		end
		def to_a
			@narray.transpose(*(0...@narray.shape.size).to_a.reverse).to_a
		end
		def transpose(*args)
			self.class.new(@narray.transpose(*args))
		end
		def method_missing(meth, *args)
			result = @narray.send(meth, *args.reverse)
			if result.kind_of? NArray
				self.class.new(result)
			else
				result
			end
		rescue NoMethodError
			self.class.new(NMath.send(meth, @narray))
		end
		def iterate(&block)
			shp = shape
			cumul = 1
			cumulshp = []
			for i in 1..shape.size
				cumulshp[shp.size-i] = cumul
				cumul *= shp[shp.size-i]
			end
			 #= shape.reverse.map{|dim| cumul*=(dim); cumul.to_i}.reverse
			#ep cumulshp; gets
			(cumulshp[0]*shp[0]).times do |n|
				#indexes = cumulshp.reverse.map{|cumul| rem = n%cumul; n-=rem; rem}.reverse
				indexes = cumulshp.map{|cumul| idx = (n/cumul).floor; n -= idx*cumul; idx}
				yield(*indexes)
			end
		end
		def iterate_row_maj(&block)
			shp = shape
			cumul = 1
			cumulshp = []
			for i in 0...shape.size
				cumulshp[i] = cumul
				cumul *= shp[i]
			end
			 #= shape.reverse.map{|dim| cumul*=(dim); cumul.to_i}.reverse
			#ep cumulshp; gets
			(cumulshp[-1]*shp[-1]).times do |n|
				#indexes = cumulshp.reverse.map{|cumul| rem = n%cumul; n-=rem; rem}.reverse
				indexes = cumulshp.reverse.map{|cumul| idx = (n/cumul).floor; n -= idx*cumul; idx}.reverse
				yield(*indexes)
			end
		end



	end
	class GSL::TensorComplex < GSL::Tensor
		attr_reader :narray
		def self.alloc(*args)
			new(NArray.complex(*args.reverse))
		end
		def real
			GSL::Tensor.new(@narray.real)
		end
		def abs
			GSL::Tensor.new(@narray.abs)
		end
	end
class CodeRunner::Gs2
	


	def gsl_tensor(name, options)
		tensor = send((name.to_s+"_gsl_tensor").to_sym , options)
	end
	module GSLTensors
		def moment_gsl_tensor(options)
			if options[:t_index]
				raise ArgumentError.new("Moments are not written out as a function of time currently")
				#ep options; gets
				raise CRFatal.new("write_phi_over_time is not enabled so this function won't work") unless @write_phi_over_time
				arr =  GSL::Tensor.new(netcdf_file.var(options[:field_name].to_s + '_t').get({'start' => [0,(options[:thetamin]||0),0,0, options[:t_index] - 1], 'end' => [-1,(options[:thetamax]||-1),(options[:nakx]||0)-1,(options[:naky]||0)-1, options[:t_index] - 1]}))
				#ep 'arr.shape', arr.shape
				arr.reshape!(*arr.shape.slice(1...arr.shape.size))
				
			else
				arr =  GSL::Tensor.new(netcdf_file.var(options[:moment_name]).get({'start' => [0,(options[:thetamin]||0),0,0,options[:species_element]], 'end' => [-1,(options[:thetamax]||-1),(options[:nakx]||0)-1,(options[:naky]||0)-1,options[:species_element]]}))
				#ep 'arr.shape', arr.shape
			end
			arr.reshape!(*arr.shape.slice(1...arr.shape.size))
			arr[0, true, true, true] = 0.0 if options[:no_zonal]
			#arr = arr[options[:nakx] ? 0...options[:nakx] : true, options[:naky] ? 0...options[:naky] : true, true, true] if options[:nakx] or options[:naky]
			return arr

		end
		def field_netcdf_name(field_name, time_varying = false)
			#p field_name.to_s
			name =  case field_name.to_s
							when /phi/
								time_varying ? 'phi_t' : 'phi'
							when /density/
								time_varying ? 'ntot_t' : 'density'
							when /apar/
								time_varying ? 'apar_t' : 'apar'
							else
								raise "Unknown field name: #{field_name}"
							end
			#p name
			return name
		end
		def field_species_element(options)
			case options[:field_name].to_s
			when /density/
				options.convert_to_index(self, :species)
				#ep 'options', options
				options[:species_index] - 1
			else
				nil
			end
		end
		def field_gsl_tensor(options)
			species_element = field_species_element(options)
			#ep 'species_element', species_element
			if options[:t_index]
				#ep options; gets
                #raise CRFatal.new("write_phi_over_time is not enabled so this function won't work") unless @write_phi_over_time
				arr =  GSL::Tensor.new(netcdf_file.var(field_netcdf_name(options[:field_name], true)).get({'start' => [0,(options[:thetamin]||0),0,0, species_element, options[:t_index] - 1].compact, 'end' => [-1,(options[:thetamax]||-1),(options[:nakx]||0)-1,(options[:naky]||0)-1, species_element, options[:t_index] - 1].compact}))
				#ep 'arr.shape', arr.shape
				arr.reshape!(*arr.shape.slice(1...arr.shape.size))
				
			else
				arr =  GSL::Tensor.new(netcdf_file.var(field_netcdf_name(options[:field_name])).get({'start' => [0,(options[:thetamin]||0),0,0, species_element].compact, 'end' => [-1,(options[:thetamax]||-1),(options[:nakx]||0)-1,(options[:naky]||0)-1, species_element].compact}))
				#ep 'arr.shape', arr.shape
			end
			if species_element
				arr.reshape!(*arr.shape.slice(1...arr.shape.size))
			end
			if options[:interpolate_x]
				shape = arr.narray.shape
				#p 'shape', shape
				shape[2] = (shape[2]-1)*options[:interpolate_x] + 1
				#p shape
				arr = GSL::Tensor.new(arr.narray.expand(*shape, 0.0))
			end
			if options[:interpolate_y]
				shape = arr.narray.shape
				#p 'shape', shape
				shape[3] = (shape[3]-1)*options[:interpolate_y] + 1
				#p shape
				arr = GSL::Tensor.new(arr.narray.expand(*shape, 0.0))
			end

			if gryfx? and options[:periodic]
				shape = arr.narray.shape
				shape[1]+=1
				arr = GSL::Tensor.new(arr.narray.expand(*shape, 0.0))
				shpe = arr.shape
				for i in 0...shpe[0]
					for j in 0...shpe[1]
						for r in 0...shpe[3]
							arr[i, j, -1, r] = arr[i, j, 0, r]
						end
					end
				end

			end

			arr[0, true, true, true] = 0.0 if options[:no_zonal]
			#arr = arr[options[:nakx] ? 0...options[:nakx] : true, options[:naky] ? 0...options[:naky] : true, true, true] if options[:nakx] or options[:naky]
			return arr

		end

		# Returns a rank 3 tensor which is the real potential (i.e. Fourier transformed from the GS2 output) as a function of the y index, the x index and the theta index.

		def phi_real_space_gsl_tensor(options)
			return field_real_space_gsl_tensor(options.absorb(field_name: :phi))
		end
																				
		def field_real_space_gsl_tensor(options)
			fieldc = field_gsl_tensor_complex(options)
			shape = fieldc.shape
			workspacex = GSL::Vector::Complex.alloc(shape[1])
			workspacey = GSL::Vector.alloc(shape[0]*2-2+shape[0]%2)
			field_real_space = GSL::Tensor.alloc(workspacey.size, shape[1], shape[2])
			for j in 0...shape[2] #theta
				for i in 0...shape[0] #ky
					#narr = fieldc[i, true, j] 
					for k in 0...shape[1]
						workspacex[k] = GSL::Complex.alloc(fieldc[i,k,j].real, fieldc[i,k,j].imag)
					end
					workspacex = workspacex.backward
					for k in 0...shape[1]
						fieldc[i,k,j] = Complex(*workspacex[k].to_a)
					end
				end
				for k in 0...shape[1] #kx
					m = 0
					for i in 0...shape[0] #ky
						workspacey[m] = fieldc[i,k,j].real
						m+=1
						next if i==0 or (shape[0]%2==0 and i == shape[0]/2 + 1)
						workspacey[m] = fieldc[i,k,j].imag
						m+=1
					end
					workspacey = workspacey.backward
					for i in 0...workspacey.size
						field_real_space[i,k,j] = workspacey[i]
					end
				end
			end
			shp = field_real_space.shape
			#ep options
			field_real_space = field_real_space[options[:ymin]||0..options[:ymax]||(shp[0]-1), options[:xmin]||0..options[:xmax]||(shp[1]-1), true] 
			if kint = options[:interpolate_theta]
				shape = field_real_space.shape
				new_shape = shape.dup
				new_shape[-1] = ((shape[-1]-1)*kint+1)
				field_real_space_new = GSL::Tensor.float(*new_shape)
				#p shape,new_shape
				for i in 0...(new_shape[0])
				for j in 0...(new_shape[1])
				field_real_space_new[i,j, new_shape[-1]-1] = field_real_space[i,j,shape[-1]-1] # set the endpoint
				for k in 0...(new_shape[-1]-1)
					km = k%kint
					frac = km.to_f/kint.to_f
					#kold = (k-km)/(new_shape[-1]-1)*(shape[-1]-1)
					kold = (k-km)/kint
					#ep ['k', k, 'kold', kold]
					field_real_space_new[i,j, k] = field_real_space[i,j, kold] * (1.0-frac) + field_real_space[i,j, kold+1] * frac
				end
				end
				end
				field_real_space = field_real_space_new
			end	

			return field_real_space

		end
		def field_real_space_gsl_tensor_2(options)
			field = field_gsl_tensor(options)
			field_narray = field.narray
			shape = field.shape
			workspacex = GSL::Vector::Complex.alloc(shape[1])
			workspacey = GSL::Vector.alloc(shape[0]*2-2+shape[0]%2)
			field_real_space = GSL::Tensor.alloc(workspacey.size, shape[1], shape[2])
			field_real_space_narray = field_real_space.narray
			for j in 0...shape[2] #theta
				for i in 0...shape[0] #ky
					#narr = fieldc[i, true, j] 
					for k in 0...shape[1]
						workspacex[k] = GSL::Complex.alloc(field_narray[0,j,k,i], field_narray[1,j,k,i])
					end
					workspacex = workspacex.backward
					for k in 0...shape[1]
						field_narray[0,j,k,i] = workspacex[k].real
						field_narray[1,j,k,i] = workspacex[k].imag
					end
				end
				for k in 0...shape[1] #kx
					m = 0
					for i in 0...shape[0] #ky
						workspacey[m] = field_narray[0,j,k,i]
						m+=1
						next if i==0 or (shape[0]%2==0 and i == shape[0]/2 + 1)
						workspacey[m] = field_narray[1,j,k,i]
						m+=1
					end
					workspacey = workspacey.backward
					for i in 0...workspacey.size
						field_real_space_narray[j,k,i] = workspacey[i]
					end
				end
			end
			shp = field_real_space.shape
			#p 'test', field_real_space[0,2,3]
			#ep options
			field_real_space = field_real_space[options[:ymin]||0..options[:ymax]||(shp[0]-1), options[:xmin]||0..options[:xmax]||(shp[1]-1), true] 
			#p 'test2', field_real_space[0,2,3]
			if kint = options[:interpolate_theta]
				shape = field_real_space.shape
				new_shape = shape.dup
				new_shape[-1] = ((shape[-1]-1)*kint+1)
				field_real_space_new = GSL::Tensor.float(*new_shape)
				field_real_space_new_narray = field_real_space_new.narray
				#p shape,new_shape
				for i in 0...(new_shape[0])
				for j in 0...(new_shape[1])
				field_real_space_new_narray[new_shape[-1]-1, j, i] = field_real_space_narray[shape[-1]-1, j, i] # set the endpoint
				for k in 0...(new_shape[-1]-1)
					km = k%kint
					frac = km.to_f._orig_div(kint.to_f)
					#kold = (k-km)/(new_shape[-1]-1)*(shape[-1]-1)
					kold = (k-km)._orig_div(kint)
					#ep ['k', k, 'kold', kold]
					field_real_space_new_narray[k,j,i] = field_real_space_narray[kold,j,i]._orig_mul(1.0-frac) + field_real_space_narray[kold+1,j,i]._orig_mul(frac)
					#if (i==0 and j==2 and k==3)
						#p ['frac', frac]
					#end
				end
				end
				end
				field_real_space = field_real_space_new
			end	
			#p field_real_space_new.shape;

			return field_real_space

		end
		def apar_gsl_tensor(options)
			return GSL::Tensor.new(netcdf_file.var('apar').get)
		end
		def bpar_gsl_tensor(options)
			return GSL::Tensor.new(netcdf_file.var('bpar').get)
		end
		# Order is R0,Z0,a0,Rprim,Zprim,aprim
		def geometric_factors_gsl_tensor(options)
			#ops = options.dup; ops.delete :phi
		#ep ops; gets
			case @equilibrium_option
			when "s-alpha"
				return geometric_factors_salpha_gsl_tensor(options)
			else
				theta_vec = gsl_vector(:theta, options)
				factors = GSL::Tensor.alloc(6,theta_vec.size)
				values = File.read("#@directory/#@run_name.g").split(/\s*\n\s*/)
				3.times{values.shift}
				values = values.map{|str| str.split(/\s+/).map{|s| s.to_f}}.transpose
				#ep values
				shape = factors.shape
				for i in 0...shape[0]
						unless options[:interpolate_theta]
							for j in 0...shape[1]
								factors[i,j] = values[i+1][j]
							end
						else
							opts = options.dup
							opts[:interpolate_theta] = nil
							theta_vec_short = gsl_vector(:theta, {})
							#p 'sizes', [theta_vec_short.size, values[i+1].to_gslv.size]
							interp = GSL::ScatterInterp.alloc(:linear, [theta_vec_short, values[i+1].to_gslv], true)
							for j in 0...theta_vec.size
								factors[i,j] = interp.eval(theta_vec[j])
							end
						end
				end
				#ep factors
				return factors
			end
		end
		# Order is R0,Z0,a0,Rprim,Zprim,aprim
		def geometric_factors_salpha_gsl_tensor(options)
			raise "Please specify options[:Rgeo]" unless options[:Rgeo]
			theta_vec = gsl_vector(:theta, options)
			factors = GSL::Tensor.alloc(6,theta_vec.size)
			q_actual = options[:q_actual]
			pka = epsl/q_actual
			#pka = @pk||2*@kp
			for i in 0...theta_vec.size
				theta = theta_vec[i]
				c = Math.cos(theta)
				s = Math.sin(theta)
				factors[0,i] = options[:Rgeo]*(1.0 + eps * c + eps * (shift||0))
				factors[1,i] = options[:Rgeo] * eps * s
				factors[2,i] = -epsl/pka * theta - eps * epsl/pka * s
				factors[3,i] = c * options[:Rgeo] * eps / 2
				factors[4,i] = s * options[:Rgeo] * eps / 2 
				factors[5,i] = - theta * epsl**2 / 2 / pka / eps * (shat||0) - epsl**2 / 2 / pka * s
			end
			return factors
		end
		private :geometric_factors_salpha_gsl_tensor
		
		# Returns a rank 2 tensor, which gives, as a function of the x index j and the theta index k, the y index nearest to a poloidal plane at angle options[:torphi] is the torus was filled with periodic copies of the flux surface. Used for making cross sections at a constant toroidal angle.

		def constant_torphi_surface_gsl_tensor(options)
			ops = options.dup
			IRRELEVANT_INDICES.each{|v| ops.delete(v)}
			return  cache[[:constant_torphi_surface_gsl_tensor, ops]] if cache[[:constant_torphi_surface_gsl_tensor, ops]]
			correct_3d_options(options)
			torphiout = options[:torphi]
			cyls = cylindrical_coordinates_gsl_tensor(options.absorb({extra_points: :y}))
			shpc = cyls.shape
			factors = geometric_factors_gsl_tensor(options)
			#ep shpc, 'shpc'
			#xsize = case shpc[2]

			yvec = gsl_vector('y', options)
			#ep yvec.to_a ; gets
			x = gsl_vector('x', options)
			dy = yvec[1] - yvec[0]
			torphi_const = GSL::Tensor.int(shpc[2], shpc[3]) # don't include extra x point
			xfac = 1.0 / options[:rho_star_actual]
			yfac = options[:rhoc_actual] / options[:q_actual] / options[:rho_star_actual]	
						#coordinates[2,i,j,k] = y[i] / yfac - factors[2,k] - x[j]/xfac*factors[5,k] # phi
			twopi = Math::PI*2
			for j in 0...shpc[2]
				for k in 0...shpc[3]
					y = yfac * (torphiout + factors[2,k] + x[j]/xfac*factors[5,k])
					if options[:no_copies]
						i = (y/dy).floor 
					else
						i = (y/dy).floor % yvec.size
					end
					torphi_const[j,k] = i
				end
			end
			return torphi_const

			#ep torphi_const; gets
		end
		#def constant_torphi_surface_gsl_tensor2(options)
			#ops = options.dup
			#IRRELEVANT_INDICES.each{|v| ops.delete(v)}
			#return  cache[[:constant_torphi_surface_gsl_tensor, ops]] if cache[[:constant_torphi_surface_gsl_tensor, ops]]
			#torphiout = options[:torphi]
			#correct_3d_options(options)
			#cyls = cylindrical_coordinates_gsl_tensor(options.absorb({extra_points: :y}))
			#shpc = cyls.shape
			##ep shpc, 'shpc'
			##xsize = case shpc[2]

			#torphi_const = GSL::Tensor.int(shpc[2], shpc[3]) # don't include extra x point
			##ep torphi_const; gets
			#y = gsl_vector('y', options)
			#lastbracketed = nil
			#lastj = -1
			#for k in 0...shpc[3] # theta index
			#for j in 0...(shpc[2]) # x index
				#deltorphi = cyls[2,shpc[1]-1,j,k] - cyls[2,0,j,k]
				#raise "Periodicity not satisfied: #{(2*Math::PI/deltorphi+1.0e-8)%1.0}, #{(2*Math::PI/deltorphi+1.0e-5)}" unless ((2.0*Math::PI/deltorphi)+1.0e-8)%1.0 < 1.0e-5
				#m3 = (torphiout)%deltorphi
			#for i in 0...shpc[1]-1 # y index, excluding periodic extra point
				##p i
				#torphi1 = cyls[2,i,j,k]
				#torphi2 = cyls[2,i+1,j,k]
				##ep cyls[2,true,j,k].to_a
				#m1 = (torphi1 )%deltorphi 
				#m2 = (torphi2 )%deltorphi 
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
					##p 'measure', measure,  [torphi1, torphi2, upp, lwr , i,j,k, y[i], y[(i+1)%y.size], deltorphi, n, a1, a2, a3, b1, b2, b3] if measure and j==0 #; gets if [j,k] == [5,8] # if measure
					#b or measure
				#end
				##bracketed = bracketed2
				#raise "Measures don't agree #{bracketed}, #{bracketed2}" unless bracketed2 == bracketed


				##d2 = torphi2 - torphiout
				##d1 = torphiout - torphi1
				#if bracketed
					#raise "Doubled up" if lastbracketed == [j,k]
					#raise "Missed: #{j},#{k}, #{lastbracketed.inspect} #{[j-1,k].inspect} #{[shpc[2]-1, k-1].inspect}" unless lastbracketed == [j-1,k] or lastbracketed == [shpc[2]-1, k-1] if lastbracketed
					#torphi_const[j,k] = i
					#lastbracketed = [j,k]
					#lastj = j
				#end
			#end # y loop
			#end # x loop
			#end # theta loop
			##torphi_const2 = constant_torphi_surface_gsl_tensor2(options)
			##p 'the same? ', torphi_const2.to_a, torphi_const.to_a, torphi_const2 == torphi_const
			##exit
			##exit
			#cache[[:constant_torphi_surface_gsl_tensor, ops]] = torphi_const
			  ## save the run to save the hard_cache
			#return torphi_const
		#end

		FIELD_VALUES = [:phi, :density, :apar, :bpar]
		TRIVIAL_INDICES = [:graphkit_name]
		TIME_VARYING_INDICES = [:t_index, :begin_element, :end_element, :frame_index, :t_index_window]
		IRRELEVANT_INDICES  = FIELD_VALUES + TRIVIAL_INDICES + TIME_VARYING_INDICES
		
		# Adjust n0, rho_star_actual and q_actual to ensure periodicity
		#
		def correct_3d_options(options)
			raise "Please specify options[:rho_star] or options[:n0]" unless options[:rho_star] or options[:n0]
			case @equilibrium_option
			when "s-alpha"
				qinp = epsl / (pk||2*kp)
				#xfac = @epsl**4/options[:rho_star]/4/pka**2/@eps**2
				#xfac_geo = 1
				#yfac = 1/options[:rho_star]/@epsl*2*pka*@eps
				#yfac_geo = 2*pka*@eps/@epsl**2
				#yfac_geo = 2*pka*@eps/@epsl**2
				options[:rhoc_actual] =rhoc =  2 * eps / epsl
			else
				options[:rhoc_actual] = rhoc = @rhoc
				qinp = @qinp
			end
			#eputs "Checking that rho_star and q satisfy periodicity..."
			rho_star_inp = options[:rho_star]
			y = gsl_vector('y', options)
			ly = (y[1]-y[0]) * (y.size) 
			n0_fac = 2.0*Math::PI * rhoc / ly
			n0_inp = options[:n0] || n0_fac / qinp / rho_star_inp 
			if n0_inp%1.0==0.0
				n0 = n0_inp
			else
				#eputs "Input n0 is equal to #{n0_inp}..."
				n0 = n0_inp.ceil
				#eputs "Set n0 to #{n0}..."
			end
			
			if (qinp*n0)%1.0==0.0
				q_actual = qinp
			else
				q_actual = (qinp*n0).round.to_f/n0
				#eputs "Set q to #{q_actual}..."
			end
			options[:q_actual] = q_actual
			unless options[:rho_star_actual] and options[:rho_star_actual] == n0_fac/n0/q_actual
				#eputs "Adjusting rho_star to satisfy periodicity ..."
				options[:rho_star_actual] = n0_fac/n0/q_actual
				#eputs "Set rhostar to #{options[:rho_star_actual]}..."
				#eputs "Note... to avoid adjustment of q specify n0 as an input rather than rho_star. Make sure that n0 is an integer and n0 * q is an integer."
			end
		end
		

		# Return a rank 4 tensor which give
		# cylindrical coordinates R,Z,torphi as a function
		# of gs2 coordinates y, x, theta.
		#
		# 	a = cylindrical_coordinates_gsl_tensor(options)
		#
		# 	# pseudocode
		# 	R(y[i], x[j], theta[k]) = a[0,i,j,k]
		# 	Z(y[i], x[j], theta[k]) = a[1,i,j,k]
		# 	torphi(y[i], x[j], theta[k]) = a[2,i,j,k]
		def cylindrical_coordinates_gsl_tensor(options)
			ops = options.dup
			(IRRELEVANT_INDICES + [:torphi, :torphi_values]).each{|v| ops.delete(v)}
			return  cache[[:cylindrical_coordinates_gsl_tensor, ops]] if cache[[:cylindrical_coordinates_gsl_tensor, ops]]
			#ep ops; gets
			#options = options.dup
			x = gsl_vector('x', options)
			y = gsl_vector('y', options)
			ly = 2*Math::PI*y0#(y[1]-y[0]) * (y.size) 
			if [true,:x].include? options[:extra_points]
				ep "Extending x..."
				x = x.connect([2*x[-1] - x[-2]].to_gslv).dup
			end
			if [true,:y].include? options[:extra_points]
				ep "Extending y..."
				y = y.connect([2.0*y[-1] - y[-2]].to_gslv).dup
				raise "ly corrected incorrectly #{ly},#{y[-1]},#{y[0]},#{y[-1]-y[0]}" unless (ly-(y[-1] - y[0])).abs / ly.abs < 1.0e-6
			end


			#if options[:xmax] 
			 #if	options[:xmin]
				 #x = x.subvector(options[:xmin], options[:xmax] - options[:xmin])
			 #else
				 #x = x[options[:xmax]].to_gslv
			 #end
			#elsif options[:xmin]
			 #x = x[options[:xmin]].to_gslv
			#end
			#if options[:ymax] 
			 #if	options[:ymin]
				 #y = y.subvector(options[:ymin], options[:ymax] - options[:ymin])
			 #else
				 #y = y[options[:ymax]].to_gslv
			 #end
			#elsif options[:ymin]
			 #y = y[options[:ymin]].to_gslv
			#end



			#ep [options, options[:xmin]||0, (options[:xmax]||x.size-1) - (options[:xmin]||0) + 1]
			x = x.subvector(options[:xmin]||0, (options[:xmax]||x.size-1) - (options[:xmin]||0) + 1).dup # if options[:xout] and options[:xin]
			y = y.subvector(options[:ymin]||0, (options[:ymax]||y.size-1) - (options[:ymin]||0) + 1).dup # if options[:yout] and options[:yin]
			###y = y.subvector(options[:ymin], options[:ymax] - options[:ymin] + 1)# if yi = options[:yout] and options[:yin]
			#	
			###ep 'ncopy', options[:ncopy]
			#y = y + options[:ncopy] * (y[-1]-y[0]) if options[:ncopy]
			y = y + options[:ncopy] * ly if options[:ncopy]
			#ep 'y', y
				#ep y; gets
			#ep options; gets
			theta = gsl_vector('theta', options)
			#ep theta; gets;
			#ep 'thsize', @ntheta, theta.size
			correct_3d_options(options)
			rhoc = options[:rhoc_actual]
			q_actual = options[:q_actual]
			xfac = 1.0 / options[:rho_star_actual]
			yfac = rhoc / q_actual / options[:rho_star_actual]
			factors = geometric_factors_gsl_tensor(options)
			
			#ep ['factors.shape', factors.shape]




			coordinates = GSL::Tensor.alloc(3, y.size, x.size, theta.size)
			for i in 0...y.size
				for j in 0...x.size
					for k in 0...theta.size
						coordinates[0,i,j,k] = factors[0,k] + x[j]/xfac*factors[3,k] # R
						coordinates[1,i,j,k] = factors[1,k] + x[j]/xfac*factors[4,k] # Z
						coordinates[2,i,j,k] = y[i] / yfac - factors[2,k] - x[j]/xfac*factors[5,k] # phi
						#ep [i,j,k], coordinates[0, false, j,k].to_a
						if gs2f = options[:gs2_coordinate_factor]
							rgs2 = (x[j]**2 + y[i]**2 )**0.5*(1.0 + 2.0 * Float::EPSILON)
							#p ['x', x[j], 'y', y[i], 'r', rgs2] if agk?
							if rgs2 < 1.0e-8
								phigs2 = 0
							else
								phigs2 = Math.acos(x[j]/rgs2)
							end
							coordinates[0,i,j,k] = rgs2 * gs2f + coordinates[0,i,j,k] * (1.0-gs2f)
							coordinates[1,i,j,k] = theta[k] * gs2f + coordinates[1,i,j,k] * (1.0-gs2f)
							coordinates[2,i,j,k] = phigs2 * gs2f + coordinates[2,i,j,k] * (1.0-gs2f)
						end


							
					end
				end
			end
			#exit
			case tp = options[:toroidal_projection]
			when Numeric
				coordinates[2, false] = tp
			end
			cache[[:cylindrical_coordinates_gsl_tensor, ops]] = coordinates
			#save  # save the run to save the hard_cache
			return coordinates
		end
	
		# Return a rank 4 tensor which give
		# cartesian coordinates X,Y,Z as a function
		# of gs2 coordinates y, x, theta.
		#
		# 	a = cartesian_coordinates_gsl_tensor(options)
		#
		# 	# pseudocode
		# 	X(y[i], x[j], theta[k]) = a[0,i,j,k]
		# 	Y(y[i], x[j], theta[k]) = a[1,i,j,k]
		# 	Z(y[i], x[j], theta[k]) = a[2,i,j,k]


		def cartesian_coordinates_gsl_tensor(options)
			cyl = cylindrical_coordinates_gsl_tensor(options)
			shape = cyl.shape
			cart = GSL::Tensor.alloc(*shape)
			for i in 0...shape[1]
				for j in 0...shape[2]
					for k in 0...shape[3]
						r = cyl[0,i,j,k]
						z = cyl[1,i,j,k]
						phi = cyl[2,i,j,k]
						#cart[0,i,j,k] = r # Y
						cart[0,i,j,k] = r*Math.cos(phi) # X
						#cart[1,i,j,k] = phi # X
						cart[1,i,j,k] = r*Math.sin(phi) # Y
						cart[2,i,j,k] = z
					end
				end
			end
			cart
		end


			




		

	end #module
	include GSLTensors
	module GSLComplexTensors
		def phi_gsl_tensor_complex(options)
			return field_gsl_tensor_complex(options.absorb({field_name: :phi}))
		end
		def field_gsl_tensor_complex(options)
			field = field_gsl_tensor(options)
			fieldc = GSL::TensorComplex.alloc(*field.shape.slice(0..2))
			nac = fieldc.narray
			na = field.narray
			for i in 0...field.shape[0]
				for j in 0...field.shape[1]
					for k in 0...field.shape[2]
						nac[k,j,i] = Complex(na[0,k,j,i],na[1,k,j,i])
					end
				end
			end
			return fieldc	
		end
	end #module
	include GSLComplexTensors
end #class
