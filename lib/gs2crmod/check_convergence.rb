class CodeRunner
class Gs2

	
def check_converged
	raise CRFatal.new("It is strongly recommended that you do not use the use_large_cache option (-U) while checking convergence. Doing so will lead to unpredictable results.") if @runner.use_large_cache
	Dir.chdir(@directory) do
		logf(:check_converged)
		return if @checked_converged and not @runner.recalc_all  

		log('@runner.class:', @runner.class)
		unless @runner.current_request == :check_converged
			@runner.requests.push :check_converged
			log 'check_converged requested recall'
			logi '@runner.requests', @runner.requests
			logi('@runner.object_id', @runner.object_id)
			return
		end	
		
		return unless @status == :Complete
		eputs @run_name
		eputs @checked_converged = true
		log("finding similar resolutions")
		@runner.generate_combined_ids(:real)
		case @grid_option
		when "box"	
			@similar_resolutions = @runner.similar_runs([:nx, :ny, :ntheta, :negrid, :naky, :ngauss, :nperiod, :delt, :jtwist], self)
		when "single"
			@similar_resolutions = @runner.similar_runs([:ntheta, :negrid, :naky, :ngauss, :nperiod], self)
		else
			raise CRFatal.new("Unknown grid option - can't get similar runs")
		end
			
		logi(@similar_resolutions)
		unless @similar_resolutions[1]
			eputs @run_name
			@converged = Feedback.get_boolean("This is is the biggest job with these params. Do you want to say it is converged?")
			return 
		end
		@similar_resolutions.sort! do |id1, id2|
			run1 = @runner.run_list[id1]
			run2 = @runner.run_list[id2]
			if @grid_option == "box" and @nonlinear_mode == "off" 
				(run1.jtwist*run1.nx*run1.negrid*run1.ngauss*run1.ntheta*run1.delt <=> run2.jtwist*run2.nx*run2.negrid*run2.ngauss*run2.ntheta*run2.delt)
			elsif @grid_option == "single" and @nonlinear_mode == "off"
				log("using nperiod: #{run1.nperiod}; #{run2.nperiod}")
				run1.negrid*run1.ngauss*run1.ntheta*run1.nperiod <=> run2.negrid*run2.ngauss*run2.ntheta*run2.nperiod

			elsif @naky	
				
				run1.nx*run1.negrid*run1.ngauss*run1.ntheta*run1.naky <=> run2.nx*run2.negrid*run2.ngauss*run2.ntheta*run2.naky
				
			else
				run1.nx*run1.negrid*run1.ngauss*run1.ntheta*run1.ny <=> run2.nx*run2.negrid*run2.ngauss*run2.ntheta*run2.ny

			end

		end

	# 	eputs @similar_resolutions
				
		log("finding my place")
		my_place = @similar_resolutions.index(@id);
	# 	eputs my_place; gets
		if my_place > 0 
			last_job = @runner.run_list[@similar_resolutions[my_place - 1]]
			unless last_job.status == :Complete
				@checked_converged = false
				return
			end
		else
			@converged = false
			return
		end

			
		log("Checking overall convergence")
		#graph = graphkit('phi2tot_vs_time_all_kys') + #last_job.graphkit('phi2tot_vs_time_all_kys')
		#graph.gnuplot
		eputs "\n \n Warning: there are no bigger jobs" unless @similar_resolutions[my_place + 1]  
		#@converged = Feedback.get_boolean("Is the plot converged?")
		#graph.close

		#(@checked_converged = true; return) unless @converged

		log("Checking convergence by ky")
		orn, last_job.runner = last_job.runner, nil
		log('last_job', last_job.pretty_inspect)
		last_job.runner = orn
# 		last_job.get_ky_graphs; last_job.get_eigenfunctions
	# 	logi(last_job.ky_graphs)
		catch(:quit_converge_check) do 
			options = {}
			list(:ky).each do |index, ky|
				options[:ky] = ky
				next if index == 1 and @grid_option == "box"
				graph = (graphkit('phi2_by_ky_vs_time', options)+last_job.graphkit('phi2_by_ky_vs_time', options))
				graph.gnuplot
				answer = Feedback.get_choice("Is the graph converged?", ["Yes", "No", "The whole run is converged, stop pestering me!"])
				graph.close
				case answer
				when /No/
					@converged = false
					throw(:quit_converge_check)
				when /stop/
					@converged = true
					throw(:quit_converge_check)
				when /Yes/
					@converged = true
				end
				cgraph = lgraph = 'efnnormmag'
				graph = (graphkit('efnnormmag', options)+last_job.graphkit('efnnormmag', options))
				
# 				graph.gnuplot

				loop do
					graph.gnuplot
					answer = Feedback.get_choice('Is the graph converged?', ['Yes', 'No', 'The whole run is converged, stop pestering me!', 'Show me the magnitude of the eigenfunctions', 'Show me the real part of the eigenfunctions again', 'Normalise the eigenfunctions', 'Denormalise the eigenfunctions', 'Reverse the axis of the current run', 'Flip the current run', 'Toggle xrange'])
					graph.close
					case answer
					when /^Yes$/		
						@converged = true
						break
					when /^No$/
						@converged = false
						throw(:quit_converge_check)
					when /stop/
						@converged = true
						throw(:quit_converge_check)
					when /magnitude/
						log 'checking convergence using magnitude'
						lgraph += 'mag'; cgraph += 'mag'
					when /Normalise/
						log 'normalising'
						lgraph += 'norm'; cgraph += 'norm'
					when /Denormalise/
						log 'denormalising'
						lgraph.gsub!(/norm/, ''); cgraph.gsub!(/norm/, '')
					when /real/
						lgraph.gsub!(/mag/, ''); cgraph.gsub!(/mag/, '')
					when /Reverse/
						cgraph = cgraph =~ /rev/ ? cgraph.sub!(/rev/, '') : cgraph + 'rev'
# 						graph = (@eigenfunctions[ky]+last_job.eigenfunctions[ky])
					when /Flip/
						cgraph = cgraph =~ /flip/ ? cgraph.sub!(/flip/, '') : cgraph + 'flip'
# 						graph = (@eigenfunctions[ky]+last_job.eigenfunctions[ky])
					when /xrange/
						if options[:range]
							options[:range] = nil
						else
							options[:range] = 0
						end
					else
						raise CRFatal.new("couldn't match choice #{answer}")
					end
					graph = graphkit(cgraph, options) + last_job.graphkit(lgraph, options)
					log graph.title
				end
				
				
			end
		end
		@checked_converged =true
		
		if last_job.checked_converged
			last_job.ky_graphs = nil
			last_job.eigenfunctions = nil
# 			last_job.t_list = nil
# 			last_job.kx_list = nil
		end
		
# 		finish_processing
	end
	ep self
end


end
end