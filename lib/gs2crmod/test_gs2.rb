class CodeRunner
	class Gs2
		
		# See TestGs2
		
		def self.test_gs2(*args)
			TestGs2.test_gs2(*args)
		end
		
		# = Gs2 Test Suite
		# This class is designed to run a set of functional tests to test the GS2 source. It is not a test suite for GS2crmod, the GS2 CodeRunner module. It takes all the input files from the folder <tt>test_cases</tt>, runs them and checks that they give the same output as the standard case. Only the NetCDF file from the standard case should be included, and should have the same name as the test case input file.
		
		module TestGs2 
			
		
			# The folder where all the test files are stored.
			
			TEST_FOLDER = Dir.pwd #File.dirname(__FILE__) + "../test_cases"
			
			
			# Give the standardised name of the test_results folder, given the name of the test case input file
			
			def self.results_folder(input_file)
				'test_results/' + input_file.sub(/\.in$/, '')
			end
			
			def self.test_case_folder(input_file)
				TEST_FOLDER + '/' + input_file.sub(/\.in$/, '') 
			end
			
			# The io object all test info is written to.
			
			TEST_OUT = STDERR
			
			# Run all the test cases. test_case_location should be the folder test_cases in the gs2 source.
			#
			# Options are
			# * restart (default: true) Delete all results and start again.
			# * submit (default: true) Submit any runs that haven't been submitted.
			# * serial (default: true) Wait till one test run has finished before starting another
			
			def self.test_gs2(test_case_location, options={})
					Gs2.send(:include, self)
					test_case_location.sub!(/~/, ENV['HOME'])
					TEST_FOLDER.gsub!(/^.*$/, test_case_location)
					raise "The first argument should be the test_cases folder" unless File.basename(TEST_FOLDER) == "test_cases"
					@tests_failed = {}

					
					serial = true unless ["false", false].include? options[:serial]
					restart = true unless ["false", false].include? options[:restart]
					submit = true unless ["false", false].include? options[:submit]
					
					
					@test_runner = CodeRunner.fetch_runner(Y: Dir.pwd, u: true)
					@submitted_tests = @test_runner.runs.map{|run| run.run_name}
					test_cases = Dir.entries(TEST_FOLDER).find_all do |entry|
# 						p entry
						File.directory?(TEST_FOLDER + '/' + entry) and not entry =~ /^\./
					end
# 					p test_cases
					
					if restart
						FileUtils.rm_r 'test_results' if FileTest.exist? 'test_results' 
						FileUtils.makedirs 'test_results'
					end
					
					# Submit the tests
					if submit
						test_cases.each do |test_case|
							next if FileTest.exist? results_folder(test_case)
							submit_run(test_case)
							run_tests if serial
						end
					end
					run_tests
					
					if @tests_failed.size == 0
						TEST_OUT.puts "All Tests Completed Successfully"
						eputs "All Tests Completed Successfully" unless TEST_OUT == STDERR
					else
						eputs "Tests Were Failed" 
						@tests_failed.each do |name, hash|
							sep = "----------------------------------------------"
							TEST_OUT.puts '', sep, "                 Test Failed", sep
							TEST_OUT.puts "Name: #{name}"
							TEST_OUT.puts "Description: \n#{hash[:description]}", ''
							TEST_OUT.puts "This test failed on these checks: #{hash[:checks_failed]}", ''
							TEST_OUT.puts sep, ''
						end
					end
								
			end
			
			# Check to see which run has completed and run tests for those runs.
			
			def self.run_tests
				eputs "Waiting for runs to complete"
				loop do
						break if @submitted_tests.size ==0
						@test_runner.update(false)	
						i = 0
						loop do
							tst = @submitted_tests[0]
							run = @test_runner.runs.find{|run| run.run_name == tst}
							if [:Complete, :Failed].include? run.status 
								run_checks run
								@submitted_tests.delete(tst)
							end
							i += 1
							break if i >= @submitted_tests.size
						end
						sleep 3
					end
			end

		
			# Run the test case with the given input file. 
			
			def self.submit_run(test_case)
				
				dir = results_folder(test_case)
				input_text = File.read("#{test_case_folder(test_case)}/#{test_case}.in")
				
				description = input_text.scan(Regexp.new("#{/^\s*!\s*description\s*=\s* /}(#{Regexp.quoted_string})")).flatten[0]
				(eputs "----Rejecting '#{file}', no description provided or description is not in the correct format: description = \" description \"";return) unless description
				run = Gs2.new(@test_runner)

				# This cryptic statement updates the run parameters from the input file 
				run.instance_eval(
					Gs2.defaults_file_text_from_input_file("#{test_case_folder(test_case)}/#{test_case}.in"))
				
				run.run_name = test_case
				run.instance_variable_set(:@dir_name, dir)
				@test_runner.test_submission = true
				@test_runner.submit(run)
				
				other_files = Dir.entries(test_case_folder(test_case)).find_all do |f|
					not (f =~ /\.in/ or f =~ /\.out\.nc/ or f =~ /\.svn/ or [".", ".."].include? f)
				end
				other_files.each do |f|
# 					p f
					FileUtils.cp(test_case_folder(test_case)+'/'+f, results_folder(test_case)+'/'+f)
				end
				
				@test_runner.test_submission = false
				@test_runner.submit(run)
				@submitted_tests.push run.run_name
				
				
			end
			
			def self.run_checks(run)
				checks_failed = []
				
				input_text = File.read("#{test_case_folder(run.run_name)}/#{run.run_name}.in")
				description = input_text.scan(Regexp.new("#{/^\s*!\s*description\s*=\s* /}(#{Regexp.quoted_string})")).flatten[0]
				description.gsub!(/(.{25,45} |.{45})/){"#$1\n"} if description
				if input_text =~ Regexp.new("#{/^\s*!\s*custom_checks\s*=\s* /}(#{Regexp.properly_nested("\\[", "\\]", false)})")
					custom_checks = eval($1)
					custom_checks.each{|check| checks_failed.push check unless run.run_check check}
				else
					checks_failed.push :standard unless run.run_check :standard
				end
				unless checks_failed.size == 0
					@tests_failed[run.run_name] = {checks_failed: checks_failed, description: description}
				end
			end
			
			# Relative error which differences cannot exceed in the standard case
			
			TOLERANCE = 1.0e-10
			
			# A list of variables that are allowed to be different in the standard test
			
			EXCLUDED_VARIABLES = ['input_file']
			
			# Check the results against standard cases using the checks described here. Any custom tests should be implemented here!
			
			def run_check(check)
				netcdf_standard_case = NumRu::NetCDF.open(TestGs2.test_case_folder(@run_name) + '/' + @run_name + '.out.nc')
				Dir.chdir(@directory) do
					case check
					when :standard
						netcdf = NumRu::NetCDF.open(@run_name + '.out.nc')
						netcdf.vars.map{|v| v.name}.each do |v|
							unless netcdf_standard_case.var(v)
								TEST_OUT.puts "Warning: variable '#{v}' is missing from the test case netcdf output for '#@run_name'. Suggest updating the test case netcdf file. This is not a GS2 fault."
							end
						end
						netcdf_standard_case.vars.map{|v| v.name}.each do |v|
							next if EXCLUDED_VARIABLES.include? v
							begin
								unless netcdf.var(v)
									TEST_OUT.puts "Error: Variable #{v} is missing from the netcdf output for #@run_name"
									return false
								end
								
								narray = netcdf.var(v).get
								standard_narray = netcdf_standard_case.var(v).get
								
								if standard_narray.abs.max < TOLERANCE 
									if narray.abs.max < TOLERANCE 
										next
									else
										TEST_OUT.puts "Error: Variable '#{v}' has failed check in '#@run_name'"
										return false
									end
								end
								
								
								difference = narray - standard_narray
								return false unless difference.abs.max / standard_narray.abs.max < TOLERANCE
	# 							ep 'var', difference.abs.max
							rescue => err
								TEST_OUT.puts "Error: #{err}"
								TEST_OUT.puts "Error: Variable '#{v}' has failed check in '#@run_name'"
								return false
							end
						end
						return true
					end
				end
			end
				                             
			
		end #module TestGs2
		
	end
end
				                  