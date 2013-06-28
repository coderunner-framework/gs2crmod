require 'helper'

CYCLONE_LOW_RES_FOLDER = 'test/cyclone_low_res'
class TestBasics < Test::Unit::TestCase
	def setup
    @runner = CodeRunner.fetch_runner(Y: 'test/slab_itg', C: 'gs2', X: '/dev/null')
  end
  def teardown
    FileUtils.rm('test/slab_itg/.code_runner_script_defaults.rb')
    FileUtils.rm('test/slab_itg/.CODE_RUNNER_TEMP_RUN_LIST_CACHE')
  end
  def test_basics
    assert_equal(@runner.run_class, CodeRunner::Gs2)
  end
end
if ENV['GS2_EXEC']
	class TestSubmission < Test::Unit::TestCase
		def setup
			CodeRunner.setup_run_class('gs2')
			CodeRunner::Gs2.make_new_defaults_file('test_gs2crmod', 'test/cyclone_low_res.in')
			FileUtils.mv('test_gs2crmod_defaults.rb', CodeRunner::Gs2.rcp.user_defaults_location + '/.')
			FileUtils.makedirs(tfolder)
		end
		def tfolder
			CYCLONE_LOW_RES_FOLDER
		end
		def teardown
			FileUtils.rm(CodeRunner::Gs2.rcp.user_defaults_location + '/' + 'test_gs2crmod_defaults.rb')
			FileUtils.rm(tfolder + '/.CODE_RUNNER_TEMP_RUN_LIST_CACHE')
			FileUtils.rm(tfolder + '/v/id_1/.code_runner_run_data')
			FileUtils.rm(tfolder + '/v/id_1/code_runner_results.rb')
			# Don't uncomment the line below unless you *really* know what you are doing!
			# Replacing the test archive will break many of the tests
			#Dir.chdir('test'){system "tar -czf cyclone_low_res.tgz cyclone_low_res/" unless FileTest.exist?('cyclone_low_res.tgz')}
			FileUtils.rm_r(tfolder)
		end
		def test_submission
			CodeRunner.submit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '4', Y: tfolder)
			CodeRunner.status(Y: tfolder)
		end
	end
else
	puts "\n************************************\nWarning: submission tests not run. Please specify the evironment variable GS2_EXEC (the path to the gs2 executable) if you wish to test submission.\n************************************\n"
	sleep 0.1
end
class TestAnalysis < Test::Unit::TestCase
	def setup
		Dir.chdir('test'){assert(system "tar -xzf cyclone_low_res.tgz")}
		@runner = CodeRunner.fetch_runner(Y: tfolder)
	end
	def test_analysis
		assert_equal(1, @runner.run_list.size)
		assert_equal(0.13066732664774272, @runner.run_list[1].max_growth_rate)
		assert_equal(0.13066732664774272, @runner.run_list[1].growth_rate_at_ky[0.5])
		assert_equal(:Complete, @runner.run_list[1].status)
	end
	def test_graphs
		kit = @runner.run_list[1].graphkit('phi2_by_ky_vs_time', {ky_index: 2})
		#kit.gnuplot
		assert_equal(51, kit.data[0].y.data.size)
		assert_equal(@runner.run_list[1].netcdf_file.var('phi2_by_ky').get('start' => [1,4], 'end' => [1,4]).to_a[0][0], kit.data[0].y.data[4])
		kit = @runner.run_list[1].graphkit('phi_real_space_surface', {n0: 3, Rgeo: 3, interpolate_theta: 4})
		kit.gnuplot
	end

	def tfolder
		CYCLONE_LOW_RES_FOLDER
	end
	def teardown
		FileUtils.rm_r(tfolder)
	end
end


