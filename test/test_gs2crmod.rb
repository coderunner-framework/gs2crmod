require 'helper'

CYCLONE_LOW_RES_FOLDER = 'test/cyclone_low_res'
class TestBasics < Test::Unit::TestCase
    def setup
    FileUtils.makedirs('test/slab_itg')
    @runner = CodeRunner.fetch_runner(Y: 'test/slab_itg', C: 'gs2', X: '/dev/null')
  end
  def teardown
    FileUtils.rm('test/slab_itg/.code_runner_script_defaults.rb')
    FileUtils.rm('test/slab_itg/.CODE_RUNNER_TEMP_RUN_LIST_CACHE')
    FileUtils.rmdir('test/slab_itg')
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
            # Don't uncomment the line below unless you *really* know what you are doing! Replacing the test archive will break many of the tests
            Dir.chdir('test'){system "tar -czf cyclone_low_res.tgz cyclone_low_res/" unless FileTest.exist?('cyclone_low_res.tgz')}
            FileUtils.rm_r(tfolder)
        end
        def test_submission
            CodeRunner.submit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '2', Y: tfolder, p: '{write_moments: ".true.", write_line: ".true.", save_for_restart: ".true.", nsave: 40}')
            CodeRunner.submit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '2', Y: tfolder, p: '{restart_id: 1, tprim_1: 4.0}')
            CodeRunner.submit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '2', Y: tfolder, p: '{restart_id: 2, nstep: 10}', gs2_options: {show_opt: true})
            CodeRunner.submit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '2', Y: tfolder, p: '{dump_response: ".true.", save_for_restart:".true.", nsave:40}')
            CodeRunner.submit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '2', Y: tfolder, p: '{restart_id: 4, read_response:".true.", nstep:10}')
            CodeRunner.resubmit(C: 'gs2', X: ENV['GS2_EXEC'], D: 'test_gs2crmod', n: '2', Y: tfolder, p: '{response_id: 4}', j:4)
            runs = CodeRunner.fetch_runner(Y: tfolder).run_list
            assert_equal(50, runs[1].completed_timesteps)
            assert_equal(4.0, runs[3].tprim_1)
            assert_equal(runs[1].gsl_vector('phi2tot_over_time')[-1].round(4), runs[2].gsl_vector('phi2tot_over_time')[0].round(4))
            CodeRunner.status(Y: tfolder)
            assert_equal(true, Dir.exists?("#{tfolder}/v/id_4/response"))
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
        @run = @runner.run_list[1]
    end
    def test_analysis
        assert_equal(1, @runner.run_list.size)
        assert_equal(0.13067, @runner.run_list[1].max_growth_rate.round(5))
        assert_equal(0.13067, @runner.run_list[1].growth_rate_at_ky[0.5].round(5))
        assert_equal(:Complete, @runner.run_list[1].status)
        p @run.frequency_at_ky_at_kx, @run.gsl_vector('kx')[1]
        assert_equal(5.8e-01, @run.frequency_at_ky_at_kx[0.5][2.5133].round(2))
        #p @run.gsl_vector('kx'); STDIN.gets
        assert_equal(5.8e-01, @run.frequency_at_ky_at_kx[0.5][@run.gsl_vector('kx')[3]].round(2))
    end
    def test_interpolation
        assert_equal(5, @run.gsl_vector('kx').size)
        assert_equal(17, @run.gsl_vector('kx', interpolate_x: 4).size)
        kxvec = @run.gsl_vector('kx', interpolate_x: 4)
        assert_equal(0.0, kxvec[(kxvec.size-1)/2])
        xvec4 = @run.gsl_vector('x', interpolate_x: 4)
        s = kxvec.size
        assert_equal((-(s-1)/2..(s-1)/2).to_a.map{|i| i.to_f * kxvec.to_box_order[1]}.to_gslv, kxvec)
        #p xvec4.to_a
        xvec = @run.gsl_vector('x')
        # The vectors don't contain the periodic point, but we can construct them 
        # and they should be equal
        assert_equal(xvec[-1] + xvec[-1] - xvec[-2], xvec4[-1] + xvec4[-1] - xvec4[-2])
        #p  @run.gsl_vector('kx', interpolate_x: 4).to_a
        kyvec4 = @run.gsl_vector('ky', interpolate_y: 4)
        #p kyvec4
        s = kyvec4.size
        assert_equal((0...s).to_a.map{|i| i.to_f * kyvec4[1]}.to_gslv, kyvec4)
        yvec4 = @run.gsl_vector('y', interpolate_y: 4)
        #p yvec4.to_a
        yvec = @run.gsl_vector('y')
        assert_equal(yvec[-1] + yvec[-1] - yvec[-2], yvec4[-1] + yvec4[-1] - yvec4[-2])
        #kit = @runner.run_list[1].graphkit('phi_real_space_poloidal_plane', {n0: 1, Rgeo: 3, interpolate_theta: 8, torphi: Math::PI/4.0, interpolate_x: 8})
        kit1 = @runner.run_list[1].graphkit('phi_real_space_surface', {n0: 1, Rgeo: 3, gs2_coordinate_factor: 1.0})
        shape = kit1.data[0].x.data.shape
        shape = shape.map{|s| (s-1)*2 +1}
        kit = @runner.run_list[1].graphkit('phi_real_space_surface', {n0: 1, Rgeo: 3, interpolate_theta: 2, interpolate_x: 2, interpolate_y: 2, gs2_coordinate_factor: 1.0})
        assert_equal(shape, kit.data[0].x.data.shape)
        assert_equal(shape, kit.data[0].f.data.shape)
        #kit.gp.view = ["equal xyz", ",,5.0"]
        #kit.gnuplot 
    end
    def test_graphs
        kit = @runner.run_list[1].graphkit('phi2_by_ky_vs_time', {ky_index: 2})
        assert_equal(51, kit.data[0].y.data.size)
        assert_equal(@runner.run_list[1].netcdf_file.var('phi2_by_ky').get('start' => [1,4], 'end' => [1,4]).to_a[0][0], kit.data[0].y.data[4])
        
        kit = @run.graphkit('tpar2_by_mode_vs_time', {ky_index:2, kx_index:1, species_index:1})
        assert_equal(@runner.run_list[1].netcdf_file.var('tpar2_by_mode').get('start' => [0,1,0,4], 'end' => [0,1,0,4]).to_a[0][0][0][0], kit.data[0].y.data[4])

        kit = @run.graphkit('tperp2_by_mode_vs_time', {ky_index:2, kx_index:1, species_index:1})
        assert_equal(@runner.run_list[1].netcdf_file.var('tperp2_by_mode').get('start' => [0,1,0,4], 'end' => [0,1,0,4]).to_a[0][0][0][0], kit.data[0].y.data[4])

    #Test heat flux as a function of kxy
        kit = @run.graphkit('es_heat_flux_vs_kx', {species_index:1})
        assert_equal(5, kit.data[0].y.data.size)
        kit = @run.graphkit('es_heat_flux_vs_ky', {species_index:1})
        assert_equal(3, kit.data[0].y.data.size)

    #Test zonal flow velocity calculation
        kit = @run.graphkit('zf_velocity_vs_x', {theta_index:4})
        assert_equal(5, kit.data[0].y.data.size)

    #Time averaged kx and ky spectra
        kit = @runner.run_list[1].graphkit('kx_spectrum_avg')
        assert_equal(5, kit.data[0].y.data.size)
        kit = @runner.run_list[1].graphkit('ky_spectrum_avg')
        assert_equal(3, kit.data[0].y.data.size)
    end
    def test_3d_graphs
        kit = @runner.run_list[1].graphkit('phi_real_space', {n0: 3, Rgeo: 3})
        assert_equal([5,5,9], kit.data[0].f.data.shape)
        assert_equal(-0.00492, kit.data[0].f.data[2,3,6].round(5))
        kit = @runner.run_list[1].graphkit('density_real_space', {n0: 3, Rgeo: 3, species_index: 1, gs2_coordinate_factor: 0.9})
        #kit.gnuplot
        assert_equal([5,5,9], kit.data[0].f.data.shape)
        assert_equal(-0.00985, kit.data[0].f.data[2,3,6].round(5))
        kit = @runner.run_list[1].graphkit('phi_real_space', {n0: 3, Rgeo: 3, t_index: 2})
        assert_equal([5,5,9], kit.data[0].f.data.shape)
        assert_equal(0.00015, kit.data[0].f.data[2,3,6].round(5))
        kit = @runner.run_list[1].graphkit('phi_real_space_surface', {n0: 3, Rgeo: 3, interpolate_theta: 2})
        assert_equal([5,5,1], kit.data[0].f.data.shape)
        assert_equal([5,1,17], kit.data[2].f.data.shape)
        assert_equal(-0.0024, kit.data[0].f.data[1,4,0].round(5))
        #kit.gnuplot
        kit = @runner.run_list[1].graphkit('phi_real_space_poloidal_plane', {n0: 1, Rgeo: 3, interpolate_theta: 8, torphi: Math::PI/4.0})
        assert_equal(-0.00176, kit.data[0].f.data[-1,1].round(5))
        assert_equal(1.707, kit.data[0].x.data[-1,1].round(3))
        kit.gp.view = ["equal xyz", ",,4.0"]
        #kit.gnuplot
        kit = @runner.run_list[1].graphkit('density_real_space_poloidal_plane', {n0: 1, Rgeo: 3, interpolate_theta: 8, torphi: Math::PI/4.0, species_index: 1})
        assert_equal(-0.00352, kit.data[0].f.data[-1,1].round(5))
        assert_equal(1.707, kit.data[0].x.data[-1,1].round(3))
        kit.gp.view = ["equal xyz", ",,4.0"]

        kit = @runner.run_list[1].graphkit('density_real_space_poloidal_plane', {n0: 1, Rgeo: 3, interpolate_theta: 8, torphi: Math::PI/4.0, species_index: 1, t_index: 50})
        assert_equal(-0.00208, kit.data[0].f.data[-1,1].round(5))
        assert_equal(1.707, kit.data[0].x.data[-1,1].round(3))
        kit.gp.view = ["equal xyz", ",,4.0"]
        #kit.gnuplot
        kit = @runner.run_list[1].graphkit('phi_real_space_standard_representation', {n0: 1, Rgeo: 3, interpolate_theta: 2, torphi_values: [Math::PI/4.0,3.0*Math::PI/4.0], interpolate_y: 2})
        assert_equal([5,17], kit.data[0].f.data.shape)
        assert_equal([3,1,17], kit.data[2].f.data.shape)
        assert_equal([3,1,17], kit.data[2].y.data.shape)
        assert_equal(0.00012, kit.data[0].f.data[-1,1].round(5))
        assert_equal(2.12132, kit.data[2].y.data[2,0,12].round(5))
        #kit.gnuplot 
    end

   #This test requires an update of the test NetCDF file to the new GS2 format!  
   #def test_new_netcdf_module
   #    assert_equal(NumRu::NetCDF, @run.new_netcdf_file.class)
   #    assert_equal(NumRu::NetCDFVar, @run.new_netcdf_file.var('phi2').class)
   #    assert_equal(@run.netcdf_file.var('phi2').get.to_a[-1], @run.new_netcdf_file.var('phi2').get.to_a[-1])
   #    assert_equal("r",  @run.netcdf_smart_reader.dimensions('phi')[0].name)
   #    assert_equal("X",  @run.netcdf_smart_reader.dimensions('phi')[2].name)
   #    assert_equal(2, @run.netcdf_smart_reader.dim_start("Y", Y_index: 3))
   #    assert_equal(@run.new_netcdf_file.var('phi2').get.to_a, @run.netcdf_smart_reader.read_variable('phi2', {}).to_a)
   #    assert_equal(@run.new_netcdf_file.var('phi2').get.to_a[3], @run.netcdf_smart_reader.read_variable('phi2', {t_index: 4})[0])
   #    assert_equal(@run.new_netcdf_file.var('phi2').get.to_a[3], @run.netcdf_smart_reader.read_variable('phi2', {tmax: 3})[-1])
   #    assert_equal(@run.new_netcdf_file.var('phi').get[0,4,3,1], @run.netcdf_smart_reader.read_variable('phi', {zmax: 5, X_index: 4, Y_element: 1})[0,-2])
   #    assert_equal("time (a/v_thr)", @run.smart_graphkit(graphkit_name: 'cdf_heat_flux_tot').xlabel)
   #    @run.smart_graphkit(graphkit_name: 'cdf_phi', r_index: 1, z_index: 5)
   #    @run.old_smart_graphkit(graphkit_name: 'nc_phi', ri_index: 1, theta_index: 5)
   #    @run.graphkit('nc_phi', ri_index: 1, theta_index: [5,6])
   #    #@runner.run_class.help_graphs
   #end

    def tfolder
        CYCLONE_LOW_RES_FOLDER
    end
    def teardown
        FileUtils.rm_rf(tfolder)
    end
end

AGK_SLAB_ITG_LOW_KPERP_FOLDER = 'test/agk_slab_itg_low_kperp'

class TestAgkAnalysis < Test::Unit::TestCase

    def setup
        Dir.chdir('test'){assert(system "tar -xzf agk_slab_itg_low_kperp.tgz")}
        @runner = CodeRunner.fetch_runner(Y: tfolder)
        @run = @runner.run_list[1]
    end
    def tfolder
        AGK_SLAB_ITG_LOW_KPERP_FOLDER
    end
    def test_graphs
        kit = @runner.run_list[1].graphkit('phi2_by_ky_vs_time', {ky_index: 2})
        #kit.gnuplot
        assert_equal(126, kit.data[0].y.data.size)
        assert_equal(@runner.run_list[1].netcdf_file.var('phi2').get('start' => [4], 'end' => [4]).to_a[0], kit.data[0].y.data[4])
    assert_equal(nil,@runner.run_list[1].instance_variable_get(:@shat))
        kit = @runner.run_list[1].graphkit('kpar_spectrum',{})
        #kit.gnuplot
    end
    def test_analysis
        assert_equal(2, @runner.run_list.size)
        # These tests were commented out when growth rate and frequency 
        # calculations moved to reading directly from the netcdf file.
        # Instead of comparing with a calculated result, compare with value
        # read directly from the NetCDF file.
        #assert_equal(0.04154, @runner.run_list[1].max_growth_rate.round(5))
        #p @runner.run_list[1].growth_rate_at_ky
        #assert_equal(0.04154, @runner.run_list[1].growth_rate_at_ky[0.01].round(5))
        assert_equal(:Complete, @runner.run_list[1].status)
    end
    def test_3d_graphs
        kit = @runner.run_list[2].graphkit('phi_real_space_surface', {rho_star: 0.1})
        #kit.gnuplot
        assert_equal([5,5,1], kit.data[0].f.data.shape)
        #assert_equal(-0.00402, kit.data[0].f.data[2,3,6].round(5))
    end
    def teardown
        FileUtils.rm_rf(tfolder)
    end
end

if ENV['AGK_EXEC']
    class TestAstrogkSubmission < Test::Unit::TestCase
        def setup
            CodeRunner.setup_run_class('gs2', modlet: 'astrogk')
            CodeRunner::Gs2::Astrogk.make_new_defaults_file('test_gs2crmod_astrogk', 'test/agk_slab_itg_low_kperp.in')
            FileUtils.mv('test_gs2crmod_astrogk_defaults.rb', CodeRunner::Gs2::Astrogk.rcp.user_defaults_location + '/.')
            FileUtils.makedirs(tfolder)
        end
        def tfolder
            AGK_SLAB_ITG_LOW_KPERP_FOLDER
        end
        def teardown
            FileUtils.rm(CodeRunner::Gs2::Astrogk.rcp.user_defaults_location + '/' + 'test_gs2crmod_astrogk_defaults.rb')
            FileUtils.rm(tfolder + '/.CODE_RUNNER_TEMP_RUN_LIST_CACHE')
            FileUtils.rm(tfolder + '/v/id_1/.code_runner_run_data')
            FileUtils.rm(tfolder + '/v/id_1/code_runner_results.rb')
             #Don't uncomment the line below unless you *really* know what you are doing! Replacing the test archive will break many of the tests
            Dir.chdir('test'){system "tar -czf agk_slab_itg_low_kperp.tgz agk_slab_itg_low_kperp/" unless FileTest.exist?('agk_slab_itg_low_kperp.tgz')}
            FileUtils.rm_r(tfolder)
        end
        def test_submission
            CodeRunner.submit(C: 'gs2', X: ENV['AGK_EXEC'], D: 'test_gs2crmod_astrogk', n: '4', Y: tfolder, m: 'astrogk')
            CodeRunner.submit(C: 'gs2', X: ENV['AGK_EXEC'], D: 'test_gs2crmod_astrogk', n: '4', Y: tfolder, m: 'astrogk', p: '{ny: 8, nx: 8, y0: 100, x0: 100, grid_option: "box"}', )
            
            CodeRunner.status(Y: tfolder)
        end
    end
else
    puts "\n************************************\nWarning: agk submission tests not run. Please specify the evironment variable AGK_EXEC (the path to the agk executable) if you wish to test submission.\n************************************\n"
    sleep 0.1
end

class TestBasicsSpectroGK < Test::Unit::TestCase
    def setup
    @runner = CodeRunner.fetch_runner(Y: 'test/slab_itg', C: 'gs2', X: '/dev/null', m: 'spectrogk')
  end
  def teardown
    FileUtils.rm('test/slab_itg/.code_runner_script_defaults.rb')
    #FileUtils.rm('test/slab_itg/.CODE_RUNNER_TEMP_RUN_LIST_CACHE')
  end
  def test_basics
    assert_equal(@runner.run_class, CodeRunner::Gs2::Spectrogk)
  end
    def test_variables
        assert_equal(Hash, @runner.run_class.rcp.namelists.class)
        assert_equal(@runner.run_class.rcp.namelists[:layouts_knobs].class, Hash)
        assert_equal(@runner.run_class.rcp.namelists[:parameters_knobs].class, Hash)
        assert_equal(@runner.run_class.rcp.namelists[:parameters_knobs][:variables][:force_5d].class, Hash)
    end
end
