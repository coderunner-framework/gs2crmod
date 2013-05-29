{:parameters=>
  {:description=>"GENERAL PARAMETERS",
   :should_include=>"true",
   :variables=>
    {:beta=>
      {:help=>
        "In GS2 <math>\\beta</math> is the ratio of the reference pressure to the reference magnetic energy density, <math>\\beta = 2 \\mu_0 n_{\\rm ref} T_{\\rm ref}/B_{\\rm ref}^2 </math>. Mainly GS2 uses <math>\\beta</math> to determine the amplitude of the perturbed magnetic fields.   \n**For electromagnetic runs, the contribution of each species to the total parallel current is weighted by a factor of <math>w_s = 2 \\beta Z_s n_s \\sqrt{T_s/m_s}</math>.\n**For electromagnetic runs that include <math>\\delta B_\\parallel</math>, this field is proportional to <math>\\beta</math>.\n**The contribution of <math>(\\delta B)^2</math> to the total gyrokinetic energy is inversely proportional to this input parameter.\n**If an antenna is set up to drive Alfven waves, then <math>\\beta</math> is used to calculate the Alfven frequency.  \n**For some collision operator models, <math>\\beta</math> is used to calculate the resistivity.  \n**For some rarely-used initial conditions, <math>\\beta</math> appears explicitly. \n**Important:  <math>\\beta</math> is not a GS2 plasma equilibrium parameter, but it must be chosen with care. <math>\\beta</math> is '''not''' automatically consistent with the value of <math>\\beta'</math> used in the local magnetic equilibrium.  The user is responsible for ensuring that the values of <math> \\beta </math> together with the densities, temperatures and gradients for all species are consistent with the local equilibrium value of <math>\\beta'</math>.",
       :should_include=>"true",
       :description=>
        "Ratio of particle to magnetic pressure (reference Beta, not total beta):  beta=n_0 T_0 /( B^2 / (8 pi))",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_eik},
     :zeff=>
      {:help=>
        "Effective ionic charge.  The parameter <math>Z_{\\rm eff}</math> appears only in the electron collision frequency, and is not automatically set to be consistent with the mix of species specified in the species namelists.",
       :should_include=>"true",
       :description=>"Effective ionic charge.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:run_parameters},
     :teti=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        "Ratio of electron to ion temperatures.  Deprecated.  Do not use unless you know what you are doing.",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:teti,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-100.0],
       :type=>:Float,
       :code_name=>:teti,
       :module=>:run_parameters},
     :parameters_k0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:k0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:k0},
     :tite=>
      {:help=>
        "Ratio of ion to electron temperatures.  This parameter is used only when there is no species called \"electron\" included.",
       :should_include=>"true",
       :description=>"Ratio of ion to electron temperatures.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:run_parameters},
     :rhostar=>
      {:should_include=>"true",
       :description=>"Rhostar, for low flow terms.",
       :help=>"Rhostar, for low flow terms.",
       :code_name=>:rhostar,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float}}},
 :kt_grids_knobs=>
  {:description=>nil,
   :help=>
    "There are two variables which may be determined at the highest level\nof the kt_grids module: grid_option and norm_option.  The first determines which\nFourier modes perpendicular to the field lines will be evolved in\ntime.  For each possible choice, there is a subsidiary module (and\nnamelist) which controls the details.  The norm_option variable\nsimply fixes the definition of the thermal velocity.",
   :should_include=>"true",
   :variables=>
    {:grid_option=>
      {:help=>
        "The general layout of the perpendicular grid.\n** 'single' Evolve a single Fourier mode.  Set up kt_grids_single_parameters.\n** 'default' Same as 'single'\n** 'range' Evolve a range of equally spaced Fourier modes. Set up kt_grids_range_parameters.\n** 'specified' Evolve an arbitrary set of Fourier modes.  Set up kt_grids_specified_parameters.\n** 'box' Simulation domain is logically rectangular.  Set up kt_grids_box_parameters.\n** 'nonlinear' Same as 'box.'\n** 'xbox' Experimental.",
       :should_include=>"true",
       :description=>"The general layout of the perpendicular grid.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "single",
         "range",
         "specified",
         "box",
         "nonlinear",
         "xbox"],
       :module=>:kt_grids}}},
 :kt_grids_box_parameters=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:nx=>
      {:help=>
        "The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.",
       :should_include=>"true",
       :description=>
        "The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[0],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:kt_grids_specified},
     :ny=>
      {:help=>
        "The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.",
       :should_include=>"true",
       :description=>
        "The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[0],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:kt_grids_specified},
     :jtwist=>
      {:help=>
        "For finite magnetic shear, determines box size in x direction according to \n** <math>L_x = L_y  jtwist / (2 \\pi \\hat{s}) </math>\n** Also affects the number of \"connections\" at each ky when linked boundary conditions are selected in the dist_fn_knobs namelist.",
       :should_include=>"true",
       :description=>"L_x = L_y  jtwist / (2 pi shat)",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[1],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:kt_grids_box},
     :y0=>
      {:help=>
        "The length of the box in the y direction (measured in the Larmour radius of species 1).  Box size in y direction is 2*pi*y0.",
       :should_include=>"true",
       :description=>
        "The length of the box in the y direction (measured in the Larmour radius of species 1)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[2.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:kt_grids_box},
     :x0=>
      {:help=>
        "The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)",
       :should_include=>"true",
       :description=>
        "The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:kt_grids_box},
     :naky=>
      {:help=>
        "The actual number of ky modes (do not use for nonlinear runs, use ny). If left as zero (recommended), automatically set to (ny-1)/3+1.",
       :should_include=>"true",
       :description=>
        "The actual number of ky modes (do not use for nonlinear runs, use ny)",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[0, 1],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:kt_grids_specified},
     :ntheta0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " If left as zero (recommended), automatically set to 2*((nx-1)/3)+1.\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:ntheta0,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0, 1],
       :type=>:Integer,
       :code_name=>:ntheta0,
       :module=>:kt_grids_specified},
     :ly=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Box size in y direction.\n** If ly=0, it is set to be 2*pi*y0 (below).\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:ly,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:ly,
       :module=>:hyper},
     :rtwist=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Experts only.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:rtwist,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:rtwist,
       :module=>:kt_grids_box},
     :nkpolar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nkpolar,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0],
       :type=>:Integer,
       :code_name=>:nkpolar,
       :module=>:kt_grids_box}}},
 :kt_grids_single_parameters=>
  {:description=>"",
   :should_include=>"@grid_option=='single'",
   :variables=>
    {:aky=>
      {:help=>
        "<math>k_y \\rho</math> for the reference species. Used only in certain modes (e.g. in single mode).",
       :should_include=>"true",
       :description=>"The actual value of ky rho",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.2, 0.4],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:kt_grids_single},
     :theta0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" theta_0 \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:theta0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:theta0,
       :module=>:kt_grids_single},
     :akx=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" kx rho (but set theta_0 instead.)\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:akx,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:akx,
       :module=>:kt_grids_single}}},
 :theta_grid_parameters=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:ntheta=>
      {:help=>
        "Number of grid points along equilibrium magnetic field between <math>\\theta=(-\\pi,\\pi)</math> (in addition to a grid point at <math>\\theta=0</math>).\n** Ignored in some cases",
       :should_include=>"true",
       :description=>
        "Number of points along field line (theta) per 2 pi segment",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[24, 32, 64],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:kt_grids_specified},
     :nperiod=>
      {:help=>
        "Sets the number of <math>2\\pi</math> segments along equilibrium magnetic field to include in simulation domain.  <math>N_{\\rm segments} = (2n_{\\rm period} - 1)</math>.\n** Ignored in some cases",
       :should_include=>"true",
       :description=>
        "Number of 2 pi segments along equilibrium magnetic field.",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[1, 2],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:theta_grid_params},
     :eps=>
      {:help=>
        "Controls particle trapping (among other things) in simple geometric models.  <math>\\epsilon = r/R</math>",
       :should_include=>"true",
       :description=>"eps=r/R",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.3, 0.5],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_gridgen},
     :epsl=>
      {:help=>
        "Sets curvature drift in simple geometric models.  <math>\\epsilon_\\ell = 2 a / R </math>, where a is the GS2 normalisation length and R is the major radius at the centre of the flux surface.",
       :should_include=>"true",
       :description=>"epsl=2 a/R",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.3],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :shat=>
      {:help=>
        "Sets value of magnetic shear in simple geometric models.\n** over-ridden by s_hat_input in theta_grid_eik_knobs for most values of bishop.",
       :should_include=>"true",
       :description=>"",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0, 0.75],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :pk=>
      {:help=>
        "<math>p_k = 2 a / q R</math> - when using high aspect ratio, <math>s-\\alpha</math> model for geometry it sets q, the magnetic safety factor, since the ratio 2a/R is set by epsl. It therefore also sets the connection length, i.e. the length of the box in the parallel direction, in terms of a, the GS2 normalization quantity.",
       :should_include=>"true",
       :description=>"pk = 2 a / q R",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.3],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_gridgen},
     :kp=>
      {:help=>
        "If kp > 0 then pk = 2*kp. Set kp rather than pk in the slab. Since in the slab limit, shat =  2a/ (L_S pk), this means that if kp = 1, shat = a/ L_S, where L_S is the magnetic shear scale length.\n\nSets parallel box length when geometry is unsheared slab.  <math>k_p = 2 \\pi /L_z</math>.",
       :should_include=>"true",
       :description=>"if kp >0 then pk = 2*kp",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[-1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:hyper},
     :rhoc=>
      {:help=>
        "rhoc is flux surface label (0< rhoc< 1). Its exact meaning depends on irho. Usually rho = diameter/diameter of LCFS\n** When irho = 1, rhoc = sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)\n** When irho = 2, rhoc =  diameter/(diameter of LCFS).  recommended\n** When irho = 3, rhoc =  poloidal flux/(poloidal flux of LCFS)",
       :should_include=>"true",
       :description=>
        "Flux surface label. Usually rho = diameter/diameter of LCFS",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.5, 0.8, 1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :qinp=>
      {:help=>
        "Sets value of the safety factor when using local toroidal equilibrium model.",
       :should_include=>"true",
       :description=>
        "Sets value of the safety factor when using local toroidal equilibrium model.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.5],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :akappa=>
      {:help=>
        "Sets local elongation when local toroidal equilibrium is specified.",
       :should_include=>"true",
       :description=>
        "Sets local elongation when local toroidal equilibrium is specified.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :akappri=>
      {:help=>
        "Sets local gradient of elongation when local toroidal equilibrium is specified.\n** akappri = <math> d\\kappa/d\\rho </math>",
       :should_include=>"true",
       :description=>"akappri = dkappa/drho",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :tri=>
      {:help=>
        "Sets local triangularity when local toroidal equilibrium is specified.\n** triangularity is tri = arcsin[(R(max(Z))-R_major)/r_mid]",
       :should_include=>"true",
       :description=>"tri = arcsin[(R(max(Z))-R_major)/r_mid]",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :tripri=>
      {:help=>
        "Sets local gradient of triangularity when local toroidal equilibrium is specified.\n** tripri =  <math>dtri/d\\rho</math>",
       :should_include=>"true",
       :description=>"tripri =  dtri/drho",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :shift=>
      {:help=>
        "Related to Shafranov shift.  Sign depends on geometric model. ?? Exact definition uncertain... please edit if you know! ?? Could be\n** <math>shift = -R q**2 d\\beta/d\\rho (>0) </math> ??\n** dR/drho (R normalized to a) (< 0) ??",
       :should_include=>"true",
       :description=>"shift = -R q**2 dbeta/drho (>0)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0, 0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :alpmhd=>
      {:should_include=>"true",
       :description=>"Do not use unless you know what you are doing.",
       :help=>"Default = 0.  Do not use unless you know what you are doing.",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:alpmhd,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[],
       :type=>:Float,
       :code_name=>:alpmhd,
       :module=>:theta_grid_salpha},
     :asym=>
      {:should_include=>"true",
       :description=>nil,
       :help=>"*asympri:",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:asym,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:asym,
       :module=>:theta_grid_params},
     :asympri=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:asympri,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:asympri,
       :module=>:theta_grid_params},
     :rmaj=>
      {:help=>"Major radius/a (Position of magnetic axis)",
       :should_include=>"true",
       :description=>"Major radius/a (Position of magnetic axis)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :r_geo=>
      {:help=>"Major radius/a (centerpoint of LCFS)",
       :should_include=>"true",
       :description=>"Major radius/a (centerpoint of LCFS)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[3],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_params},
     :btor_slab=>
      {:should_include=>"true",
       :description=>
        "Determines the angle between the field and the flow in the slab limit",
       :help=>
        "In the slab limit, determines the angle between the field and the background flow (which flows in an imaginary toroidal direction). It is effectively equal to <math>\\frac{B_t}{B}=\\frac{u_{\\parallel}}{u}</math>.",
       :tests=>["Tst::FLOAT"],
       :code_name=>:btor_slab,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn}}},
 :theta_grid_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:equilibrium_option=>
      {:help=>
        "The equilibrium_option variable controls which geometric assumptions are used in the run.  Additional namelists determine the geometric parameters according to the choice made here. Allowed values are:\n** 'eik' Use routines from the geometry module, controlled mainly by the subsidiary namelists theta_grid_parameters and theta_grid_eik_knob.\n** 'default' Same as 'eik'                                                                                                                          \n** 's-alpha' Use high aspect-ratio toroidal equilbrium (which can be simplified to slab or cylinder), controlled by the subsidiary namelist theta_grid_parameters and  theta_grid_salpha_knob\n** 'file'  Use output from rungridgen code directly.  Controlled by theta_grid_file_knobs. Note: in this case, the variables nperiod and ntheta (from the theta_grid_parameters namelist) are ignored.\n** 'grid.out'  Same as 'file'",
       :should_include=>"true",
       :description=>
        "Controls which geometric assumptions are used in the run.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>["default", "eik", "s-alpha", "grid.out", "file"],
       :module=>:theta_grid},
     :gb_to_cv=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:gb_to_cv,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:gb_to_cv,
       :module=>:theta_grid}}},
 :theta_grid_salpha_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:model_option=>
      {:help=>
        " \n** 's-alpha' High aspect ratio toroidal equilibrium.  (Note that the curvature and grad-B drifts are equal.)\n** 'default' Same as 's-alpha'\n**\n'alpha1','rogers','b2','normal_only',const-curv',no-curvature': See output of ingen until this page is improved.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "s-alpha",
         "alpha1",
         "rogers",
         "b2",
         "normal_only",
         "const-curv",
         "no-curvature"],
       :module=>:theta_grid_salpha},
     :alpmhdfac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:alpmhdfac,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:alpmhdfac,
       :module=>:theta_grid_salpha},
     :alpha1=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:alpha1,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:alpha1,
       :module=>:theta_grid_salpha}}},
 :theta_grid_file_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:gridout_file=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Name of file with output from rungridgen.\n",
       :tests=>["Tst::STRING"],
       :gs2_name=>:gridout_file,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :autoscanned_defaults=>["grid.out"],
       :type=>:String,
       :code_name=>:gridout_file,
       :module=>:theta_grid_file},
     :no_geo_info=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:no_geo_info,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[],
       :type=>:Fortran_Bool,
       :code_name=>:no_geo_info,
       :module=>:theta_grid_file}}},
 :theta_grid_gridgen_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:npadd=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::INT"],
       :gs2_name=>:npadd,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[2],
       :type=>:Integer,
       :code_name=>:npadd,
       :module=>:theta_grid_gridgen},
     :alknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:alknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 0.1],
       :type=>:Float,
       :code_name=>:alknob,
       :module=>:theta_grid_gridgen},
     :epsknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:epsknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[],
       :type=>:Float,
       :code_name=>:epsknob,
       :module=>:theta_grid_gridgen},
     :bpknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:bpknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0e-08],
       :type=>:Float,
       :code_name=>:bpknob,
       :module=>:theta_grid_gridgen},
     :extrknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:extrknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 10.0],
       :type=>:Float,
       :code_name=>:extrknob,
       :module=>:theta_grid_gridgen},
     :tension=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tension,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:tension,
       :module=>:theta_grid_gridgen},
     :thetamax=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:thetamax,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:thetamax,
       :module=>:theta_grid_gridgen},
     :deltaw=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:deltaw,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:deltaw,
       :module=>:theta_grid_gridgen},
     :widthw=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:widthw,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:widthw,
       :module=>:theta_grid_gridgen}}},
 :theta_grid_eik_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:itor=>
      {:help=>" Do not change.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[1],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:theta_grid_eik},
     :iflux=>
      {:help=>
        " Choose mode of operation: \n** 0: Use Miller parameterization of local toroidal MHD equilibrium.\n** 1: Read equilibrium data from output of MHD code\n** 2: Running inside a transport code without numerical equilibrium\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[0, 1],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:theta_grid_eik},
     :irho=>
      {:help=>
        "Choose definition of flux surface coordinate\n** 1: rho == sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)\n** 2: rho == midplane diameter/LCFS diameter     Recommended\n** 3: rho == poloidal flux/poloidal flux of LCFS",
       :should_include=>"true",
       :description=>"Chooses definition of flux surface coordinate.",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[1, 2],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:theta_grid_eik},
     :local_eq=>
      {:help=>" Use Miller-style local equilibrium\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".false.", ".true."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:theta_grid_eik},
     :eqfile=>
      {:help=>" Name of file with numerical equilibrium data (if required)\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["dskeq.cdf", "test8_efit.dat"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:theta_grid_eik},
     :bishop=>
      {:help=>
        " Use Bishop relations to generate metric coefficients.\n** 0: Use high-aspect ratio coefficients (only for debugging)\n** 1: Use actual equilibrium values of shat, p'  Recommended \n** 2: Use numerical equilibrium + s_hat_input and p_prime_input\n** 3: Use numerical equilibrium + s_hat_input and inv_Lp_input\n** 4: Use numerical equilibrium + s_hat_input and beta_prime_input \n** 5: Use numerical equilibrium + s_hat_input and alpha_input\n** 6: Use numerical equilibrium + beta_prime_input\n** 7: Use numerical equilibrium and multiply pressure gradient by dp_mult\n** 8: Use numerical equilibrium + s_hat_input and multiply pressure gradient by dp_mult\n** 9: Use numerical equilibrium + s_hat_input and beta_prime_input\n** Otherwise: Use magnetic shear and pressure gradient as set elsewhere.\n** \n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[1, 5, 6],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:theta_grid_eik},
     :s_hat_input=>
      {:help=>" Overrides shat.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0, 0.8, 1],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_eik},
     :beta_prime_input=>
      {:help=>
        "The gradient of the pressure. Strictly speaking this parameter is not <math>\\frac{\\partial \\beta}{\\partial \\rho} </math> but <math>\\beta \\frac{1}{p}\\frac{\\partial p}{\\partial \\rho}</math>: in other words, the gradient of the magnetic field is ignored. Used only if bishop = 4 or 9.",
       :should_include=>"true",
       :description=>"The gradient of the pressure.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[-1, -0.0, 0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_eik},
     :delrho=>
      {:help=>" Default usually okay\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.01],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_eik},
     :isym=>
      {:help=>" \n** 0:  Recommended \n** 1: Force up-down symmetry.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[0],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:theta_grid_eik},
     :ppl_eq=>
      {:help=>" Use Menard-style NetCDF equilibrium (JSOLVER)\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false.", ".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:theta_grid_eik},
     :gen_eq=>
      {:help=>" Use Toq-style NetCDF equilibrium (TOQ)\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:theta_grid_eik},
     :efit_eq=>
      {:help=>" Use EFIT equilibrium (EFIT, codes with eqdsk format)\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false.", ".t."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:theta_grid_eik},
     :equal_arc=>
      {:help=>" Change field-line coordinate.  Recommended value: F\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false.", ".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:theta_grid_eik},
     :writelots=>
      {:help=>" Write a little extra about geometry to the screen.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false.", ".t."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:theta_grid_eik},
     :dfit_eq=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Vacuum magnetic dipole geometry\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:dfit_eq,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:dfit_eq,
       :module=>:theta_grid_eik},
     :idfit_eq=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:idfit_eq,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:idfit_eq,
       :module=>:theta_grid_eik},
     :gs2d_eq=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Read Colin Roach's GS2D equilibrium file\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:gs2d_eq,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:gs2d_eq,
       :module=>:theta_grid_eik},
     :transp_eq=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Use new PPL NetCDF equilibrium (TRANSP)\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:transp_eq,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false.", "F"],
       :type=>:Fortran_Bool,
       :code_name=>:transp_eq,
       :module=>:theta_grid_eik},
     :alpha_input=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Used only if bishop = 5.  Needs to be checked.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:alpha_input,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>["alpmhd_in"],
       :type=>:Float,
       :code_name=>:alpha_input,
       :module=>:theta_grid_eik},
     :dp_mult=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Used only if bishop = 7 or 8.  \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:dp_mult,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:dp_mult,
       :module=>:theta_grid_eik},
     :rmin=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Never used\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:rmin,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.01, 0.05],
       :type=>:Float,
       :code_name=>:rmin,
       :module=>:theta_grid_eik},
     :rmax=>
      {:should_include=>"true",
       :description=>nil,
       :help=>"  Never used\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:rmax,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 1.0],
       :type=>:Float,
       :code_name=>:rmax,
       :module=>:theta_grid_eik},
     :invlp_input=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Used only if bishop = 3.  \n",
       :tests=>["Tst::INT"],
       :gs2_name=>:invLp_input,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[5],
       :type=>:Integer,
       :code_name=>:invLp_input}}},
 :le_grids_knobs=>
  {:description=>"PITCH ANGLE/ENERGY GRID SETUP",
   :should_include=>"true",
   :variables=>
    {:vcut=>
      {:help=>
        "No. of standard deviations from the standard Maxwellian beyond which the distribution function will be set to 0",
       :should_include=>"true",
       :description=>
        "No. of standard deviations from the standard Maxwellian beyond which the distribution function will be set to 0",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[4.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:le_grids},
     :ngauss=>
      {:help=>
        "Number of untrapped pitch-angles moving in one direction along field line.",
       :should_include=>"true",
       :description=>
        "Number of untrapped pitch-angles moving in one direction along field line.",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[5],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :allowed_values=>[1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48],
       :module=>:le_grids},
     :negrid=>
      {:help=>"Total number of energy grid points",
       :should_include=>"true",
       :description=>"Total number of energy grid points",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[-10],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:le_grids},
     :bouncefuzz=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Should not have to change this\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:bouncefuzz,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[],
       :type=>:Float,
       :code_name=>:bouncefuzz,
       :module=>:le_grids},
     :nesuper=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Only used if advanced_egrid = F; sets number of energy grid points above ecut\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:nesuper,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[2],
       :type=>:Integer,
       :allowed_values=>[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15],
       :code_name=>:nesuper,
       :module=>:le_grids},
     :nesub=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Only used if advanced_egrid = F; sets number of energy grid points below ecut\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:nesub,
       :allowed_values=>
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 32, 48, 64],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[8],
       :type=>:Integer,
       :code_name=>:nesub,
       :module=>:le_grids},
     :le_grids_test=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Write some data to the screen and stop before doing much calculation. (Debugging only)\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:test,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:test},
     :trapped_particles=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:trapped_particles,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[],
       :type=>:Fortran_Bool,
       :code_name=>:trapped_particles,
       :module=>:le_grids},
     :nmax=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nmax,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>["size(dv)", "size(h)"],
       :type=>:Integer,
       :code_name=>:nmax,
       :module=>:le_grids},
     :wgt_fac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:wgt_fac,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[],
       :type=>:Float,
       :code_name=>:wgt_fac,
       :module=>:le_grids},
     :new_trap_int=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:new_trap_int,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:new_trap_int,
       :module=>:le_grids},
     :nterp=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nterp,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[],
       :type=>:Integer,
       :code_name=>:nterp,
       :module=>:le_grids}},
   :help=>
    "A concise description of the original GS2 pitch angle and energy grid choices may be found in the article entitled \"Comparison of Initial Value and Eigenvalue Codes for Kinetic Toroidal Plasma Instabilities\" by M. Kotschenreuther, et al., in Computer Physics Communications, Vol. 88, page 128, 1995.\n\nSince then the energy grid has been updated to use the Candy/Waltz grid, but the pitch angle grid remains the same."},
 :dist_fn_knobs=>
  {:description=>"BOUNDARY CONDITIONS",
   :help=>
    "Two important choices are made in this namelist: \n(a) The boundary conditions along the field line\n(b) The form of the adiabatic response (if a species is being modeled as adiabatic)",
   :should_include=>"true",
   :variables=>
    {:gridfac=>
      {:help=>
        "Affects boundary condition at end of theta grid.   Recommended value: 1.",
       :should_include=>"true",
       :description=>"Affects boundary condition at end of theta grid.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :omprimfac=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :boundary_option=>
      {:help=>
        "Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi). Possible values are: \n'zero', 'default', 'unconnected' - Use Kotschenreuther boundary condition at endpoints of theta grid\n'self-periodic', 'periodic', 'kperiod=1' - Each mode is periodic in theta with itself\n'linked' - Twist and shift boundary conditions\n'alternate-zero' - Ignored",
       :should_include=>"true",
       :description=>
        "Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi).",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "zero",
         "unconnected",
         "self-periodic",
         "periodic",
         "kperiod=1",
         "linked",
         "alternate-zero"],
       :module=>:dist_fn},
     :adiabatic_option=>
      {:help=>
        "The form of the adiabatic response (if a species is being modeled as adiabatic). Ignored if there are electrons in the species list.\n**  'no-field-line-average-term'  Adiabatic species has n = Phi.  Appropriate for single-species ETG simulations. \n**  'default'  Same as 'no-field-line-average-term'\n**  'iphi00=0' Same as 'no-field-line-average-term'\n**  'iphi00=1' Same as 'no-field-line-average-term'\n**  'field-line-average-term'  Adiabatic species has n=Phi-< Phi >.  Appropriate for single-species ITG simulations.\n**  'iphi00=2' Same as field-line-average-term'\n**  'iphi00=3' Adiabatic species has n=Phi-< Phi >_y.  Incorrect implementation of field-line-average-term.",
       :should_include=>"true",
       :description=>
        "The form of the adiabatic response (if a species is being modeled as adiabatic).",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "no-field-line-average-term",
         "field-line-average-term",
         "iphi00=0",
         "iphi00=1",
         "iphi00=2",
         "iphi00=3",
         "dimits"],
       :module=>:dist_fn},
     :test=>
      {:help=>"For debugging only.",
       :should_include=>"true",
       :description=>"For debugging",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:dist_fn},
     :g_exb=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :btor_slab=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :apfac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Leave as unity.  For debugging.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:apfac,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:apfac,
       :module=>:dist_fn},
     :driftknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Leave as unity.  For debugging/\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:driftknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:driftknob,
       :module=>:dist_fn},
     :tpdriftknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tpdriftknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-9900000000.0],
       :type=>:Float,
       :code_name=>:tpdriftknob,
       :module=>:dist_fn},
     :poisfac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " If non-zero, quasineutrality is not enforced; poisfac=  (lambda_Debye/rho)**2 \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:poisfac,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:poisfac,
       :module=>:dist_fn},
     :kfilter=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" For debugging only.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:kfilter,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:kfilter,
       :module=>:dist_fn},
     :afilter=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" For debugging only.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:afilter,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:afilter,
       :module=>:dist_fn},
     :nonad_zero=>
      {:should_include=>"true",
       :description=>"If true, new boundary conditions g_wesson=0",
       :help=>
        "If true, new boundary conditions g_wesson=0. If false, g_gs2 = 0",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:nonad_zero,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:nonad_zero,
       :module=>:dist_fn},
     :mult_imp=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Allow different species to have different values of bakdif and fexpr.   Not allowed for nonlinear runs. \n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:mult_imp,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:mult_imp,
       :module=>:dist_fn},
     :def_parity=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" True only allows solutions of single parity.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:def_parity,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:def_parity,
       :module=>:dist_fn},
     :dist_fn_even=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" If def_parity=T, determines allowed parity.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:even,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:even},
     :wfb=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:wfb,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:wfb,
       :module=>:dist_fn},
     :mach=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:mach,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :g_exbfac=>
      {:should_include=>"true",
       :description=>"Multiplies g_exb in the perpendicular shearing term.",
       :help=>
        "Multiplies g_exb in the perpendicular shearing term ''but not'' in the parallel drive term. Use for simulations with purely parallel flow.",
       :tests=>["Tst::FLOAT"],
       :code_name=>:g_exbfac,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :nonad_zero=>
      {:should_include=>"true",
       :description=>
        "If true switches on new parallel boundary condition where h=0 at incoming boundary instead of g=0.",
       :help=>
        "If true switches on new parallel boundary condition where h=0 at incoming boundary instead of g=0.",
       :code_name=>:nonad_zero,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool},
     :g_exb_start_time=>
      {:should_include=>"true",
       :description=>"Flow shear switched on at this time.",
       :help=>"Flow shear switched on at this time.",
       :code_name=>:g_exb_start_time,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :g_exb_start_timestep=>
      {:should_include=>"true",
       :description=>"Flow shear is switched on at this time step.",
       :help=>"Flow shear is switched on at this time step.",
       :code_name=>:g_exb_start_timestep,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :g_exb_error_limit=>
      {:should_include=>"true",
       :description=>
        "With flow shear in single mode, constrains relative error in phi^2.",
       :help=>
        "With flow shear in single mode, constrains relative error in phi^2.",
       :code_name=>:g_exb_error_limit,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float}}},
 :fields_knobs=>
  {:description=>"ALGORITHMIC CHOICES",
   :should_include=>"true",
   :variables=>
    {:field_option=>
      {:help=>
        "The field_option variable controls which time-advance algorithm is used for the linear terms. Allowed values are:                                          \n** 'implicit' Advance linear terms with Kotschenreuther's implicit algorithm.                                                                   \n** 'default'  Same as 'implicit'                                                                                                                     \n** 'explicit' Use second-order Runge-Kutta.  Experimental.                                                                                          \n** 'test' Use for debugging.",
       :should_include=>"true",
       :description=>
        "Controls which time-advance algorithm is used for the linear terms.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>["default", "implicit", "explicit", "test"],
       :module=>:fields},
     :remove_zonal_flows_switch=>
      {:should_include=>"true",
       :description=>"Delete zonal flows at every timestep.",
       :help=>"Delete zonal flows at every timestep.",
       :code_name=>:remove_zonal_flows_switch,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:fields}}},
 :knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:wstar_units=>
      {:help=>
        "Default = .false.  Make the timestep proportional to <math>k_y \\rho</math>.  This can be useful for linear stability calculations that have a wide range of <math>k_y</math> values.  Do not set to true for nonlinear runs.  Be aware that the units of some output quantities can change when wstar_units = .true.",
       :should_include=>"true",
       :description=>
        "For linear runs only. Evolves each k_y with a different timestep.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:run_parameters},
     :fphi=>
      {:help=>
        "Multiplies <math>\\Phi</math> (electrostatic potential) throughout.  Useful for debugging. Non-experts use 1.0",
       :should_include=>"true",
       :description=>"Multiplies Phi (electrostatic potential).",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:run_parameters},
     :fapar=>
      {:help=>
        "Multiplies <math>A_\\parallel</math> throughout.  Set to zero for electrostatic calculations, or unity for electromagnetic.  Set to zero if <math>\\beta = 0</math>.",
       :should_include=>"true",
       :description=>
        "Multiplies A_par. Use 1 for finite beta (electromagnetic), 0 otherwise (electrostatic)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0, 1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:run_parameters},
     :faperp=>
      {:help=>
        "Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead. Defaults to zero.  Do not change this value unless you know what you are doing.",
       :should_include=>"true",
       :description=>
        "Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:run_parameters},
     :fbpar=>
      {:help=>
        "Multiplies <math>B_\\parallel</math> throughout.  Set to zero or unity, depending on whether you want to include physics related to <math>\\delta B_\\parallel</math>. Set to zero if <math>\\beta = 0</math>.",
       :should_include=>"true",
       :description=>
        "Multiplies B_parallel. Use 1 for high beta, 0 otherwise.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[-1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:run_parameters},
     :delt_option=>
      {:help=>
        "Deprecated.  Do not use.  (Use 'check_restart' to get initial timestep from restart file, 'default' otherwise.)",
       :should_include=>"true",
       :description=>"Deprecated.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>["default", "set_by_hand", "check_restart"],
       :module=>:run_parameters},
     :delt=>
      {:help=>
        "Timestep, in units of <math>a/v_{t0}</math>.  For linear runs, this value does not change.  For nonlinear runs, the timestep used by the code will not be allowed to exceed this value.",
       :should_include=>"true",
       :description=>"Time step",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:theta_grid_gridgen},
     :nstep=>
      {:help=>
        "Number of timesteps that will be taken, unless the code stops for some (usually user-specified) reason.",
       :should_include=>"true",
       :description=>"Maximum number of timesteps",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:parameter_scan},
     :eqzip=>
      {:should_include=>"true",
       :description=>
        "True only for secondary/tertiary instability calculations.",
       :help=>
        "Default = .false.   Do not evolve certain <math>k</math> modes in time.  Set this to true only if you know what you are doing. True only for secondary/tertiary instability calculations.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:eqzip,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:eqzip,
       :module=>:run_parameters},
     :margin=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Obsolete.  Fraction of T3E batch job used for finishing up.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:margin,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.05],
       :type=>:Float,
       :code_name=>:margin,
       :module=>:run_parameters},
     :secondary=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        "Default = .false.  Do not set to true unless you know what you are doing.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:secondary,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false.", ".true."],
       :type=>:Fortran_Bool,
       :code_name=>:secondary,
       :module=>:run_parameters},
     :tertiary=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        "Default = .false.  Do not set to true unless you know what you are doing.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:tertiary,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:tertiary,
       :module=>:run_parameters},
     :harris=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        "Default = .false.  Do not set to true unless you know what you are doing.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:harris,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:harris,
       :module=>:run_parameters},
     :avail_cpu_time=>
      {:should_include=>"true",
       :description=>
        "Specify the available wall clock time in seconds. GS2 will exit before this time.",
       :help=>
        "Specify the available wall clock time in seconds. GS2 will exit before this time. This ensures that all the output files are written correctly.  CodeRunner automatically sets this quantity unless it is given the value false.",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:avail_cpu_time,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[10000000000.0],
       :type=>:Float,
       :code_name=>:avail_cpu_time,
       :module=>:run_parameters},
     :eqzip_option=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :code_name=>:eqzip_option,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:run_parameters},
     :neo_test=>
      {:should_include=>"true",
       :description=>
        "Exit after writing out neoclassical quantities of interest.",
       :help=>
        "If true, GS2 exits after writing out neoclassical quantities of interest.",
       :code_name=>:neo_test,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool}}},
 :reinit_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:delt_adj=>
      {:help=>"When the time step needs to be changed, it is adjusted ",
       :should_include=>"true",
       :description=>"When the time step needs to be changed, it is adjusted ",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[2.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:gs2_reinit},
     :delt_minimum=>
      {:help=>"The minimum time step is delt_minimum.",
       :should_include=>"true",
       :description=>"The minimum time step is delt_minimum.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0e-05],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:gs2_reinit},
     :delt_cushion=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.5],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:gs2_reinit},
     :abort_rapid_time_step_change=>
      {:should_include=>"true",
       :description=>"If true (default), exit if time step changes rapidly.",
       :help=>
        "If true (default), exit if time step changes rapidly, that is, if the time step changes at four consecutive time steps.",
       :code_name=>:abort_rapid_time_step_change,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool}}},
 :layouts_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:layout=>
      {:help=>
        "Determines the way the grids are laid out in memory. Rightmost is parallelised first. \n** Can be 'yxles', 'lxyes', 'lyxes', 'lexys' \n** Strongly affects performance on parallel computers\n** In general avoid parallelizing over x. For this reason 'lxyes' is often a good choice.",
       :should_include=>"true",
       :description=>
        "'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["lxyes"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:gs2_layouts},
     :local_field_solve=>
      {:help=>
        "Strongly affects initialization time on some parallel computers. \n**  Recommendation:  Use T on computers with slow communication networks.",
       :should_include=>"true",
       :description=>
        "Strongly affects initialization time on some parallel computers.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:gs2_layouts}}},
 :collisions_knobs=>
  {:description=>"COLLISIONS",
   :should_include=>"true",
   :variables=>
    {:collision_model=>
      {:help=>
        "Collision model used in the simulation. \n\n** ''default'' = pitch angle scattering and energy diffusion\n** ''collisionless'',''none'' = collisionless\n** ''lorentz'' =  pitch angle scattering only\n** ''ediffuse'' = energy diffusion only\n** ''krook'' = use home made krook operator (no reason to use this!)",
       :should_include=>"true",
       :description=>
        "Collision model used in the simulation. Options: 'default', 'none', 'lorentz', 'ediffuse'",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "lorentz",
         "ediffuse",
         "krook",
         "krook-test",
         "lorentz-test",
         "none",
         "collisionless"],
       :module=>:collisions},
     :heating=>
      {:help=>"Set to .true. to compute collisional heating.",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool},
     :conservative=>
      {:help=>"Set to .true. to guarantee exact conservation properties.",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool},
     :conserve_moments=>
      {:help=>
        "Set to .true. to guarantee collision operator conserves momentum and energy.",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool},
     :resistivity=>
      {:should_include=>"true",
       :description=>"",
       :help=>"",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:resistivity,
       :autoscanned_defaults=>[".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :code_name=>:resistivity},
     :adjust=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:adjust,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:adjust},
     :const_v=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:const_v,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:const_v},
     :cfac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        "Factor multiplying FLR terms in collision operator.  1.0 by default.  Set to 0.0 to turn off FLR corrections.",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:cfac,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 1.0],
       :type=>:Float,
       :code_name=>:cfac},
     :hypermult=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:hypermult,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:hypermult},
     :ei_coll_only=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:ei_coll_only,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:ei_coll_only},
     :lorentz_scheme=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :gs2_name=>:lorentz_scheme,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :autoscanned_defaults=>["default"],
       :type=>:String,
       :text_options=>["default", "old"],
       :code_name=>:lorentz_scheme},
     :ediff_scheme=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :gs2_name=>:ediff_scheme,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :autoscanned_defaults=>["default"],
       :type=>:String,
       :text_options=>["default", "old"],
       :code_name=>:ediff_scheme},
     :test_collisions=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:test,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:test},
     :vnfac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:vnfac,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:collisions},
     :etol=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:etol,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:collisions},
     :ewindow=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:ewindow,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:collisions},
     :ncheck=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:ncheck,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:collisions},
     :vnslow=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:vnslow,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:collisions},
     :vary_vnew=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:vary_vnew,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:collisions},
     :etola=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:etola,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:collisions},
     :ewindowa=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:ewindowa,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:collisions}}},
 :hyper_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:hyper_option=>
      {:help=>"'default' is 'none'",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>["default", "none", "visc_only", "res_only", "both"],
       :module=>:hyper},
     :isotropic_shear=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".true."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:hyper},
     :const_amp=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:hyper},
     :include_kpar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:include_kpar,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:include_kpar,
       :module=>:hyper},
     :omega_osc=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:omega_osc,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.4],
       :type=>:Float,
       :code_name=>:omega_osc,
       :module=>:hyper},
     :gridnorm=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:gridnorm,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:gridnorm,
       :module=>:hyper},
     :nexp=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nexp,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[2],
       :type=>:Integer,
       :code_name=>:nexp,
       :module=>:hyper},
     :d_hypervisc=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[-10.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :d_hyperres=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:d_hyperres,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-10.0],
       :type=>:Float,
       :code_name=>:d_hyperres},
     :d_hyper=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:d_hyper,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-10.0],
       :type=>:Float,
       :code_name=>:d_hyper},
     :damp_zonal_only=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :code_name=>:damp_zonal_only,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool}},
   :should_pass=>
    [{:test=>
       "\ncase @hyper_option\nwhen 'none', 'default', nil\n  (not @d_hyperres or @d_hyperres <= 0.0) and (not @d_hypervisc or @d_hypervisc <= 0.0)\nwhen 'visc_only'\t\n  (not @d_hyperres or @d_hyperres <= 0.0) \nwhen 'res_only'\n  (not @d_hypervisc or @d_hypervisc <= 0.0)\nelse true\nend",
      :explanation=>
       "\nYou have entered positive values for the hyper resistivity or hyper diffusivity that are incompatible with hyper_option. They will be set to 0 where appropriate."}],
   :must_pass=>
    [{:test=>
       "\ncase @hyper_option\nwhen 'visc_only'\t\n   @d_hypervisc and @d_hypervisc > 0.0\nwhen 'res_only'\n  @d_hyperres and @d_hyperres > 0.0\nwhen 'both'\n  @d_hyperres and @d_hyperres > 0.0 and @d_hypervisc and @d_hypervisc > 0.0\nelse true\nend",
      :explanation=>
       "\nYou have entered values for the hyper resistivity or hyper diffusivity that are incompatible with hyper_option. If hyper_option specifies hyper diffusivity, hyper resistivity or both, the values for d_hypervisc or d_hyperres should be set positive accordingly."}]},
 :nonlinear_terms_knobs=>
  {:description=>"NONLINEARITY",
   :should_include=>"true",
   :variables=>
    {:nonlinear_mode=>
      {:help=>
        "Should the nonlinear terms be calculated?\n \n** 'none', 'default', 'off':  Do not include nonlinear terms, i.e. run a linear calculation.\n** 'on' Include nonlinear terms.",
       :should_include=>"true",
       :description=>"Include nonlinear terms? ('on','off')",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>["default", "none", "off", "on"],
       :module=>:nonlinear_terms},
     :flow_mode=>
      {:help=>" Experimental\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>["default", "off", "on"],
       :module=>:nonlinear_terms},
     :cfl=>
      {:help=>"The maximum delt < cfl * min(Delta_perp/v_perp)",
       :should_include=>"true",
       :description=>"The maximum delt < cfl * min(Delta_perp/v_perp)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.1],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:nonlinear_terms},
     :p_x=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Ignored.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:p_x,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[6.0],
       :type=>:Float,
       :code_name=>:p_x},
     :p_y=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Ignored.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:p_y,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[6.0],
       :type=>:Float,
       :code_name=>:p_y},
     :p_z=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Ignored.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:p_z,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[6.0],
       :type=>:Float,
       :code_name=>:p_z},
     :zip=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Experts only (for secondary/tertiary calculations). \n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:zip,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[],
       :type=>:Fortran_Bool,
       :code_name=>:zip,
       :module=>:run_parameters},
     :c_par=>
      {:should_include=>"true",
       :description=>nil,
       :help=>"  Ignored.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:C_par,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.1],
       :type=>:Float,
       :code_name=>:C_par},
     :c_perp=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Ignored.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:C_perp,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.1],
       :type=>:Float,
       :code_name=>:C_perp}}},
 :additional_linear_terms_knobs=>
  {:description=>"ADDITIONAL LINEAR TERMS",
   :should_include=>"true",
   :variables=>{}},
 :species_knobs=>
  {:description=>"EVOLVED SPECIES",
   :should_include=>"true",
   :variables=>
    {:nspec=>
      {:help=>"Number of kinetic species evolved.",
       :should_include=>"true",
       :description=>"Number of kinetic species evolved.",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[2],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:species}}},
 :species_parameters=>
  {:description=>"SPECIES PARAMETERS",
   :help=>
    "There should be a separate namelist for each species.  For example, if\nthere are two species, there will be namelists called\nspecies_parameters_1 and species_parameters_2. Charge, mass, density and temperature for each species are relative to some reference species.",
   :enumerator=>{:name=>:nspec, :estimated_value=>5},
   :should_include=>"true",
   :variables=>
    {:z=>
      {:help=>"Charge",
       :should_include=>"true",
       :description=>"Charge",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:fields},
     :mass=>
      {:help=>"Mass",
       :should_include=>"true",
       :description=>"Mass",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :dens=>
      {:help=>"Density\t",
       :should_include=>"true",
       :description=>"Density\t",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :temp=>
      {:help=>"Temperature",
       :should_include=>"true",
       :description=>"Temperature",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:file_utils},
     :tprim=>
      {:help=>
        "Normalised inverse temperature gradient: <math>-1/T (dT/d\\rho)</math>",
       :should_include=>"true",
       :description=>"-1/T (dT/drho)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[6.9],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :fprim=>
      {:help=>
        "Normalised inverse density gradient: <math>-1/n (dn/d\\rho)</math>",
       :should_include=>"true",
       :description=>"-1/n (dn/drho)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[2.2],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :uprim=>
      {:help=>"?",
       :should_include=>"true",
       :description=>"?",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :vnewk=>
      {:help=>"Collisionality parameter: collisionality normalized to v_th/a",
       :should_include=>"true",
       :description=>"collisionality parameter",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :type=>
      {:help=>
        "Type of species:\n** 'ion' Thermal ion species\n** 'default' Same as 'ion'\n** 'electron' Thermal electron species\n** 'e' Same as 'electron'\n** 'beam' Slowing down distribution (Requires advanced_egrid = F)\n** 'slowing_down' Same as 'beam'\n** 'fast' Same as 'beam'\n** 'alpha' Same as 'beam'",
       :should_include=>"true",
       :description=>"Type of species, e.g. 'ion', 'electron', 'beam'",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "ion",
         "electron",
         "e",
         "beam",
         "fast",
         "alpha",
         "slowing-down",
         "trace"],
       :module=>:parameter_scan},
     :dens0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:dens0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:dens0,
       :autoscanned_defaults=>[1.0]},
     :u0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:u0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:u0,
       :autoscanned_defaults=>[1.0]},
     :uprim2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:uprim2,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:uprim2,
       :autoscanned_defaults=>[0.0]},
     :nustar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:nustar,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:nustar,
       :autoscanned_defaults=>[-1.0]},
     :nu=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:nu,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:nu,
       :autoscanned_defaults=>[-1.0]},
     :nu_h=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:nu_h,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:nu_h,
       :autoscanned_defaults=>[0.0]},
     :tperp0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tperp0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:tperp0,
       :autoscanned_defaults=>[0.0]},
     :tpar0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tpar0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :code_name=>:tpar0,
       :autoscanned_defaults=>[0.0]}}},
 :dist_fn_species_knobs=>
  {:description=>"",
   :should_include=>"true",
   :enumerator=>{:name=>:nspec, :estimated_value=>5},
   :variables=>
    {:fexpr=>
      {:help=>
        "Temporal implicitness parameter. Any value other than 0.5 adds numerical dissipation.\n** Recommended value: 0.48",
       :should_include=>"true",
       :description=>
        "Temporal implicitness parameter. Recommended value: 0.48",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>["real(fexp_out)"],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :fexpi=>
      {:help=>nil,
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>["aimag(fexp_out)"],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :bakdif=>
      {:help=>
        "Spatial implicitness parameter. Any value greater than 0 adds numerical dissipation (usually necessary).\n** Recommended value: 0.05",
       :should_include=>"true",
       :description=>"Spatial implicitness parameter. Recommended value: 0.05",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>["bakdif_out"],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :bd_exp=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:bd_exp,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :code_name=>:bd_exp,
       :autoscanned_defaults=>["bd_exp_out"],
       :module=>:dist_fn}}},
 :init_g_knobs=>
  {:description=>"INITIAL CONDITIONS",
   :should_include=>"true",
   :variables=>
    {:refac=>
      {:help=>"Used in rare cases.",
       :should_include=>"true",
       :description=>"Used in rare cases.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :imfac=>
      {:help=>"Used in rare cases.",
       :should_include=>"true",
       :description=>"Used in rare cases.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:dist_fn},
     :den0=>
      {:help=>"Parameters for setting up special initial conditions.",
       :should_include=>"true",
       :description=>"Parameters for setting up special initial conditions.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :den1=>
      {:help=>"Parameters for setting up special initial conditions.",
       :should_include=>"true",
       :description=>"Parameters for setting up special initial conditions.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :upar0=>
      {:help=>"Parameters for setting up special initial conditions.",
       :should_include=>"true",
       :description=>"Parameters for setting up special initial conditions.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :upar1=>
      {:help=>"Parameters for setting up special initial conditions.",
       :should_include=>"true",
       :description=>"Parameters for setting up special initial conditions.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :tpar0=>
      {:help=>"  Parameters for setting up special initial conditions.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :tperp0=>
      {:help=>"  Parameters for setting up special initial conditions.\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[0.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :chop_side=>
      {:help=>"Rarely needed.  Forces asymmetry into initial condition.",
       :should_include=>"true",
       :description=>"Rarely needed. Forces asymmetry into initial condition.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[".true."],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String},
     :width0=>
      {:help=>
        "Initial perturbation has Gaussian envelope in theta, with width width0",
       :should_include=>"true",
       :description=>
        "Initial perturbation has Gaussian envelope in theta with width width0,",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[-3.5],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:init_g},
     :phiinit=>
      {:help=>
        "Average amplitude of initial perturbation of each Fourier mode.",
       :should_include=>"true",
       :description=>
        "Average amplitude of initial perturbation of each Fourier mode.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :restart_file=>
      {:help=>"Base of filenames with restart data.",
       :should_include=>"true",
       :description=>"Base of filenames with restart data.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>
        ["file",
         "trim(restart_dir)//trim(restart_file)",
         "trim(run_name)//\".nc\""],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:init_g},
     :ginit_option=>
      {:help=>
        "Sets the way that the distribution function is initialized. There are many possible choices.\n**  'default'\n**  'noise'  This is the  recommended selection.  Pretty random.\n**  'test1'\n**  'xi'\n**  'xi2'\n**  'zero'\n**  'test3'\n**  'convect'\n**  'rh'\n**  'many'\n**  'small'\n**  'file'\n**  'cont'\n**  'kz0'  initialise only with k_parallel=0\n**  'nl'\n**  'nl2'\n**  'nl3'\n**  'nl4'\n**  'nl5'\n**  'nl6'\n**  'gs'\n**  'kpar'\n**  'zonal_only'  Restart but set all non-zonal components of the potential and the distribution function to 0. Noise can be added to these other components by setting iphiinit > 0.\n**  'single_parallel_mode'  Initialise only with a single parallel mode specified by either ikpar_init for periodic boundary conditions or kpar_init for linked boundary conditions. Intended for linear calculations.\n**  'all_modes_equal'  Initialise with every single parallel and perpendicular mode given the same amplitude. Intended for linear calculations.",
       :should_include=>"true",
       :description=>
        "Sets the way that the distribution function is initialized.",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>["default"],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :text_options=>
        ["default",
         "noise",
         "test1",
         "xi",
         "xi2",
         "zero",
         "test3",
         "convect",
         "rh",
         "many",
         "small",
         "file",
         "cont",
         "kz0",
         "nl",
         "nl2",
         "nl3",
         "nl3r",
         "nl4",
         "nl5",
         "nl6",
         "nl7",
         "alf",
         "gs",
         "kpar",
         "smallflat",
         "harris",
         "recon",
         "recon3",
         "zonal_only",
         "single_parallel_mode",
         "all_modes_equal"]},
     :zf_init=>
      {:help=>
        " Amplitude of initial zonal flow perturbations relative to other modes\n",
       :should_include=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :left=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Chop out left side in theta. \n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:left,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:left},
     :ikk=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Used only for secondary/tertiary calculations.\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:ikk,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[],
       :type=>:Integer,
       :code_name=>:ikk},
     :itt=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Used only for secondary/tertiary calculations.\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:itt,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[],
       :type=>:Integer,
       :code_name=>:itt},
     :scale=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Allows rescaling of amplitudes for restarts.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:scale,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:scale},
     :tstart=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Force t=tstart at beginning of run.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tstart,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-1.0, 0.0],
       :type=>:Float,
       :code_name=>:tstart,
       :module=>:init_g},
     :even=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Sometimes initial conditions have definite parity; this picks the parity in those cases.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:even,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:even,
       :module=>:dist_fn},
     :tpar1=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Parameters for setting up special initial conditions.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tpar1,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:tpar1},
     :tperp1=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Parameters for setting up special initial conditions.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tperp1,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:tperp1},
     :den2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Parameters for setting up special initial conditions.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:den2,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:den2},
     :upar2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Parameters for setting up special initial conditions.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:upar2,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:upar2},
     :tpar2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Parameters for setting up special initial conditions.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tpar2,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:tpar2},
     :tperp2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Parameters for setting up special initial conditions.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tperp2,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:tperp2},
     :dphiinit=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:dphiinit,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:dphiinit},
     :apar0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:apar0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:apar0},
     :restart_dir=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :code_name=>:restart_dir,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String},
     :new_field_init=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:new_field_init,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:init_g},
     :phiinit0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:phiinit0,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :a0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:a0,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:kt_grids_single},
     :b0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:b0,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :null_phi=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:null_phi,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool},
     :null_bpar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:null_bpar,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool},
     :null_apar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:null_apar,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool},
     :adj_spec=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:adj_spec,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :eq_type=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :code_name=>:eq_type,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String},
     :prof_width=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:prof_width,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :eq_mode_u=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:eq_mode_u,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :eq_mode_n=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:eq_mode_n,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :input_check_recon=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:input_check_recon,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool},
     :nkxy_pt=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["true"],
       :code_name=>:nkxy_pt,
       :must_pass=>
        [{:test=>"kind_of? Complex",
          :explanation=>"This variable must be a complex number."}],
       :type=>:Complex},
     :ukxy_pt=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["true"],
       :code_name=>:ukxy_pt,
       :must_pass=>
        [{:test=>"kind_of? Complex",
          :explanation=>"This variable must be a complex number."}],
       :type=>:Complex},
     :ikkk=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:ikkk,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :ittt=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:ittt,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :phifrac=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:phifrac,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :ikpar_init=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :code_name=>:ikpar_init,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :kpar_init=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :code_name=>:kpar_init,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :ikx_init=>
      {:should_include=>"true",
       :description=>
        "Only initialise noise for the kx mode indexed by ikx_index.",
       :help=>
        "Only initialise a single kx with noise. This input parameter is used when noise is being initialised. If specified, i.e. if it is set greater than zero, noise will only be initialised for itheta0 = ikx_index, i.e. for the mode indexed by ikx_index. this is useful for linear runs with flow shear, to track the evolution of a single Lagrangian mode.",
       :tests=>["Tst::INT"],
       :code_name=>:ikx_init,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer},
     :phiinit_rand=>
      {:should_include=>"true",
       :description=>"Amplitude of random perturbation for ginit_recon3",
       :help=>"Amplitude of random perturbation for ginit_recon3 (R Numata)",
       :tests=>["Tst::FLOAT"],
       :code_name=>:phiinit_rand,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float},
     :phiamp=>
      {:should_include=>"true",
       :description=>"For the Orszag-Tang 2D vortex problem.",
       :help=>"Used in initialization for the Orszag-Tang 2D vortex problem.",
       :tests=>["true"],
       :code_name=>:phiamp,
       :must_pass=>
        [{:test=>"kind_of? Complex",
          :explanation=>"This variable must be a complex number."}],
       :type=>:Complex},
     :aparamp=>
      {:should_include=>"true",
       :description=>
        "Used in initialization for the Orszag-Tang 2D vortex problem",
       :help=>"Used in initialization for the Orszag-Tang 2D vortex problem",
       :tests=>["true"],
       :code_name=>:aparamp,
       :must_pass=>
        [{:test=>"kind_of? Complex",
          :explanation=>"This variable must be a complex number."}],
       :type=>:Complex},
     :force_single_kpar=>
      {:should_include=>"true",
       :description=>
        "Set other parallel modes to zero when kpar_init specified.",
       :help=>
        "When initialised with single parallel mode, sets other modes to zero at every time step.",
       :code_name=>:force_single_kpar,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool}}},
 :gs2_diagnostics_knobs=>
  {:description=>"DIAGNOSTICS",
   :should_include=>"true",
   :variables=>
    {:print_flux_line=>
      {:help=>"Instantaneous fluxes output to screen every nwrite timesteps",
       :should_include=>"true",
       :description=>"Instantaneous fluxes output to screen",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_nl_flux=>
      {:help=>"Write nonlinear fluxes as a function of time.",
       :should_include=>"true",
       :description=>"Write nonlinear fluxes as a function of time.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :print_line=>
      {:help=>
        "Estimated frequencies and output to the screen/stdout every nwrite timesteps",
       :should_include=>"true",
       :description=>
        "Estimated frequencies and growth rates to the screen/stdout",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false.", ".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_verr=>
      {:help=>"Write velocity space diagnostics to '.lpc' and '.verr' files",
       :should_include=>"true",
       :description=>
        "Write velocity space diagnostics to '.lpc' and '.verr' files",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_g=>
      {:help=>"Write the distribution function to the '.dist' (NetCDF?)",
       :should_include=>"true",
       :description=>
        "Write the distribution function to the '.dist' (NetCDF?)",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_line=>
      {:help=>
        "If (write_ascii = T) write estimated frequencies and growth rates to the output file (usually runname.out) every nwrite steps.",
       :should_include=>"true",
       :description=>
        "If (write_ascii = T) write estimated frequencies and growth rates to the output file",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".true."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_gyx=>
      {:help=>"Write dist fn at a given physical spacial point to a file",
       :should_include=>"true",
       :description=>
        "Write dist fn at a given physical spacial point to a file",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_hrate=>
      {:help=>
        "Write heating rate, collisonal entropy generation etc to '.heat'",
       :should_include=>"true",
       :description=>
        "Write heating rate, collisonal entropy generation etc to '.heat'",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_final_epar=>
      {:help=>
        "If (write_ascii = T) E_parallel(theta) written to runname.eigenfunc\n** Write to runname.out.nc even if (write_ascii = F)",
       :should_include=>"true",
       :description=>
        "If (write_ascii = T) E_parallel(theta) written to runname.eigenfunc",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_avg_moments=>
      {:help=>
        "Ignored unless grid_option='box'\n** Flux surface averaged low-order moments of g written to runname.out.nc\n** If (write_ascii = T) flux surface averaged low-order moments of g written to runname.moments",
       :should_include=>"true",
       :description=>
        "Write flux surface averaged low-order moments of g to runname.out.nc and runname.moments (if write_ascii = T)",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_lorentzian=>
      {:help=>"Frequency Sweep Data",
       :should_include=>"true",
       :description=>"Frequency Sweep Data",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_omega=>
      {:help=>
        "If (write_ascii = T) instantaneous omega to output file every nwrite timesteps. Very heavy output.",
       :should_include=>"true",
       :description=>
        "If (write_ascii = T) instantaneous omega to output file. Very heavy output",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_omavg=>
      {:help=>
        "If (write_ascii = T) time-averaged frequencies written to runname.out every nwrite timesteps.\n** Average is over navg steps.",
       :should_include=>"true",
       :description=>
        "If (write_ascii = T) time-averaged growth rate and frequency to the output file.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_eigenfunc=>
      {:help=>
        "If (write_ascii = T) Normalized Phi(theta) written to runname.eigenfunc\n** Write to runname.out.nc even if (write_ascii = F)",
       :should_include=>"true",
       :description=>
        "If (write_ascii = T) Normalized phi written to runname.eigenfunc",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_final_fields=>
      {:help=>
        "If (write_ascii = T) Phi(theta) written to runname.fields\n** Write to runname.out.nc even if (write_ascii = F)",
       :should_include=>"true",
       :description=>"If (write_ascii = T) Phi(theta) written to '.fields'",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_final_moments=>
      {:help=>
        "If (write_ascii = T) low-order moments of g written to runname.moments and int dl/B averages of low-order moments of g written to  runname.amoments\n** Write to runname.out.nc even if (write_ascii = F)",
       :should_include=>"true",
       :description=>"write final n, T",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_parity=>
      {:help=>"Writes parities in dist fn and particle fluxes",
       :should_include=>"true",
       :description=>"Writes parities in dist fn and particle fluxes",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :nsave=>
      {:help=>"Write restart files every nsave timesteps",
       :should_include=>"true",
       :description=>"Write restart files every nsave timesteps",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[-1],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:gs2_diagnostics},
     :nwrite=>
      {:help=>"Output diagnostic data every nwrite timesteps.",
       :should_include=>"true",
       :description=>"Output diagnostic data every nwrite",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[100],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:gs2_diagnostics},
     :navg=>
      {:help=>"Any time averages performed over navg timesteps.",
       :should_include=>"true",
       :description=>"Any time averages performed over navg",
       :tests=>["Tst::INT"],
       :autoscanned_defaults=>[100],
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:gs2_diagnostics},
     :omegatol=>
      {:help=>
        "In linear runs GS2 will exit if the growth rate has converged to an accuracy of one part in 1/omegatol. Set negative to switch off this feature.",
       :should_include=>"true",
       :description=>
        "The convergence has to be better than one part in 1/omegatol",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:gs2_diagnostics},
     :omegatinst=>
      {:help=>"Recommended value: 500.",
       :should_include=>"true",
       :description=>"Recommended value: 500.",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[1.0],
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:gs2_diagnostics},
     :save_for_restart=>
      {:help=>
        "If true then restart files written to the local folder and the simulation can be restarted from the point it ended.\n** Restart files written to restart_file.PE#.  \n** Recommended for nonlinear runs.",
       :should_include=>"true",
       :description=>"Write restart files.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :autoscanned_defaults=>[".false."],
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a Fortran boolean. (In Ruby this is represented as a string: e.g. '.true.') "}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_flux_line=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " If (write_ascii = T) instantaneous fluxes output to runname.out every nwrite timesteps\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_flux_line,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:write_flux_line,
       :module=>:gs2_diagnostics},
     :write_ascii=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" If true, some data is written to runname.out\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_ascii,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:write_ascii,
       :module=>:gs2_diagnostics},
     :write_kpar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Spectrum in k_parallel calculated and written.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_kpar,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_kpar,
       :module=>:gs2_diagnostics},
     :write_gs=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_gs,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_gs,
       :module=>:gs2_diagnostics},
     :write_gg=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_gg,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_gg,
       :module=>:gs2_diagnostics},
     :write_lpoly=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_lpoly,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_lpoly,
       :module=>:gs2_diagnostics},
     :write_fields=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_fields,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_fields,
       :module=>:gs2_diagnostics},
     :write_final_antot=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " If (write_ascii = T) Sources for Maxwell eqns. written to runname.antot\n** Write to runname.out.nc even if (write_ascii = F)\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_final_antot,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_final_antot,
       :module=>:gs2_diagnostics},
     :write_cerr=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_cerr,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_cerr,
       :module=>:gs2_diagnostics},
     :write_max_verr=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_max_verr,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_max_verr,
       :module=>:gs2_diagnostics},
     :nmovie=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nmovie,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[1000],
       :type=>:Integer,
       :code_name=>:nmovie,
       :module=>:gs2_diagnostics},
     :igomega=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Theta index at which frequencies are calculated.\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:igomega,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0],
       :type=>:Integer,
       :code_name=>:igomega,
       :module=>:gs2_diagnostics},
     :exit_when_converged=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " When the frequencies for each k have converged, the run will stop.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:exit_when_converged,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:exit_when_converged,
       :module=>:gs2_diagnostics},
     :write_full_moments_notgc=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_full_moments_notgc,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_full_moments_notgc,
       :module=>:gs2_diagnostics},
     :write_cross_phase=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_cross_phase,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[],
       :type=>:Fortran_Bool,
       :code_name=>:write_cross_phase,
       :module=>:gs2_diagnostics},
     :dump_check1=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Field-line avg of Phi written to dump.check1\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:dump_check1,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:dump_check1,
       :module=>:gs2_diagnostics},
     :dump_check2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Apar(kx, ky, igomega) written to dump.check2\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:dump_check2,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:dump_check2,
       :module=>:gs2_diagnostics},
     :dump_fields_periodically=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Phi, Apar, Bpar written to dump.fields.t=(time).  This is expensive!\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:dump_fields_periodically,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:dump_fields_periodically,
       :module=>:gs2_diagnostics},
     :make_movie=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:make_movie,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:make_movie,
       :module=>:gs2_diagnostics},
     :write_phi_over_time=>
      {:should_include=>"true",
       :description=>"Write entire Phi field to NetCDF file every nwrite.",
       :help=>
        "If this variable is set to true then the entire field Phi will be written to the NetCDF file every nwrite. Useful for making films. This can cause the NetCDF file to be huge, if resolution is large or nwrite is small.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:write_phi_over_time,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_apar_over_time=>
      {:should_include=>"true",
       :description=>
        "Write entire A_parallel field to NetCDF file every nwrite.",
       :help=>
        "If this variable is set to true then the entire field A_parallel will be written to the NetCDF file every nwrite. This can cause the NetCDF file to be huge, if resolution is large or nwrite is small.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:write_apar_over_time,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_bpar_over_time=>
      {:should_include=>"true",
       :description=>
        "Write entire B_parallel field to NetCDF file every nwrite.",
       :help=>
        "If this variable is set to true then the entire field B_parallel will be written to the NetCDF file every nwrite. Useful for making films. This can cause the NetCDF file to be huge, if resolution is large or nwrite is small.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:write_bpar_over_time,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_symmetry=>
      {:should_include=>"true",
       :description=>"Test the symmetry properties of the GK eqn.",
       :help=>
        "Switch on a diagnostic to test the symmetry properties of the GK eqn.  It calculates the momentum flux as a function of vpar, theta, and time.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:write_symmetry,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :save_distfn=>
      {:should_include=>"true",
       :description=>"Save dist_fn with lots of detail.",
       :help=>
        "If true, saves the restart files with name 'rootname.nc.dfn.<proc>' with lots of extra detail about the dist function --- velocity space grids and so on, when GS2 exits.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:save_distfn,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_correlation_extend=>
      {:should_include=>"true",
       :description=>"Extend domain of correlation function calculation.",
       :help=>
        "If used in conjunction with write_correlation, extends the length of <math>\\Delta \\theta</math> for which the correlation function is calculated.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:write_correlation_extend,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :nwrite_mult=>
      {:should_include=>"true",
       :description=>
        "Large datasets written every nwrite_mult * nwrite timesteps.",
       :help=>
        "Multiplies nwrite to determine when large/expensive to calculate datasets such as the parallel correlation function are written to file.",
       :tests=>["Tst::INT"],
       :code_name=>:nwrite_mult,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:gs2_diagnostics},
     :write_correlation=>
      {:should_include=>"true",
       :description=>"Write parallel correlation.",
       :help=>
        "Write correlation function diagnostic... shows parallel correlation as a function of ky. See arXiv 1104.4514.",
       :tests=>["Tst::FORTRAN_BOOL"],
       :code_name=>:write_correlation,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:gs2_diagnostics},
     :write_moments=>
      {:should_include=>"true",
       :description=>"",
       :help=>"",
       :code_name=>:write_moments,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool},
     :write_final_db=>
      {:should_include=>"true",
       :description=>"Write final delta B.",
       :help=>"Write final delta B.",
       :code_name=>:write_final_db,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool}},
   :help=>
    "Controls what information is output by GS2 during and at the end of a simulation."},
 :testgridgen=>
  {:description=>"",
   :should_include=>true,
   :variables=>
    {:source=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:source,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:source,
       :module=>:dist_fn},
     :gsource=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :gs2_name=>:gsource,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :autoscanned_defaults=>["eik6.out"],
       :type=>:String,
       :code_name=>:gsource},
     :nthetaout=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nthetaout,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[32],
       :type=>:Integer,
       :code_name=>:nthetaout},
     :nlambdaout=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nlambdaout,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[20],
       :type=>:Integer,
       :code_name=>:nlambdaout},
     :nperiodout=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nperiodout,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[2],
       :type=>:Integer,
       :code_name=>:nperiodout},
     :npadd=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:npadd,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[2],
       :type=>:Integer,
       :code_name=>:npadd,
       :module=>:theta_grid_gridgen},
     :alknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:alknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 0.1],
       :type=>:Float,
       :code_name=>:alknob,
       :module=>:theta_grid_gridgen},
     :epsknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:epsknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[],
       :type=>:Float,
       :code_name=>:epsknob,
       :module=>:theta_grid_gridgen},
     :bpknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:bpknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0e-08],
       :type=>:Float,
       :code_name=>:bpknob,
       :module=>:theta_grid_gridgen},
     :extrknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:extrknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 10.0],
       :type=>:Float,
       :code_name=>:extrknob,
       :module=>:theta_grid_gridgen},
     :smoothknob=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:smoothknob,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:smoothknob},
     :nfinegrid=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nfinegrid,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[200],
       :type=>:Integer,
       :code_name=>:nfinegrid},
     :thetamax=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:thetamax,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:thetamax,
       :module=>:theta_grid_gridgen},
     :deltaw=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:deltaw,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:deltaw,
       :module=>:theta_grid_gridgen},
     :widthw=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:widthw,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:widthw,
       :module=>:theta_grid_gridgen},
     :tension=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tension,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:tension,
       :module=>:theta_grid_gridgen},
     :gingrid=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::STRING"],
       :gs2_name=>:gingrid,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :autoscanned_defaults=>["gingrid"],
       :type=>:String,
       :code_name=>:gingrid},
     :screenout=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:screenout,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:screenout},
     :auto_width=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:auto_width,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:auto_width},
     :cv_fraction=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:cv_fraction,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.6],
       :type=>:Float,
       :code_name=>:cv_fraction},
     :delth_max=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:delth_max,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.5],
       :type=>:Float,
       :code_name=>:delth_max},
     :max_autoiter=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:max_autoiter,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[3],
       :type=>:Integer,
       :code_name=>:max_autoiter},
     :three_dim=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:three_dim,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:three_dim},
     :iperiod=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:iperiod,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[1],
       :type=>:Integer,
       :code_name=>:iperiod}}},
 :driver=>
  {:description=>"",
   :should_include=>true,
   :variables=>
    {:amplitude=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Amplitude of Langevin antenna.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:amplitude,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:amplitude,
       :module=>:antenna},
     :w_antenna=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Frequency of Langevin antenna.\n",
       :tests=>["true"],
       :gs2_name=>:w_antenna,
       :must_pass=>
        [{:test=>"kind_of? Complex",
          :explanation=>"This variable must be a complex number."}],
       :autoscanned_defaults=>[Complex(1.0,0.0)],
       :type=>:Complex,
       :code_name=>:w_antenna},
     :nk_stir=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Number of independent Fourier modes driven by antenna.\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:nk_stir,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[1],
       :type=>:Integer,
       :code_name=>:nk_stir},
     :write_antenna=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Write antenna amplitudes to ASCII file for debugging.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:write_antenna,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:write_antenna},
     :ant_off=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Overrides all and turns off antenna if true.\n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:ant_off,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".false."],
       :type=>:Fortran_Bool,
       :code_name=>:ant_off},
     :w_dot=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:w_dot,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:w_dot},
     :driver_t0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:t0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-1.0, 100.0],
       :type=>:Float,
       :code_name=>:t0},
     :restarting=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:restarting,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[],
       :type=>:Fortran_Bool,
       :code_name=>:restarting}}},
 :stir=>
  {:description=>"",
   :should_include=>true,
   :variables=>
    {:stir_kx=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Mode number for stirring\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:kx,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[1],
       :type=>:Integer,
       :code_name=>:kx},
     :stir_ky=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Mode number for stirring\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:ky,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[1],
       :type=>:Integer,
       :code_name=>:ky},
     :stir_kz=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Mode number for stirring\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:kz,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[1],
       :type=>:Integer,
       :code_name=>:kz},
     :stir_travel=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Launches traveling wave (or standing wave if F). \n",
       :tests=>["Tst::FORTRAN_BOOL"],
       :gs2_name=>:travel,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :autoscanned_defaults=>[".true."],
       :type=>:Fortran_Bool,
       :code_name=>:travel},
     :stir_a=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Initial amplitude of right-moving component. It is not necessary to set a and b unless you are\ndoing restarts, which are rather clunky at the moment with the antenna included. \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:a,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-1.0, 0.0],
       :type=>:Float,
       :code_name=>:a},
     :stir_b=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Initial amplitude of left-moving component. It is not necessary to set a and b unless you are\ndoing restarts, which are rather clunky at the moment with the antenna included. \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:b,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-1.0, 0.0],
       :type=>:Float,
       :code_name=>:b}}},
 :source_knobs=>
  {:description=>"",
   :should_include=>true,
   :variables=>
    {:t0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Turn on any artificial sources after time t0.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:t0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[-1.0, 100.0],
       :type=>:Float,
       :code_name=>:t0,
       :module=>:dist_fn},
     :omega0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Frequency of non-standard source (if selected above). \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:omega0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:omega0,
       :module=>:dist_fn},
     :gamma0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Growth rate of non-standard source (if selected above). \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:gamma0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:gamma0,
       :module=>:dist_fn},
     :source0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Amplitude of non-standard source (if selected above). \n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:source0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[1.0],
       :type=>:Float,
       :code_name=>:source0,
       :module=>:dist_fn},
     :phi_ext=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Amplitude of external Phi added as source term.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:phi_ext,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:phi_ext,
       :module=>:dist_fn},
     :source_option=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " \n** 'source_option_full' Solve GK equation in standard form (with no artificial sources)\n** 'default' Same as 'source_option_full' \n** 'zero' The GK distribution function will be advanced non-self-consistently.\n** 'sine' The GK distribution function will be advanced non-self-consistently.\n** 'cosine'The GK distribution function will be advanced non-self-consistently.\n** 'test1' The GK distribution function will be advanced non-self-consistently.\n** 'phiext_full' Solve GK equation with additional source proportional to phi_ext*F_0.  \n** 'test2_full' Solve GK equation with additional developmental sources included.  Experts only.\n** 'convect_full' Solve GK equation with additional developmental sources included.  Experts only.\n** 'test1' The GK distribution function will be advanced non-self-consistently.\n",
       :tests=>["Tst::STRING"],
       :gs2_name=>:source_option,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :autoscanned_defaults=>["default"],
       :type=>:String,
       :text_options=>
        ["default",
         "full",
         "zero",
         "sine",
         "cosine",
         "test1",
         "hm",
         "phiext_full",
         "test2_full",
         "convect_full",
         "neo"],
       :code_name=>:source_option,
       :module=>:dist_fn}}},
 :kt_grids_range_parameters=>
  {:description=>"",
   :should_include=>"@grid_option=='range'",
   :variables=>
    {:naky=>
      {:should_include=>"true",
       :description=>"The number of 'actual' ky modes.",
       :help=>"The number of 'actual' ky modes.",
       :tests=>["Tst::INT"],
       :gs2_name=>:naky,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0, 1],
       :type=>:Integer,
       :code_name=>:naky,
       :module=>:kt_grids_specified},
     :ntheta0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Number of theta_0 (kx) modes\n",
       :tests=>["Tst::INT"],
       :gs2_name=>:ntheta0,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0, 1],
       :type=>:Integer,
       :code_name=>:ntheta0,
       :module=>:kt_grids_specified},
     :aky_min=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Lower limit of (ky rho) range.  Should set to something other than zero.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:aky_min,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:aky_min,
       :module=>:kt_grids_range},
     :aky_max=>
      {:should_include=>"true",
       :description=>nil,
       :help=>
        " Upper limit of (ky rho) range.  Should set to something other than zero.\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:aky_max,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0, 1.0],
       :type=>:Float,
       :code_name=>:aky_max,
       :module=>:kt_grids_range},
     :theta0_min=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Lower limit of theta_0 range\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:theta0_min,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:theta0_min,
       :module=>:kt_grids_range},
     :theta0_max=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" Upper limit of theta_0 range\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:theta0_max,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:theta0_max,
       :module=>:kt_grids_range},
     :akx_min=>
      {:should_include=>"true",
       :description=>
        "Min kx for periodic finite kx ballooning space runs with shat=0",
       :help=>
        "Min kx for periodic finite kx ballooning space runs with <math>\\hat{s}=0</math>.",
       :tests=>["Tst::FLOAT"],
       :code_name=>:akx_min,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:kt_grids_range},
     :akx_max=>
      {:should_include=>"true",
       :description=>
        "Max kx for periodic finite kx ballooning space runs with shat=0.",
       :help=>
        "Max kx for periodic finite kx ballooning space runs with shat=0.",
       :tests=>["Tst::FLOAT"],
       :code_name=>:akx_max,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:kt_grids_range}}},
 :kt_grids_specified_parameters=>
  {:description=>"",
   :should_include=>"@grid_option=='specified'",
   :variables=>
    {:naky=>
      {:should_include=>"true",
       :description=>"Number of ky values evolved",
       :help=>
        "Number of ky values evolved. Total number of modes evolved = max(naky, ntheta0). Also set up the appropriate number of kt_grids_specified_element_i namelists.",
       :tests=>["Tst::INT"],
       :gs2_name=>:naky,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0, 1],
       :type=>:Integer,
       :code_name=>:naky,
       :module=>:kt_grids_specified},
     :ntheta0=>
      {:should_include=>"true",
       :description=>"Number of theta0 values.",
       :help=>
        "Number of theta0 values. Total number of modes evolved = max(naky, ntheta0). Also set up the appropriate number of kt_grids_specified_element_i namelists.",
       :tests=>["Tst::INT"],
       :gs2_name=>:ntheta0,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0, 1],
       :type=>:Integer,
       :code_name=>:ntheta0,
       :module=>:kt_grids_specified},
     :nx=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:nx,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0],
       :type=>:Integer,
       :code_name=>:nx,
       :module=>:kt_grids_specified},
     :ny=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:ny,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :autoscanned_defaults=>[0],
       :type=>:Integer,
       :code_name=>:ny,
       :module=>:kt_grids_specified}}},
 :kt_grids_specified_element=>
  {:description=>"SPECIFIX FOURIER MODE",
   :should_include=>"@grid_option=='specified'",
   :variables=>
    {:aky=>
      {:should_include=>"true",
       :description=>nil,
       :help=>" ky rho\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:aky,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.2, 0.4],
       :type=>:Float,
       :code_name=>:aky,
       :module=>:kt_grids_single},
     :theta0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>"  theta_0\n",
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:theta0,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:theta0,
       :module=>:kt_grids_single},
     :akx=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:akx,
       :must_pass=>
        [{:test=>"kind_of? Float or kind_of? Integer",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :autoscanned_defaults=>[0.0],
       :type=>:Float,
       :code_name=>:akx,
       :module=>:kt_grids_single}},
   :help=>
    "There should be a separate namelist for each Fourier mode. For example, if there are two modes, there will be namelists called kt_grids_specified_element_1 and kt_grids_specified_element_2."},
 :kt_grids_xbox_parameters=>
  {:description=>"", :should_include=>"@grid_option=='xbox'", :variables=>{}},
 :gs2_flux_knobs=>{:description=>"", :should_include=>true, :variables=>{}},
 :flux_target=>{:description=>"", :should_include=>true, :variables=>{}},
 :parameter_scan_knobs=>
  {:description=>"",
   :should_include=>"true",
   :variables=>
    {:scan_type=>
      {:should_include=>"true",
       :description=>"Specifies the way that the parameter scan is conducted.",
       :help=>
        "Specifies the way that the parameter scan is conducted. Possible values are:\n** 'none' -- do not conduct a parameter scan (default)\n** 'range' --  vary parameter in constant increments between 2 values: par_start and par_end. The step size is given by par_inc.\n** 'target' --  start with the parameter at par_start, and then change the parameter by par_inc until the target parameter has reached the target value\n** 'root_finding' -- the same as target, but the increment is changed intelligently using a Newton-like method.",
       :code_name=>:scan_type,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:parameter_scan},
     :scan_par=>
      {:should_include=>"true",
       :description=>"Specify the parameter to be varied.",
       :help=>
        "Specify the parameter to be varied.  If the parameter pertains to a species, the scan_spec must be specified as well.",
       :code_name=>:scan_par,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:parameter_scan},
     :target_par=>
      {:should_include=>"true",
       :description=>
        "If the scan is being run in 'target' or 'root_finding' mode, specifies the target parameter.",
       :help=>
        "If the scan is being run in 'target' or 'root_finding' mode, specifies the target parameter. \n**Possible values are 'hflux_tot', 'momflux_tot', 'phi2_tot'.",
       :code_name=>:target_par,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:parameter_scan},
     :par_start=>
      {:should_include=>"true",
       :description=>"Specifies the starting value for the parameter scan.",
       :help=>"Specifies the starting value for the parameter scan.",
       :code_name=>:par_start,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:parameter_scan},
     :par_end=>
      {:should_include=>"true",
       :description=>
        "If the scan is being run in 'range' mode, specifies the value of the parameter that will be reached.",
       :help=>
        "If the scan is being run in 'range' mode, specifies the value of the parameter that will be reached.",
       :code_name=>:par_end,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:parameter_scan},
     :par_inc=>
      {:should_include=>"true",
       :description=>
        "Specifies the amount by which the parameter is varied at one go.",
       :help=>
        "If the parameter scan is being run in 'range' or 'target' modes, specifies the amount by which the parameter is varied at one go.",
       :code_name=>:par_inc,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:parameter_scan},
     :inc_con=>
      {:should_include=>"true",
       :description=>"Specifies the condition for incrementing the parameter.",
       :help=>
        "Specifies the condition for incrementing the parameter. Possible values are:\n** 'n_timesteps' -- change the parameter after a given number of time steps\n** 'delta_t' -- change the parameter after an elapsed time\n** 'saturated' -- change the parameter after the simulation has reached a saturated state (determined using the target parameter) at the current value of the parameter",
       :code_name=>:inc_con,
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}],
       :type=>:String,
       :module=>:parameter_scan},
     :nstep_init=>
      {:should_include=>"true",
       :description=>
        "The parameter will not be changed until nstep_init have elapsed from the beginning of the simulation.",
       :help=>
        "When the increment condition is 'n_timesteps' or 'saturated',  the parameter will not be changed until nstep_init have elapsed from the beginning of the simulation. Note that if the simulation is restarted, this parameter will measure from the restart.",
       :code_name=>:nstep_init,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:parameter_scan},
     :nstep_inc=>
      {:should_include=>"true",
       :description=>"The parameter will be changed every nstep_inc.",
       :help=>
        "When the increment condition is 'n_timesteps', the parameter will be changed every nstep_inc.",
       :code_name=>:nstep_inc,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:parameter_scan},
     :delta_t_init=>
      {:should_include=>"true",
       :description=>
        "The parameter will not be changed until delta_t_init time has elapsed from the beginning of the simulation.",
       :help=>
        "When the increment condition is 'delta_t', the parameter will not be changed until delta_t_init time has elapsed from the beginning of the simulation. Note, that if the simulation is restarted, this parameter will measure from beginning of original simulation.",
       :code_name=>:delta_t_init,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:parameter_scan},
     :delta_t_inc=>
      {:should_include=>"true",
       :description=>
        "The parameter will be changed every time delta_t time has elapsed.",
       :help=>
        "When the increment condition is 'delta_t', the parameter will be changed every time delta_t time has elapsed.",
       :code_name=>:delta_t_inc,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:parameter_scan},
     :scan_spec=>
      {:should_include=>"true",
       :description=>
        "When parameter pertains to a species, specifies the index of the species.",
       :help=>
        "When parameter pertains to a species, specifies the index of the species.",
       :code_name=>:scan_spec,
       :must_pass=>
        [{:test=>"kind_of? Integer",
          :explanation=>"This variable must be an integer."}],
       :type=>:Integer,
       :module=>:parameter_scan},
     :scan_restarted=>
      {:should_include=>"true",
       :description=>
        "If true current value of the scan parameter will be read from the restart files.",
       :help=>
        "Must be set to true if the current value of the scan parameter must be read from the restart files.  Otherwise, the scan will start from the beginning.",
       :code_name=>:scan_restarted,
       :must_pass=>
        [{:test=>"kind_of? String and FORTRAN_BOOLS.include? self",
          :explanation=>
           "This variable must be a fortran boolean. (In Ruby this is represented as a string: e.g. '.true.')"}],
       :type=>:Fortran_Bool,
       :module=>:parameter_scan},
     :target_val=>
      {:should_include=>"true",
       :description=>"Specifies the value to be targeted.",
       :help=>
        "If the scan is being run in 'target' or 'root_finding'  mode, specifies the value to be targeted. The scan will complete when this target value is reached.",
       :code_name=>:target_val,
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}],
       :type=>:Float,
       :module=>:parameter_scan}}}}
