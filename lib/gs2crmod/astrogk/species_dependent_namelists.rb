{:species_parameters=>
  {:description=>"SPECIES ",
   :include_conditions=>"true",
   :variables=>
    {:z=>
      {:help=>"Charge",
       :include_conditions=>"true",
       :description=>"Charge",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :mass=>
      {:help=>"Mass",
       :include_conditions=>"true",
       :description=>"Mass",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :dens=>
      {:help=>"Density\t",
       :include_conditions=>"true",
       :description=>"Density\t",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :temp=>
      {:help=>"Temperature",
       :include_conditions=>"true",
       :description=>"Temperature",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :tprim=>
      {:help=>"6.0   !-1/T (dT/drho)",
       :include_conditions=>"true",
       :description=>"6.0   !-1/T (dT/drho)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :fprim=>
      {:help=>"2.22   !-1/n (dn/drho)",
       :include_conditions=>"true",
       :description=>"2.22   !-1/n (dn/drho)",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :uprim=>
      {:help=>"?",
       :include_conditions=>"true",
       :description=>"?",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :vnewk=>
      {:help=>"1.e-2   !collisionality parameter",
       :include_conditions=>"true",
       :description=>"1.e-2   !collisionality parameter",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :type=>
      {:help=>"'ion' Thermal ion species ",
       :include_conditions=>"true",
       :description=>"'ion' Thermal ion species ",
       :tests=>["Tst::STRING"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? String",
          :explanation=>"This variable must be a string."}]},
     :dens0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:dens0,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[1.0]},
     :u0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:u0,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[1.0]},
     :uprim2=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:uprim2,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[0.0]},
     :nustar=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:nustar,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[-1.0]},
     :nu=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:nu,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[-1.0]},
     :nu_h=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:nu_h,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[0.0]},
     :tperp0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tperp0,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[0.0]},
     :tpar0=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::FLOAT"],
       :gs2_name=>:tpar0,
       :must_pass=>
        {:test=>"kind_of? Numeric",
         :explanation=>
          "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."},
       :autoscanned_defaults=>[0.0]}}},
 :dist_fn_species_knobs=>
  {:description=>"",
   :include_conditions=>"true",
   :variables=>
    {:fexpr=>
      {:help=>
        "0.5   !Temporal implicitness parameter. Recommended value: 0.48",
       :include_conditions=>"true",
       :description=>
        "0.5   !Temporal implicitness parameter. Recommended value: 0.48",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :fexpi=>
      {:help=>nil,
       :include_conditions=>"true",
       :description=>nil,
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :bakdif=>
      {:help=>"0.0   !Spatial implicitness parameter. Recommended value: 0.05",
       :include_conditions=>"true",
       :description=>
        "0.0   !Spatial implicitness parameter. Recommended value: 0.05",
       :tests=>["Tst::FLOAT"],
       :autoscanned_defaults=>[],
       :must_pass=>
        [{:test=>"kind_of? Numeric",
          :explanation=>
           "This variable must be a floating point number (an integer is also acceptable: it will be converted into a floating point number)."}]},
     :bd_exp=>
      {:should_include=>"true",
       :description=>nil,
       :help=>nil,
       :tests=>["Tst::INT"],
       :gs2_name=>:bd_exp,
       :must_pass=>
        {:test=>"kind_of? Integer",
         :explanation=>"This variable must be an integer."},
       :autoscanned_defaults=>["bd_exp_out"]}}}}
