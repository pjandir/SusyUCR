


   void setup_ws_files( float min_signal_frac = 0.00,
                        bool skip_testfit = true,
                        bool skip_modelconfig = true,
                        float saveall_below_N = 0.,
                        bool no_rounding = true) {

      gROOT -> LoadMacro("wrapper_build_ra2b_ws2.c") ;

      wrapper_build_ra2b_ws2( "outputfiles/lhbuilder-input-t1bbbbH.txt", "outputfiles/ws-t1bbbbH.root",
             min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;

      wrapper_build_ra2b_ws2( "outputfiles/lhbuilder-input-t1bbbbC.txt", "outputfiles/ws-t1bbbbC.root",
             min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;


      wrapper_build_ra2b_ws2( "outputfiles/lhbuilder-input-t1ttttH.txt", "outputfiles/ws-t1ttttH.root",
             min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;

      wrapper_build_ra2b_ws2( "outputfiles/lhbuilder-input-t1ttttC.txt", "outputfiles/ws-t1ttttC.root",
             min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;


      wrapper_build_ra2b_ws2( "outputfiles/lhbuilder-input-t1qqqqH.txt", "outputfiles/ws-t1qqqqH.root",
             min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;

      wrapper_build_ra2b_ws2( "outputfiles/lhbuilder-input-t1qqqqC.txt", "outputfiles/ws-t1qqqqC.root",
             min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;


   } // setup_ws_files




