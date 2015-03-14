



void wrapper_build_ra2b_ws2(
          const char* infile = "outputfiles/lhbuilder-input-t1bbbbH.txt",
          const char* outfile = "outputfiles/ws-t1bbbbH.root",
          float min_signal_frac = 0.0,
          bool skip_testfit = true,
          bool skip_modelconfig = true,
          float saveall_below_N = 0.,
          bool no_rounding = true ) {

   gROOT->LoadMacro("RooPoissonLogEval.cxx+") ;
   gROOT->LoadMacro("RooProdPdfLogSum.cxx+") ;
   gROOT->LoadMacro("build_ra2b_ws2.c+") ;

   build_ra2b_ws2( infile, outfile, min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding ) ;

}

