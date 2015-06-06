



void wrapper_build_ra2b_ws4b(
          const char* infile = "outputfiles/lhbuilder-input-v4b-t1bbbbH.txt",
          const char* outfile = "outputfiles/ws-v4b-t1bbbbH.root",
          float min_signal_frac = 0.0,
          bool skip_testfit = true,
          bool skip_modelconfig = true,
          float saveall_below_N = 0.,
          bool no_rounding = true,
          bool use_mht_ratios = false ) {

   gROOT->LoadMacro("RooPoissonLogEval.cxx+") ;
   gROOT->LoadMacro("RooProdPdfLogSum.cxx+") ;
   gROOT->LoadMacro("build_ra2b_ws4b.c+") ;

   build_ra2b_ws4b( infile, outfile, min_signal_frac, skip_testfit, skip_modelconfig, saveall_below_N, no_rounding, use_mht_ratios ) ;

}

