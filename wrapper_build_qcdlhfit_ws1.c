



void wrapper_build_qcdlhfit_ws1(
                            const char* outfile = "outputfiles/ws-kqcd-lhfit-perfect-qcd-closure-random-nobs-with-constraints.root",
                            const char* fname_fitconfig = "kqcd-fitconfig.txt",
                            const char* fname_data    = "outputfiles/kqcd-input-fakedata-perfect-qcd-closure-random-nobs.txt",
                            bool  skip_testfit = false,
                            bool  skip_modelconfig = true,
                            const char* fname_lostlep = "outputfiles/kqcd-input-lostlep.txt",
                            const char* fname_hadtau  = "outputfiles/kqcd-input-hadtau.txt",
                            const char* fname_znunu   = "outputfiles/kqcd-input-znunu.txt",
                            const char* fname_sigmc   = "outputfiles/kqcd-input-sigmc-t1bbbbH.txt"
           ) {

   gROOT->LoadMacro("RooPoissonLogEval.cxx+") ;
   gROOT->LoadMacro("RooProdPdfLogSum.cxx+") ;
   gROOT->LoadMacro("build_qcdlhfit_ws1.c+") ;

   build_qcdlhfit_ws1( outfile, fname_fitconfig, fname_data, skip_testfit, skip_modelconfig, fname_lostlep, fname_hadtau, fname_znunu, fname_sigmc ) ;

}

