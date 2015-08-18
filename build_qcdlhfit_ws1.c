
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooUniform.h"
#include "RooStats/ModelConfig.h"

#include "RooProdPdfLogSum.h"
#include "RooPoissonLogEval.h"

#include <fstream>

   using namespace RooFit ;
   using namespace RooStats ;

   RooArgSet* globalObservables ;
   RooArgSet* allNuisances ;
   RooArgSet* allNuisancePdfs ;
   RooArgSet* pdf_sbIndexList ;
   RooArgSet* allBGmuPars ;

   bool  find_line( ifstream& ifs, const char* key ) ;
   float find_line_val( ifstream& ifs, const char* key ) ;
   bool  find_line_val_err( ifstream& ifs, const char* key, float& val, float& err ) ;
   float get_par_max( float mu ) ;

   RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err, int sbi = -1 ) ;
   RooAbsReal* makeCorrelatedLognormalConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false ) ;

   RooWorkspace* ws_pointer ;

  //=================================================================================

   void build_qcdlhfit_ws1(
                            const char* outfile = "outputfiles/ws-kqcd-lhfit-test.root",
                            const char* fname_fitconfig = "outputfiles/kqcd-fitconfig-with-constraints1.txt",
                            const char* fname_data    = "outputfiles/kqcd-input-fakedata.txt",
                            bool  skip_testfit = false,
                            bool  skip_modelconfig = true,
                            const char* fname_lostlep = "outputfiles/kqcd-input-lostlep.txt",
                            const char* fname_hadtau  = "outputfiles/kqcd-input-hadtau.txt",
                            const char* fname_znunu   = "outputfiles/kqcd-input-znunu.txt",
                            const char* fname_sigmc   = "outputfiles/kqcd-input-sigmc-t1bbbbH.txt"
                          ) {


      bool no_rounding(true) ;

      bool correlate_nonqcd_errors(true) ;
      //bool correlate_nonqcd_errors(false) ;

      char pname[100] ;
      char pname2[100] ;

      ifstream ifs ;
      ifs.open( fname_fitconfig ) ;
      if ( !ifs.good() ) {
         printf("\n\n *** Problem with fitconfig file : %s\n\n", fname_fitconfig ) ; return ;
      }

      printf("\n QCD model parameters\n") ;
      int n_qcd_kht_pars = find_line_val( ifs, "N-QCD-Kht-pars" ) ;
      float qcd_kht_val[10] ;
      float qcd_kht_err[10] ;
      for ( int i=1; i<=n_qcd_kht_pars; i++ ) {
         sprintf( pname, "QCD-Kht%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_kht_val[i], qcd_kht_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_kht_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_kht_val[i] ) ;
         } else if ( qcd_kht_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_kht_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_kht_val[i], qcd_kht_err[i] ) ;
         }
      } // i

      int n_qcd_kmht_pars = find_line_val( ifs, "N-QCD-Kmht-pars" ) ;
      float qcd_kmht_val[10] ;
      float qcd_kmht_err[10] ;
      for ( int i=1; i<=n_qcd_kmht_pars; i++ ) {
         sprintf( pname, "QCD-Kmht%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_kmht_val[i], qcd_kmht_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_kmht_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_kmht_val[i] ) ;
         } else if ( qcd_kmht_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_kmht_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_kmht_val[i], qcd_kmht_err[i] ) ;
         }
      } // i

      int n_qcd_knjet_pars = find_line_val( ifs, "N-QCD-Knjet-pars" ) ;
      float qcd_knjet_val[10] ;
      float qcd_knjet_err[10] ;
      for ( int i=1; i<=n_qcd_knjet_pars; i++ ) {
         sprintf( pname, "QCD-Knjet%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_knjet_val[i], qcd_knjet_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_knjet_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_knjet_val[i] ) ;
         } else if ( qcd_knjet_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_knjet_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_knjet_val[i], qcd_knjet_err[i] ) ;
         }
      } // i


      int n_qcd_knb_pars = find_line_val( ifs, "N-QCD-Knb-pars" ) ;
      float qcd_knb_val[10] ;
      float qcd_knb_err[10] ;
      for ( int i=1; i<=n_qcd_knb_pars; i++ ) {
         sprintf( pname, "QCD-Knb%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_knb_val[i], qcd_knb_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_knb_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_knb_val[i] ) ;
         } else if ( qcd_knb_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_knb_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_knb_val[i], qcd_knb_err[i] ) ;
         }
      } // i




    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      printf("\n\n Creating workspace.\n\n") ;

      RooWorkspace workspace("ws") ;
      ws_pointer = &workspace ;

      workspace.autoImportClassCode(true) ;

      globalObservables      = new RooArgSet("globalObservables");
      allNuisances           = new RooArgSet("allNuisances");
      allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      pdf_sbIndexList        = new RooArgSet("pdf_sbIndexList") ;
      allBGmuPars            = new RooArgSet("allBGmuPars") ;
      RooArgSet* observedParametersList = new RooArgSet("observables") ;

      RooArgSet pdflist ;

      sprintf( pname, "sig_strength" ) ;
      RooRealVar* rv_sig_strength = new RooRealVar( pname, pname, 0.0, 0., 10. ) ; // nominal value is zero
      rv_sig_strength -> setConstant(kTRUE) ; // fixed
      rv_sig_strength -> Print() ;
      printf("  %s\n\n", pname ) ;


      RooArgSet sbIndexList("sbIndexList") ;






      printf("\n  QCD Model parameters:\n") ;

      RooAbsReal* rv_qcd_kht[10] ;
      RooAbsReal* rv_qcd_kmht[10] ;
      RooAbsReal* rv_qcd_knjet[10] ;
      RooAbsReal* rv_qcd_knb[10] ;

      for ( int i=1; i<=n_qcd_kht_pars; i++ ) {
         sprintf( pname, "Kqcd_ht%d", i ) ;
         if ( qcd_kht_err[i] < 0 ) {
            rv_qcd_kht[i] = new RooRealVar( pname, pname, qcd_kht_val[i], 0., 0.50 ) ;
         } else {
            rv_qcd_kht[i] = makeLognormalConstraint( pname, qcd_kht_val[i], qcd_kht_err[i] ) ;
         }
      } // i

      for ( int i=1; i<=n_qcd_kmht_pars; i++ ) {
         sprintf( pname, "Kqcd_mht%d", i ) ;
         if ( qcd_kmht_err[i] < 0 ) {
            rv_qcd_kmht[i] = new RooRealVar( pname, pname, qcd_kmht_val[i], 0.0, 6.0 ) ;
         } else {
            rv_qcd_kmht[i] = makeLognormalConstraint( pname, qcd_kmht_val[i], qcd_kmht_err[i] ) ;
         }
      } // i

      for ( int i=1; i<=n_qcd_knjet_pars; i++ ) {
         sprintf( pname, "Kqcd_njet%d", i ) ;
         if ( qcd_knjet_err[i] < 0 ) {
            rv_qcd_knjet[i] = new RooRealVar( pname, pname, qcd_knjet_val[i], 0.0, 6.0 ) ;
         } else {
            rv_qcd_knjet[i] = makeLognormalConstraint( pname, qcd_knjet_val[i], qcd_knjet_err[i] ) ;
         }
      } // i

      for ( int i=1; i<=n_qcd_knb_pars; i++ ) {
         sprintf( pname, "Kqcd_nb%d", i ) ;
         if ( qcd_knb_err[i] < 0 ) {
            rv_qcd_knb[i] = new RooRealVar( pname, pname, qcd_knb_val[i], 0., 3.0 ) ;
         } else {
            rv_qcd_knb[i] = makeLognormalConstraint( pname, qcd_knb_val[i], qcd_knb_err[i] ) ;
         }
      } // i







      ifstream ifs_data ;
      ifstream ifs_lostlep ;
      ifstream ifs_hadtau ;
      ifstream ifs_znunu ;
      ifstream ifs_sigmc ;


      ifs_data.open( fname_data ) ;
      if ( !ifs_data.good() ) { printf("\n\n *** Bad input data file: %s\n\n", fname_data ) ; return ; }

      ifs_lostlep.open( fname_lostlep ) ;
      if ( !ifs_lostlep.good() ) { printf("\n\n *** Bad input lostlep file: %s\n\n", fname_lostlep ) ; return ; }

      ifs_hadtau.open( fname_hadtau ) ;
      if ( !ifs_hadtau.good() ) { printf("\n\n *** Bad input hadtau file: %s\n\n", fname_hadtau ) ; return ; }

      ifs_znunu.open( fname_znunu ) ;
      if ( !ifs_znunu.good() ) { printf("\n\n *** Bad input znunu file: %s\n\n", fname_znunu ) ; return ; }

      ifs_sigmc.open( fname_sigmc ) ;
      if ( !ifs_sigmc.good() ) { printf("\n\n *** Bad input sigmc file: %s\n\n", fname_sigmc ) ; return ; }

      printf("\n\n Reading input event counts.\n\n" ) ;

      while ( ifs_data.good() ) {

         printf("\n") ;

         TString ts ;

       //-- data
         ts.ReadLine( ifs_data ) ;
         if ( !ifs_data.good() ) { printf("\n\n  Reached end of data input file.\n\n" ) ; break ; }
         int bin_index ;
         char bin_name[100] ;
         int  data_ldp_val(0), data_zl_val(0) ;
         sscanf( ts.Data(), "%d %s %d %d", &bin_index, bin_name, &data_ldp_val, &data_zl_val ) ;
         int fb_nji(-1), fb_mbi(-1), fb_hbi(-1) ;
         sscanf( bin_name, "FB-Njet%d-Nbsum-MHT%d-HT%d", &fb_nji, &fb_mbi, &fb_hbi ) ;
         int fb_nbi(1) ;

         printf( " bin %2d : %30s,  Njet index %d,  MHT index %d,  HT index %d\n", bin_index, bin_name, fb_nji, fb_mbi, fb_hbi ) ;
         printf( "      data    : LDP %5d            ,   ZL %4d\n", data_ldp_val, data_zl_val ) ;





         sprintf( pname, "R_qcd_ldp_%s", bin_name ) ;
         RooFormulaVar* rv_R_qcd_ldp = new RooFormulaVar( pname, "@0 * @1 * @2 * @3", 
            RooArgSet(
               *(rv_qcd_kht[ fb_hbi ]),
               *(rv_qcd_kmht[ fb_mbi ]),
               *(rv_qcd_knjet[ fb_nji ]),
               *(rv_qcd_knb[ fb_nbi ])
                     ) ) ;
         rv_R_qcd_ldp -> Print() ;




         sprintf( pname, "Nldp_%s", bin_name ) ;
         RooRealVar* rv_Nldp = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rv_Nldp -> setVal( data_ldp_val ) ;
         rv_Nldp -> setConstant( kTRUE ) ;
         rv_Nldp -> Print() ;
         observedParametersList -> add( *rv_Nldp ) ;

         sprintf( pname, "Nzl_%s", bin_name ) ;
         RooRealVar* rv_Nzl = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rv_Nzl -> setVal( data_zl_val ) ;
         rv_Nzl -> setConstant( kTRUE ) ;
         rv_Nzl -> Print() ;
         observedParametersList -> add( *rv_Nzl ) ;




       //-- lostlep
         ts.ReadLine( ifs_lostlep ) ;
         int lostlep_bin_index ;
         char lostlep_bin_name[100] ;
         float lostlep_ldp_val, lostlep_ldp_err, lostlep_zl_val, lostlep_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &lostlep_bin_index, lostlep_bin_name, &lostlep_ldp_val, &lostlep_ldp_err, &lostlep_zl_val, &lostlep_zl_err ) ;
         if ( lostlep_bin_index != bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", bin_index, lostlep_bin_index ) ; return ; }
         if ( strcmp( bin_name, lostlep_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, lostlep_bin_name ) ; return ; }
         printf("      lostlep : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", lostlep_ldp_val, lostlep_ldp_err, lostlep_zl_val, lostlep_zl_err ) ;
         RooAbsReal* rv_mu_lostlep_ldp(0x0) ;
         RooAbsReal* rv_mu_lostlep_zl(0x0) ;
         if ( correlate_nonqcd_errors ) {
            sprintf( pname, "mu_lostlep_ldp_%s", bin_name ) ;
            rv_mu_lostlep_ldp = makeCorrelatedLognormalConstraint( pname, lostlep_ldp_val, lostlep_ldp_err, "lostlep_ldp" ) ;
            sprintf( pname, "mu_lostlep_zl_%s", bin_name ) ;
            rv_mu_lostlep_zl  = makeCorrelatedLognormalConstraint( pname, lostlep_zl_val, lostlep_zl_err, "lostlep_zl" ) ;
         } else {
            sprintf( pname, "mu_lostlep_ldp_%s", bin_name ) ;
            rv_mu_lostlep_ldp = makeLognormalConstraint( pname, lostlep_ldp_val, lostlep_ldp_err ) ;
            sprintf( pname, "mu_lostlep_zl_%s", bin_name ) ;
            rv_mu_lostlep_zl  = makeLognormalConstraint( pname, lostlep_zl_val, lostlep_zl_err ) ;
         }

       //-- hadtau
         ts.ReadLine( ifs_hadtau ) ;
         int hadtau_bin_index ;
         char hadtau_bin_name[100] ;
         float hadtau_ldp_val, hadtau_ldp_err, hadtau_zl_val, hadtau_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &hadtau_bin_index, hadtau_bin_name, &hadtau_ldp_val, &hadtau_ldp_err, &hadtau_zl_val, &hadtau_zl_err ) ;
         if ( hadtau_bin_index != bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", bin_index, hadtau_bin_index ) ; return ; }
         if ( strcmp( bin_name, hadtau_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, hadtau_bin_name ) ; return ; }
         printf("      hadtau  : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", hadtau_ldp_val, hadtau_ldp_err, hadtau_zl_val, hadtau_zl_err ) ;
         RooAbsReal* rv_mu_hadtau_ldp(0x0) ;
         RooAbsReal* rv_mu_hadtau_zl(0x0) ;
         if ( correlate_nonqcd_errors ) {
            sprintf( pname, "mu_hadtau_ldp_%s", bin_name ) ;
            rv_mu_hadtau_ldp = makeCorrelatedLognormalConstraint( pname, hadtau_ldp_val, hadtau_ldp_err, "hadtau_ldp" ) ;
            sprintf( pname, "mu_hadtau_zl_%s", bin_name ) ;
            rv_mu_hadtau_zl  = makeCorrelatedLognormalConstraint( pname, hadtau_zl_val, hadtau_zl_err, "hadtau_zl" ) ;
         } else {
            sprintf( pname, "mu_hadtau_ldp_%s", bin_name ) ;
            rv_mu_hadtau_ldp = makeLognormalConstraint( pname, hadtau_ldp_val, hadtau_ldp_err ) ;
            sprintf( pname, "mu_hadtau_zl_%s", bin_name ) ;
            rv_mu_hadtau_zl  = makeLognormalConstraint( pname, hadtau_zl_val, hadtau_zl_err ) ;
         }

       //-- znunu
         ts.ReadLine( ifs_znunu ) ;
         int znunu_bin_index ;
         char znunu_bin_name[100] ;
         float znunu_ldp_val, znunu_ldp_err, znunu_zl_val, znunu_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &znunu_bin_index, znunu_bin_name, &znunu_ldp_val, &znunu_ldp_err, &znunu_zl_val, &znunu_zl_err ) ;
         if ( znunu_bin_index != bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", bin_index, znunu_bin_index ) ; return ; }
         if ( strcmp( bin_name, znunu_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, znunu_bin_name ) ; return ; }
         printf("      znunu   : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", znunu_ldp_val, znunu_ldp_err, znunu_zl_val, znunu_zl_err ) ;
         RooAbsReal* rv_mu_znunu_ldp(0x0) ;
         RooAbsReal* rv_mu_znunu_zl(0x0) ;
         if ( correlate_nonqcd_errors ) {
            sprintf( pname, "mu_znunu_ldp_%s", bin_name ) ;
            rv_mu_znunu_ldp = makeCorrelatedLognormalConstraint( pname, znunu_ldp_val, znunu_ldp_err, "znunu_ldp" ) ;
            sprintf( pname, "mu_znunu_zl_%s", bin_name ) ;
            rv_mu_znunu_zl  = makeCorrelatedLognormalConstraint( pname, znunu_zl_val, znunu_zl_err, "znunu_zl" ) ;
         } else {
            sprintf( pname, "mu_znunu_ldp_%s", bin_name ) ;
            rv_mu_znunu_ldp = makeLognormalConstraint( pname, znunu_ldp_val, znunu_ldp_err ) ;
            sprintf( pname, "mu_znunu_zl_%s", bin_name ) ;
            rv_mu_znunu_zl  = makeLognormalConstraint( pname, znunu_zl_val, znunu_zl_err ) ;
         }



       //-- sigmc
         ts.ReadLine( ifs_sigmc ) ;
         int sigmc_bin_index ;
         char sigmc_bin_name[100] ;
         float sigmc_ldp_val, sigmc_ldp_err, sigmc_zl_val, sigmc_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &sigmc_bin_index, sigmc_bin_name, &sigmc_ldp_val, &sigmc_ldp_err, &sigmc_zl_val, &sigmc_zl_err ) ;
         if ( sigmc_bin_index != bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", bin_index, sigmc_bin_index ) ; return ; }
         if ( strcmp( bin_name, sigmc_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", bin_name, sigmc_bin_name ) ; return ; }
         printf("      sigmc   : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", sigmc_ldp_val, sigmc_ldp_err, sigmc_zl_val, sigmc_zl_err ) ;

         sprintf( pname, "mu_sig0_ldp_%s", bin_name ) ;
         RooRealVar* rv_mu_sig0_ldp = new RooRealVar( pname, pname, sigmc_ldp_val, 0., 1.e6 ) ;
         rv_mu_sig0_ldp -> setConstant( kTRUE ) ;
         rv_mu_sig0_ldp -> Print() ;
         sprintf( pname, "mu_sig_ldp_%s", bin_name ) ;
         RooFormulaVar* rv_mu_sig_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_ldp ) ) ;
         rv_mu_sig_ldp -> Print() ;

         sprintf( pname, "mu_sig0_zl_%s", bin_name ) ;
         RooRealVar* rv_mu_sig0_zl = new RooRealVar( pname, pname, sigmc_zl_val, 0., 1.e6 ) ;
         rv_mu_sig0_zl -> setConstant( kTRUE ) ;
         rv_mu_sig0_zl -> Print() ;
         sprintf( pname, "mu_sig_zl_%s", bin_name ) ;
         RooFormulaVar* rv_mu_sig_zl = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_zl ) ) ;
         rv_mu_sig_zl -> Print() ;











         double initial_qcd_ldp_val =  data_ldp_val - ( lostlep_ldp_val + hadtau_ldp_val + znunu_ldp_val ) ;
         if ( initial_qcd_ldp_val < 0. ) initial_qcd_ldp_val = 0. ;

         sprintf( pname, "mu_qcd_ldp_%s", bin_name ) ;
         RooRealVar* rv_mu_qcd_ldp = new RooRealVar( pname, pname, initial_qcd_ldp_val, 0., get_par_max( data_ldp_val ) ) ;
         rv_mu_qcd_ldp -> setConstant( kFALSE ) ;
         rv_mu_qcd_ldp -> Print() ;
         allBGmuPars -> add( *rv_mu_qcd_ldp ) ;


         sprintf( pname, "mu_qcd_zl_%s", bin_name ) ;
         RooFormulaVar* rv_mu_qcd_zl = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_R_qcd_ldp, *rv_mu_qcd_ldp ) ) ;
         rv_mu_qcd_zl -> Print() ;





         sprintf( pname, "n_ldp_%s", bin_name ) ;
         RooFormulaVar* rv_n_ldp = new RooFormulaVar( pname, "@0 + @1 + @2 + @3 + @4", RooArgSet(
             *rv_mu_qcd_ldp, *rv_mu_lostlep_ldp, *rv_mu_hadtau_ldp, *rv_mu_znunu_ldp, *rv_mu_sig_ldp ) ) ;
         rv_n_ldp -> Print() ;

         sprintf( pname, "n_zl_%s", bin_name ) ;
         RooFormulaVar* rv_n_zl = new RooFormulaVar( pname, "@0 + @1 + @2 + @3 + @4", RooArgSet(
             *rv_mu_qcd_zl, *rv_mu_lostlep_zl, *rv_mu_hadtau_zl, *rv_mu_znunu_zl, *rv_mu_sig_zl ) ) ;
         rv_n_zl -> Print() ;





         sprintf( pname, "pdf_ldp_%s", bin_name ) ;
         RooPoissonLogEval* rv_pdf_ldp = new RooPoissonLogEval( pname, pname, *rv_Nldp, *rv_n_ldp, no_rounding ) ;
         rv_pdf_ldp -> Print() ;

         pdflist.add( *rv_pdf_ldp ) ;


         sprintf( pname, "pdf_zl_%s", bin_name ) ;
         RooPoissonLogEval* rv_pdf_zl = new RooPoissonLogEval( pname, pname, *rv_Nzl, *rv_n_zl, no_rounding ) ;
         rv_pdf_zl -> Print() ;

         pdflist.add( *rv_pdf_zl ) ;



      } // reading bins from data file ( ifs_data )



      printf("\n\n Creating and importing dataset into workspace.\n\n") ;

      RooDataSet* dsObserved = new RooDataSet("observed_rds", "observed data values", *observedParametersList ) ;
      dsObserved -> add( *observedParametersList ) ;
      workspace.import( *dsObserved ) ;



      pdflist.add( *allNuisancePdfs ) ;

      printf("\n List of all PDFs\n") ;
      pdflist.Print() ;
      printf("\n") ;

      RooProdPdfLogSum* likelihood = new RooProdPdfLogSum( "likelihood", "likelihood", pdflist ) ;
      likelihood->Print() ;







      if ( !skip_testfit ) {

         printf("\n\n Running a test fit.\n\n") ;

         printf("\n\n =============================================\n\n") ;
         ////likelihood -> fitTo( *dsObserved, PrintLevel(3), Hesse(0), Minos(0) ) ;
         ///likelihood -> fitTo( *dsObserved, Optimize(0), PrintLevel(3), Hesse(0), Minos(0) ) ;
         likelihood -> fitTo( *dsObserved, Optimize(0), PrintLevel(3), Hesse(0), Minos(0) ) ;
         printf("\n\n =============================================\n\n") ;

      }


      if ( !skip_modelconfig ) {

        //-- Set up RooStats models.

         printf("\n\n Setting up S+B model.\n\n") ;

         RooArgSet poi( *rv_sig_strength, "poi" ) ;
         RooUniform signal_prior( "signal_prior", "signal_prior", *rv_sig_strength ) ;

         ModelConfig sbModel ("SbModel");
         sbModel.SetWorkspace( workspace ) ;
         sbModel.SetPdf( *likelihood ) ;
         sbModel.SetParametersOfInterest( poi );
         sbModel.SetPriorPdf(signal_prior);
         sbModel.SetObservables( *observedParametersList );
         sbModel.SetNuisanceParameters( *allNuisances );
         sbModel.SetGlobalObservables( *globalObservables );

         workspace.Print() ;

         printf("\n\n Doing fit for S+B model.\n" ) ; fflush(stdout) ;

         RooAbsReal* pNll = sbModel.GetPdf()->createNLL(*dsObserved);
         RooAbsReal* pProfile = pNll->createProfile(RooArgSet());
         pProfile->getVal();
         RooArgSet* pPoiAndNuisance = new RooArgSet();
         pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
         if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
         printf("\n\n Will save these parameter points that correspond to the fit to data.\n\n") ; fflush(stdout) ;
         pPoiAndNuisance->Print("v");
         sbModel.SetSnapshot(*pPoiAndNuisance);
         workspace.import (sbModel);

         delete pProfile ;
         delete pNll ;
         delete pPoiAndNuisance ;

         printf("\n\n Setting up BG-only model.\n\n") ;

         ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
         bModel.SetName("BModel");
         bModel.SetWorkspace(workspace);

         printf("\n\n Doing fit for BG-only model.\n" ) ; fflush(stdout) ;
         pNll = bModel.GetPdf()->createNLL(*dsObserved);
         pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
         ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.);
         pProfile->getVal();
         pPoiAndNuisance = new RooArgSet();
         pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
         if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
         printf("\n\n Should use these parameter points to generate pseudo data for bkg only.\n\n") ; fflush(stdout) ;
         pPoiAndNuisance->Print("v");
         bModel.SetSnapshot(*pPoiAndNuisance);
         workspace.import (bModel);

         delete pProfile ;
         delete pNll ;
         delete pPoiAndNuisance ;

      } else {

         workspace.import( *likelihood ) ;

      }

      workspace.defineSet( "all_nuisance_pars", *allNuisances, kFALSE ) ;  // for convenience. do not import if missing.  should already be in there.
      workspace.defineSet( "all_nuisance_pdfs", *allNuisancePdfs, kFALSE ) ;  // for convenience. do not import if missing.  should already be in there.
      workspace.defineSet( "all_bg_mu_pars", *allBGmuPars, kFALSE ) ;  // for convenience. do not import if missing.  should already be in there.

      workspace.Print() ;

      printf("\n\n Saving workspace in : %s\n\n", outfile ) ;

      gSystem->Exec(" mkdir -p outputfiles " ) ;

      workspace.writeToFile( outfile ) ;














   } // build_qcdlhfit_ws1

  //=================================================================================

   bool find_line( ifstream& ifs, const char* key ) {

      ifs.seekg(0) ;
      TString ts ;

      int lc(0) ;
      while ( ifs.good() ) {
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         /////////// printf(" %4d : %s\n", lc++, ts.Data() ) ;
         if ( strcmp( first->GetString().Data(), key ) == 0 ) {
            /////////// printf(" Found %s\n", key ) ;
            return true ;
         }
      }

      return false ;


   } // find_line

  //=================================================================================

   float find_line_val( ifstream& ifs, const char* key ) {

      ifs.seekg(0) ;
      TString ts ;

      int lc(0) ;
      while ( ifs.good() ) {
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         /////////// printf(" %4d : %s\n", lc++, ts.Data() ) ;
         if ( strcmp( first->GetString().Data(), key ) == 0 ) {
            /////////// printf(" Found %s\n", key ) ;
            float val(-1) ;
            if ( tokens->GetEntries() > 1 ) {
               TObjString* second = (TObjString*) tokens->At(1) ;
               sscanf( second->GetString().Data(), "%f", &val ) ;
               return val ;
            } else {
               return -2 ;
            }
         }
      }

      return -1 ;

   } // find_line_val

  //=================================================================================

   bool find_line_val_err( ifstream& ifs, const char* key, float& val, float& err ) {

      val = 0. ;
      err = 0. ;

      ifs.seekg(0) ;
      TString ts ;

      int lc(0) ;
      while ( ifs.good() ) {
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         /////////// printf(" %4d : %s\n", lc++, ts.Data() ) ;
         if ( strcmp( first->GetString().Data(), key ) == 0 ) {
            /////////// printf(" Found %s\n", key ) ;
            if ( tokens->GetEntries() >= 3 ) {
               TObjString* second = (TObjString*) tokens->At(1) ;
               sscanf( second->GetString().Data(), "%f", &val ) ;
               TObjString* third = (TObjString*) tokens->At(2) ;
               sscanf( third->GetString().Data(), "%f", &err ) ;
               return true ;
            }
         }
      }

      return false ;


   } // find_line

  //=================================================================================

    float get_par_max( float mu ) {

   /// if ( mu <= 1. ) return 10. ;
   /// if ( mu <= 2. ) return 15. ;
   /// if ( mu <= 50. ) return mu + 4*sqrt(mu) ;
   /// if ( mu <= 100. ) return mu + 3*sqrt(mu) ;
   /// return mu + 2*sqrt(mu) ;

       if ( mu <= 1. ) return 10. ;
       if ( mu <= 2. ) return 15. ;
       if ( mu <= 50. ) return mu + 6*sqrt(mu) ;
       if ( mu <= 100. ) return mu + 5*sqrt(mu) ;
       return mu + 4*sqrt(mu) ;


    } // get_par_max

  //=================================================================================


    RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err, int sbi ) {

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       if ( NP_err <= 0. ) {
          printf(" makeLognormalConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          ws_pointer -> import( *g_mean ) ;
          ws_pointer -> import( *g_sigma ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       char pname[1000] ;
       sprintf( pname, "prim_%s", NP_name ) ;

       printf(" makeLognormalConstraint : creating primary log-normal variable %s\n", pname ) ;
       RooRealVar* np_prim_rrv = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_rrv -> setVal( 0. ) ;
       np_prim_rrv -> setConstant( kFALSE ) ;

       sprintf( pname, "prim_mean_%s", NP_name ) ;
       RooRealVar* np_prim_mean = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_mean->setConstant(kTRUE) ;

       sprintf( pname, "prim_sigma_%s", NP_name ) ;
       RooConstVar* np_prim_sigma = new RooConstVar( pname, pname, 1. ) ;


       char pdfname[1000] ;
       sprintf( pdfname, "pdf_prim_%s", NP_name ) ;
       RooGaussian* np_prim_pdf = new RooGaussian( pdfname, pdfname, *np_prim_rrv, *np_prim_mean, *np_prim_sigma ) ;

       allNuisances -> add( *np_prim_rrv ) ;
       allNuisancePdfs -> add( *np_prim_pdf ) ;
       globalObservables -> add( *np_prim_mean ) ;

       sprintf( pname, "%s_sb_index", np_prim_pdf->GetName() ) ;
       RooConstVar* rv_pdf_sb_index = new RooConstVar( pname, pname, sbi ) ;
       pdf_sbIndexList -> add( *rv_pdf_sb_index ) ;



       //-- compute the log-normal-distributed parameter from the primary parameter.

       //--- This is the new way.  RMS of lognormal is much closer to sigma when sigma is
       //    large, doing it this way.  When sigma/mean is small, they are about the same.
       //    That is, exp(sigma/mean) is close to (sigma/mean + 1).  This one is better when
       //    sigma/mean is not small.  The high-side tail is not as strong.
       //
        RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( ( @1/@0 + 1. ), @2)",
                  RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;
       //------------------------------------------------------------------------------------------



       printf("  makeLognormalConstraint : created log-normal nuisance parameter %s : val = %g\n", NP_name, np_rfv -> getVal() ) ;

       return np_rfv ;


    } // makeLognormalConstraint.


   //==============================================================================================================

    RooAbsReal* makeCorrelatedLognormalConstraint(
            const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign ) {

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* ln_mean  = new RooConstVar( vname, vname, NP_val ) ;

       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* ln_sigma = new RooConstVar( vname, vname, NP_err ) ;


       if ( NP_err <= 0. ) {
          printf("  makeCorrelatedLognormalConstraint: Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          ws_pointer -> import( *ln_mean ) ;
          ws_pointer -> import( *ln_sigma ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       char prim_name[1000] ;
       sprintf( prim_name, "prim_%s", NP_base_name ) ;
       RooRealVar* rrv_np_base_par = (RooRealVar*) allNuisances -> find( prim_name ) ;

       if ( rrv_np_base_par == 0x0 ) {

          printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter - %s\n\n", prim_name ) ;
          rrv_np_base_par = new RooRealVar( prim_name, prim_name, -6.0, 6.0 ) ;
          rrv_np_base_par -> setVal( 0. ) ;
          rrv_np_base_par -> setConstant( kFALSE ) ;
          allNuisances -> add( *rrv_np_base_par ) ;

          char vname[1000] ;
          sprintf( vname, "prim_mean_%s", NP_base_name ) ;
          RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-10.,10. ) ;
          g_mean->setConstant(kTRUE);
          sprintf( vname, "prim_sigma_%s", NP_base_name ) ;
          RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

          char pdfname[100] ;
          sprintf( pdfname, "pdf_prim_%s", NP_base_name ) ;
          printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
          RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;

          allNuisancePdfs -> add( *base_np_pdf ) ;
          globalObservables -> add( *g_mean ) ;

       }


       RooAbsReal* rar(0x0) ;


       char formula[1000] ;

       if ( !changeSign ) {
          sprintf( formula, "@0 * pow( ( @1/@0 + 1.), @2 )" ) ;
       } else {
          sprintf( formula, "@0 * pow( ( @1/@0 + 1.), -1.0 * @2 )" ) ;
       }

       rar = new RooFormulaVar( NP_name, formula, RooArgSet( *ln_mean, *ln_sigma, *rrv_np_base_par ) ) ;

       printf(" makeCorrelatedLognormalConstraint : creating correlated log-normal NP with formula : %s,  %s, val = %g, mean=%g, sigma=%g\n", formula, NP_name, rar->getVal(), NP_val, NP_err ) ;


       return rar ;

    } // makeCorrelatedLognormalConstraint.

   //==============================================================================================================




