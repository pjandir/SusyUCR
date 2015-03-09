
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

  //=================================================================================

   void build_ra2b_ws2( const char* infile = "outputfiles/lhbuilder-input-t1bbbbH.txt",
                        const char* outfile = "outputfiles/ws-t1bbbbH.root",
                        float min_signal_frac = 0.15,
                        bool skip_testfit = true,
                        bool skip_modelconfig = true,
                        float saveall_below_N = 0.,
                        bool no_rounding = true ) {

      char pname[100] ;
      char pname2[100] ;

      ifstream ifs ;
      ifs.open( infile ) ;
      if ( !ifs.good() ) {
         printf("\n\n *** Problem with input file : %s\n\n", infile ) ; return ;
      }

      int n_fine_njet_bins = find_line_val( ifs, "FB-Njet-Binning" ) ;
      int n_fine_nb_bins   = find_line_val( ifs, "FB-Nb-Binning" ) ;
      int n_fine_mht_bins  = find_line_val( ifs, "FB-MHT-Binning" ) ;
      int n_fine_ht_bins   = find_line_val( ifs, "FB-HT-Binning" ) ;

      const int MAX_FB(1000) ;
      char fb_name[MAX_FB][100] ;
      float fb_rsl_zl_val[MAX_FB] ;
      float fb_rsl_zl_err[MAX_FB] ;
      float fb_rsl_ldp_val[MAX_FB] ;
      float fb_rsl_ldp_err[MAX_FB] ;
      float fb_nzl[MAX_FB] ;
      float fb_nsl[MAX_FB] ;
      float fb_nldp[MAX_FB] ;
      float fb_sig_nzl[MAX_FB] ;
      float fb_sig_nsl[MAX_FB] ;
      float fb_sig_nldp[MAX_FB] ;

      int n_fine_bins = find_line_val( ifs, "N-fine-bins" ) ;
      printf("\n Number of fine bins is %d\n", n_fine_bins ) ;
      for ( int fbi=0; fbi<n_fine_bins; fbi++ ) {
         TString line ;
         line.ReadLine( ifs ) ;
         printf("  Fine bin %3d : %s\n", fbi, line.Data() ) ;
         sscanf( line.Data(), "%s %f %f  %f %f  %f %f %f  %f %f %f",
           fb_name[fbi],
           &(fb_rsl_zl_val[fbi]),
           &(fb_rsl_zl_err[fbi]),
           &(fb_rsl_ldp_val[fbi]),
           &(fb_rsl_ldp_err[fbi]),
           &(fb_nzl[fbi]),
           &(fb_nsl[fbi]),
           &(fb_nldp[fbi]),
           &(fb_sig_nzl[fbi]),
           &(fb_sig_nsl[fbi]),
           &(fb_sig_nldp[fbi]) ) ;
      } // fbi.

      for ( int fbi=0; fbi<n_fine_bins; fbi++ ) {
         printf( " %3d : %20s  %7.1f %7.1f %7.1f\n", fbi, fb_name[fbi], fb_nzl[fbi], fb_nsl[fbi], fb_nldp[fbi] ) ;
      } // fbi.


      const int MAX_SB(1000) ;
      char sb_name[MAX_SB][100] ;

      int n_search_bins = find_line_val( ifs, "N-search-bins" ) ;
      printf("\n Number of search bins is %d\n", n_search_bins ) ;

      for ( int sbi=0; sbi<n_search_bins; sbi++ ) {
         TString line ;
         line.ReadLine( ifs ) ;
         sprintf( sb_name[sbi], "%s", line.Data() ) ;
      } // sbi

      for ( int sbi=0; sbi<n_search_bins; sbi++ ) {
         printf( " %3d : %30s\n", sbi, sb_name[sbi] ) ;
      } // sbi



      int sb_nfb[MAX_SB] ;
      int sb_fbi[MAX_SB][10] ;
      for ( int sbi=0; sbi<n_search_bins; sbi++ ) {
         sb_nfb[sbi] = 0 ;
         for ( int i=0; i<10; i++ ) sb_fbi[sbi][i] = -1 ;
      } // sbi


      if ( !find_line( ifs, "Fine-bin-search-bin-map" ) ) {
         printf("\n\n *** can't find map.\n\n") ; return ;
      }
      int fb_sbi[MAX_SB] ;
      for ( int fbi=0; fbi<n_fine_bins; fbi++ ) {
         TString line ;
         line.ReadLine( ifs ) ;
         char fbname[100] ;
         char sbname[100] ;
         sscanf( line, "%s %s", fbname, sbname ) ;
         if ( strcmp( fbname, fb_name[fbi] ) != 0 ) { printf("\n\n *** Inconsistency %s %s\n", fbname, fb_name[fbi] ) ; return ; }
         if ( strcmp( sbname, "X" ) == 0 ) {
            fb_sbi[fbi] = -1 ;
         } else {
            for ( int sbi=0; sbi<n_search_bins; sbi++ ) {
               bool found(false) ;
               if ( strcmp( sbname, sb_name[sbi] ) == 0 ) {
                  fb_sbi[fbi] = sbi ;
                  sb_fbi[sbi][sb_nfb[sbi]] = fbi ;
                  sb_nfb[sbi] ++ ;
                  found = true ;
                  break ;
               }
               if ( !found ) {
                  fb_sbi[fbi] = -1 ;
               }
            } // sbi
         }
         printf( " %3d %30s - %3d ", fbi, fb_name[fbi], fb_sbi[fbi] ) ;
         if ( fb_sbi[fbi] >= 0 ) {
            printf( "%30s\n", sb_name[fb_sbi[fbi]] ) ;
         } else {
            printf( "not used\n") ;
         }
      } // fbi.


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




      int fb_qcd_ht_par_ind[10] ;
      int fb_qcd_mht_par_ind[10] ;
      int fb_qcd_njet_par_ind[10] ;
      int fb_qcd_nb_par_ind[10] ;

      for ( int i=1; i<=n_fine_ht_bins; i++ ) {
         sprintf( pname, "FB-HT%d-QCD-HT-par-ind", i ) ;
         fb_qcd_ht_par_ind[i] = find_line_val( ifs, pname ) ;
      }

      for ( int i=1; i<=n_fine_mht_bins; i++ ) {
         sprintf( pname, "FB-MHT%d-QCD-MHT-par-ind", i ) ;
         fb_qcd_mht_par_ind[i] = find_line_val( ifs, pname ) ;
      }

      for ( int i=1; i<=n_fine_njet_bins; i++ ) {
         sprintf( pname, "FB-Njet%d-QCD-Njet-par-ind", i ) ;
         fb_qcd_njet_par_ind[i] = find_line_val( ifs, pname ) ;
      }

      for ( int i=0; i<n_fine_nb_bins; i++ ) {
         sprintf( pname, "FB-Nb%d-QCD-Nb-par-ind", i ) ;
         fb_qcd_nb_par_ind[i] = find_line_val( ifs, pname ) ;
      }



      printf( "\n\n ======== Done reading input from %s\n\n\n", infile ) ;









      printf("\n\n Creating workspace.\n\n") ;

      RooWorkspace workspace("ws") ;
      workspace.autoImportClassCode(true) ;

      globalObservables      = new RooArgSet("globalObservables");
      allNuisances           = new RooArgSet("allNuisances");
      allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      pdf_sbIndexList        = new RooArgSet("pdf_sbIndexList") ;
      allBGmuPars            = new RooArgSet("allBGmuPars") ;
      RooArgSet* observedParametersList = new RooArgSet("observables") ;

      RooArgSet pdflist ;

      sprintf( pname, "sig_strength" ) ;
      RooRealVar* rv_sig_strength = new RooRealVar( pname, pname, 1.0, 0., 10. ) ;
      rv_sig_strength -> setConstant(kFALSE) ;
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
            rv_qcd_kht[i] = new RooRealVar( pname, pname, qcd_kht_val[i], 0., 10. ) ;
         } else {
            rv_qcd_kht[i] = makeLognormalConstraint( pname, qcd_kht_val[i], qcd_kht_err[i] ) ;
         }
      } // i

      for ( int i=1; i<=n_qcd_kmht_pars; i++ ) {
         sprintf( pname, "Kqcd_mht%d", i ) ;
         if ( qcd_kmht_err[i] < 0 ) {
            rv_qcd_kmht[i] = new RooRealVar( pname, pname, qcd_kmht_val[i], 0., 10. ) ;
         } else {
            rv_qcd_kmht[i] = makeLognormalConstraint( pname, qcd_kmht_val[i], qcd_kmht_err[i] ) ;
         }
      } // i

      for ( int i=1; i<=n_qcd_knjet_pars; i++ ) {
         sprintf( pname, "Kqcd_njet%d", i ) ;
         if ( qcd_knjet_err[i] < 0 ) {
            rv_qcd_knjet[i] = new RooRealVar( pname, pname, qcd_knjet_val[i], 0., 10. ) ;
         } else {
            rv_qcd_knjet[i] = makeLognormalConstraint( pname, qcd_knjet_val[i], qcd_knjet_err[i] ) ;
         }
      } // i

      for ( int i=1; i<=n_qcd_knb_pars; i++ ) {
         sprintf( pname, "Kqcd_nb%d", i ) ;
         if ( qcd_knb_err[i] < 0 ) {
            rv_qcd_knb[i] = new RooRealVar( pname, pname, qcd_knb_val[i], 0., 10. ) ;
         } else {
            rv_qcd_knb[i] = makeLognormalConstraint( pname, qcd_knb_val[i], qcd_knb_err[i] ) ;
         }
      } // i



      printf("\n\n ========== Begin loop over search bins. ==============\n\n") ;

      for ( int sbi=0; sbi<n_search_bins; sbi++ ) {


         printf(" %3d : %30s :\n", sbi, sb_name[sbi] ) ;

         sprintf( pname, "sb_index_%s", sb_name[sbi] ) ;
         RooConstVar* rv_sb_index = new RooConstVar( pname, pname, sbi ) ;
         sbIndexList.add( *rv_sb_index ) ;

         float Nzl_sum(0.) ;

         RooArgList ral_mu_qcd_zl ;
         RooArgList ral_mu_ll_zl ;
         RooArgList ral_mu_sig0_zl ;

         {
            float nobs(0.) ;
            float nsig(0.) ;
            float nbg(0.) ;
            for ( int i=0; i<sb_nfb[sbi]; i++ ) {
               int fbi = sb_fbi[sbi][i] ;
               nobs += fb_nzl[fbi] ;
               nsig += fb_sig_nzl[fbi] ;
               nbg  += nobs - nsig ; // assuming input has signal strength set to 1.
            }
            float sig_frac(0.) ;
            if ( nobs > 0 ) sig_frac = nsig / nobs ;
            float q = 2.*(sqrt(nbg+nsig)-sqrt(nbg)) ;
            if ( sig_frac < min_signal_frac && nobs > saveall_below_N ) {
               printf("\n %35s : Signal fraction low:  Q = %5.3f,   Nsig = %8.2f,  Nobs = %8.2f\n\n", sb_name[sbi], q, nsig, nobs ) ;
               continue ;
            }
            printf("\n %35s : Keeping this one:  Q = %5.3f,   Nsig = %8.2f,  Nobs = %8.2f\n\n", sb_name[sbi], q, nsig, nobs ) ;

         }

         for ( int i=0; i<sb_nfb[sbi]; i++ ) {

            int fbi = sb_fbi[sbi][i] ;
            int fb_nji(-1), fb_nbi(-1), fb_mbi(-1), fb_hbi(-1) ;
            sscanf( fb_name[fbi], "FB-Njet%d-Nb%d-MHT%d-HT%d", &fb_nji, &fb_nbi, &fb_mbi, &fb_hbi ) ;

            printf( "\n ------------------------------------------------------------------------\n") ;
            printf( "      %s : fb_nji=%d, fb_nbi=%d, fb_mbi=%d, fb_hbi=%d\n", fb_name[fbi], fb_nji, fb_nbi, fb_mbi, fb_hbi ) ;
            printf( "      %s : observables : Nldp = %8.1f, Nsl = %8.1f, Nzl = %8.1f\n", 
                          fb_name[fbi], fb_nldp[fbi], fb_nsl[fbi], fb_nzl[fbi] ) ;
            printf( "      %s : QCD model : Knjet%d * Knb%d * Kmht%d * Kht%d = %4.2f * %4.2f * %4.2f * %4.2f = %5.3f\n",
                fb_name[fbi],
                fb_qcd_njet_par_ind[fb_nji], fb_qcd_nb_par_ind[fb_nbi], fb_qcd_mht_par_ind[fb_mbi], fb_qcd_ht_par_ind[fb_hbi],
                qcd_knjet_val[ fb_qcd_njet_par_ind[fb_nji] ],
                qcd_knb_val[   fb_qcd_nb_par_ind[fb_nbi] ],
                qcd_kmht_val[  fb_qcd_mht_par_ind[fb_mbi] ],
                qcd_kht_val[   fb_qcd_ht_par_ind[fb_hbi] ],
                ( qcd_knjet_val[ fb_qcd_njet_par_ind[fb_nji] ] * qcd_knb_val[   fb_qcd_nb_par_ind[fb_nbi] ] * qcd_kmht_val[  fb_qcd_mht_par_ind[fb_mbi] ] * qcd_kht_val[   fb_qcd_ht_par_ind[fb_hbi] ] )
                ) ;
            printf("\n") ;

            Nzl_sum += fb_nzl[fbi] ;





            printf("\n --- Single Lepton:\n") ;

            sprintf( pname, "mu_ll_sl_%s", fb_name[fbi] ) ;
            RooRealVar* rv_mu_ll_sl = new RooRealVar( pname, pname, fb_nsl[fbi], 0., get_par_max( fb_nsl[fbi] ) ) ;
            rv_mu_ll_sl -> setConstant( kFALSE ) ;
            rv_mu_ll_sl -> Print() ;
            allBGmuPars -> add( *rv_mu_ll_sl ) ;

            sprintf( pname, "mu_sig0_sl_%s", fb_name[fbi] ) ;
            RooRealVar* rv_mu_sig0_sl = new RooRealVar( pname, pname, fb_sig_nsl[fbi], 0., 1.e6 ) ;
            rv_mu_sig0_sl -> setConstant( kTRUE ) ;
            rv_mu_sig0_sl -> Print() ;

            sprintf( pname, "mu_sig_sl_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_mu_sig_sl = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_sl ) ) ;
            rv_mu_sig_sl -> Print() ;

            sprintf( pname, "n_sl_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_n_sl = new RooFormulaVar( pname, "@0 + @1", RooArgSet( *rv_mu_ll_sl, *rv_mu_sig_sl ) ) ;
            rv_n_sl -> Print() ;

            sprintf( pname, "Nsl_%s", fb_name[fbi] ) ;
            RooRealVar* rv_Nsl = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_Nsl -> setVal( fb_nsl[fbi] ) ;
            rv_Nsl -> setConstant( kTRUE ) ;
            rv_Nsl -> Print() ;
            observedParametersList -> add( *rv_Nsl ) ;

            sprintf( pname, "pdf_sl_%s", fb_name[fbi] ) ;
            //////RooPoisson* rv_pdf_sl = new RooPoisson( pname, pname, *rv_Nsl, *rv_n_sl, no_rounding ) ;
            RooPoissonLogEval* rv_pdf_sl = new RooPoissonLogEval( pname, pname, *rv_Nsl, *rv_n_sl, no_rounding ) ;
            rv_pdf_sl -> Print() ;

            pdflist.add( *rv_pdf_sl ) ;

            sprintf( pname2, "%s_sb_index", pname ) ;
            RooConstVar* rv_pdf_sb_index_sl = new RooConstVar( pname2, pname2, sbi ) ;
            pdf_sbIndexList -> add( *rv_pdf_sb_index_sl ) ;




            printf("\n --- Low Delta Phi:\n") ;

            float rsl_ldp_val = fb_rsl_ldp_val[fbi] ;
            float rsl_ldp_err = fb_rsl_ldp_err[fbi] ;
            if ( rsl_ldp_val <= 0. ) {
               rsl_ldp_val = 1.0 ;
               rsl_ldp_err = 1.0 ;
            } else if ( rsl_ldp_err > 3 ) {
               rsl_ldp_val = 2.0 ;
               rsl_ldp_err = 2.0 ;
            }

            sprintf( pname, "R_sl_ldp_%s", fb_name[fbi] ) ;
            RooAbsReal* rv_R_sl_ldp = makeLognormalConstraint( pname, rsl_ldp_val, rsl_ldp_err, sbi ) ;
            rv_R_sl_ldp -> Print() ;

            sprintf( pname, "mu_ll_ldp_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_mu_ll_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_R_sl_ldp, *rv_mu_ll_sl ) ) ;
            rv_mu_ll_ldp -> Print() ;

            float initial_qcd_ldp_val = 0. ;
            if ( fb_nldp[fbi] > rv_mu_ll_ldp -> getVal() ) {
               initial_qcd_ldp_val = fb_nldp[fbi] - rv_mu_ll_ldp -> getVal() ;
            }
            sprintf( pname, "mu_qcd_ldp_%s", fb_name[fbi] ) ;
            RooRealVar* rv_mu_qcd_ldp = new RooRealVar( pname, pname, initial_qcd_ldp_val, 0., get_par_max( fb_nldp[fbi] ) ) ;
            rv_mu_qcd_ldp -> setConstant( kFALSE ) ;
            rv_mu_qcd_ldp -> Print() ;
            allBGmuPars -> add( *rv_mu_qcd_ldp ) ;

            sprintf( pname, "mu_sig0_ldp_%s", fb_name[fbi] ) ;
            RooRealVar* rv_mu_sig0_ldp = new RooRealVar( pname, pname, fb_sig_nldp[fbi], 0., 1.e6 ) ;
            rv_mu_sig0_ldp -> setConstant( kTRUE ) ;
            rv_mu_sig0_ldp -> Print() ;

            sprintf( pname, "mu_sig_ldp_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_mu_sig_ldp = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_ldp ) ) ;
            rv_mu_sig_ldp -> Print() ;

            sprintf( pname, "n_ldp_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_n_ldp = new RooFormulaVar( pname, "@0 + @1 + @2", RooArgSet( *rv_mu_qcd_ldp, *rv_mu_ll_ldp, *rv_mu_sig_ldp ) ) ;
            rv_n_ldp -> Print() ;

            sprintf( pname, "Nldp_%s", fb_name[fbi] ) ;
            RooRealVar* rv_Nldp = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_Nldp -> setVal( fb_nldp[fbi] ) ;
            rv_Nldp -> setConstant( kTRUE ) ;
            rv_Nldp -> Print() ;
            observedParametersList -> add( *rv_Nldp ) ;

            sprintf( pname, "pdf_ldp_%s", fb_name[fbi] ) ;
            /////////RooPoisson* rv_pdf_ldp = new RooPoisson( pname, pname, *rv_Nldp, *rv_n_ldp, no_rounding ) ;
            RooPoissonLogEval* rv_pdf_ldp = new RooPoissonLogEval( pname, pname, *rv_Nldp, *rv_n_ldp, no_rounding ) ;
            rv_pdf_ldp -> Print() ;

            pdflist.add( *rv_pdf_ldp ) ;

            sprintf( pname2, "%s_sb_index", pname ) ;
            RooConstVar* rv_pdf_sb_index_ldp = new RooConstVar( pname2, pname2, sbi ) ;
            pdf_sbIndexList -> add( *rv_pdf_sb_index_ldp ) ;




            printf("\n --- Zero Lepton:\n") ;

            float rsl_zl_val = fb_rsl_zl_val[fbi] ;
            float rsl_zl_err = fb_rsl_zl_err[fbi] ;
            if ( rsl_zl_val <= 0. ) {
               rsl_zl_val = 1.0 ;
               rsl_zl_err = 1.0 ;
            } else if ( rsl_zl_err > 3 ) {
               rsl_zl_val = 1.0 ;
               rsl_zl_err = 1.0 ;
            }

            sprintf( pname, "R_sl_zl_%s", fb_name[fbi] ) ;
            RooAbsReal* rv_R_sl_zl = makeLognormalConstraint( pname, rsl_zl_val, rsl_zl_err, sbi ) ;
            rv_R_sl_zl -> Print() ;

            sprintf( pname, "mu_ll_zl_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_mu_ll_zl = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_R_sl_zl, *rv_mu_ll_sl ) ) ;
            rv_mu_ll_zl -> Print() ;

            ral_mu_ll_zl.add( *rv_mu_ll_zl ) ;

            sprintf( pname, "R_qcd_ldp_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_R_qcd_ldp = new RooFormulaVar( pname, "@0 * @1 * @2 * @3", 
               RooArgSet(
                  *(rv_qcd_kht[ fb_qcd_ht_par_ind[fb_hbi] ]),
                  *(rv_qcd_kmht[ fb_qcd_mht_par_ind[fb_mbi] ]),
                  *(rv_qcd_knjet[ fb_qcd_njet_par_ind[fb_nji] ]),
                  *(rv_qcd_knb[ fb_qcd_nb_par_ind[fb_nbi] ])
                        ) ) ;
            rv_R_qcd_ldp -> Print() ;

            sprintf( pname, "mu_qcd_zl_%s", fb_name[fbi] ) ;
            RooFormulaVar* rv_mu_qcd_zl = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_R_qcd_ldp, *rv_mu_qcd_ldp ) ) ;
            rv_mu_qcd_zl -> Print() ;

            ral_mu_qcd_zl.add( *rv_mu_qcd_zl ) ;


            sprintf( pname, "mu_sig0_zl_%s", fb_name[fbi] ) ;
            RooRealVar* rv_mu_sig0_zl = new RooRealVar( pname, pname, fb_sig_nzl[fbi], 0., 1.e6 ) ;
            rv_mu_sig0_zl -> setConstant( kTRUE ) ;
            rv_mu_sig0_zl -> Print() ;

            ral_mu_sig0_zl.add( *rv_mu_sig0_zl ) ;


         } // i

         printf("\n Search bin BG sums and Zero Lepton PDF.\n") ;

        //--- Zero Lepton section

         char sumstring[1000] ;
         sprintf( sumstring, "@0" ) ;
         for ( int i=1; i<sb_nfb[sbi]; i++ ) {
            char tmpstr[1000] ;
            sprintf( tmpstr, "%s", sumstring ) ;
            sprintf( sumstring, "%s + @%d", tmpstr, i ) ;
         } // i

         sprintf( pname, "mu_qcd_zl_%s", sb_name[sbi] ) ;
         RooFormulaVar* rv_mu_qcd_zl = new RooFormulaVar( pname, sumstring, ral_mu_qcd_zl ) ;
         rv_mu_qcd_zl -> Print() ;

         sprintf( pname, "mu_ll_zl_%s", sb_name[sbi] ) ;
         RooFormulaVar* rv_mu_ll_zl = new RooFormulaVar( pname, sumstring, ral_mu_ll_zl ) ;
         rv_mu_ll_zl -> Print() ;

         sprintf( pname, "mu_allbg_zl_%s", sb_name[sbi] ) ;
         RooFormulaVar* rv_mu_allbg_zl = new RooFormulaVar( pname, "@0 + @1", RooArgSet( *rv_mu_qcd_zl, *rv_mu_ll_zl ) ) ;
         rv_mu_allbg_zl -> Print() ;

         sprintf( pname, "mu_sig0_zl_%s", sb_name[sbi] ) ;
         RooFormulaVar* rv_mu_sig0_zl = new RooFormulaVar( pname, sumstring, ral_mu_sig0_zl ) ;
         rv_mu_sig0_zl -> Print() ;

         sprintf( pname, "mu_sig_zl_%s", sb_name[sbi] ) ;
         RooFormulaVar* rv_mu_sig_zl = new RooFormulaVar( pname, "@0 * @1", RooArgSet( *rv_sig_strength, *rv_mu_sig0_zl ) ) ;
         rv_mu_sig_zl -> Print() ;

         sprintf( pname, "n_zl_%s", sb_name[sbi] ) ;
         RooFormulaVar* rv_n_zl = new RooFormulaVar( pname, "@0 + @1", RooArgSet( *rv_mu_allbg_zl, *rv_mu_sig_zl ) ) ;
         rv_n_zl -> Print() ;

         sprintf( pname, "Nzl_%s", sb_name[sbi] ) ;
         RooRealVar* rv_Nzl = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rv_Nzl -> setVal( Nzl_sum ) ;
         rv_Nzl -> setConstant( kTRUE ) ;
         rv_Nzl -> Print() ;
         observedParametersList -> add( *rv_Nzl ) ;

         sprintf( pname, "pdf_zl_%s", sb_name[sbi] ) ;
         /////////////RooPoisson* rv_pdf_zl = new RooPoisson( pname, pname, *rv_Nzl, *rv_n_zl, no_rounding ) ;
         RooPoissonLogEval* rv_pdf_zl = new RooPoissonLogEval( pname, pname, *rv_Nzl, *rv_n_zl, no_rounding ) ;
         rv_pdf_zl -> Print() ;

         pdflist.add( *rv_pdf_zl ) ;

         sprintf( pname2, "%s_sb_index", pname ) ;
         RooConstVar* rv_pdf_sb_index_zl = new RooConstVar( pname2, pname2, sbi ) ;
         pdf_sbIndexList -> add( *rv_pdf_sb_index_zl ) ;

         printf("\n") ;

      } // sbi


      printf("\n\n ========== End loop over search bins. ==============\n\n") ;


      printf("\n List of search bins:\n") ;
      //workspace.import( sbIndexList ) ;
      workspace.defineSet( "sbIndexList", sbIndexList, kTRUE ) ;

      printf("\n List of pdf-to-SB index const pars:\n") ;
      //workspace.import( *pdf_sbIndexList ) ;
      workspace.defineSet( "pdf_sbIndexList", *pdf_sbIndexList, kTRUE ) ;
      printf("\n") ;





      printf("\n\n Creating and importing dataset into workspace.\n\n") ;

      RooDataSet* dsObserved = new RooDataSet("observed_rds", "observed data values", *observedParametersList ) ;
      dsObserved -> add( *observedParametersList ) ;
      workspace.import( *dsObserved ) ;



      pdflist.add( *allNuisancePdfs ) ;

      printf("\n List of all PDFs\n") ;
      pdflist.Print() ;
      printf("\n") ;

      ///////RooProdPdf* likelihood = new RooProdPdf( "likelihood", "likelihood", pdflist ) ;
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








   } // build_ra2b_ws2

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


    RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err, int sbi ) {

       if ( NP_err <= 0. ) {
          printf(" makeLognormalConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
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



       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

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


  //=================================================================================

    float get_par_max( float mu ) {

       if ( mu <= 1. ) return 10. ;
       if ( mu <= 2. ) return 15. ;
       if ( mu <= 50. ) return mu + 4*sqrt(mu) ;
       if ( mu <= 100. ) return mu + 3*sqrt(mu) ;
       return mu + 2*sqrt(mu) ;


    } // get_par_max

  //=================================================================================





