
#include "TSystem.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>


   int fbhi[300] ;
   int fbmi[300] ;
   int fbji[300] ;
   int fbbi[300] ;
   int nfb ;

   float kqcd_ht_val[300] ;
   float kqcd_ht_err[300] ;
   float kqcd_mht_val[300] ;
   float kqcd_mht_err[300] ;
   float kqcd_njets_val[300] ;
   float kqcd_njets_err[300] ;

   int   nobs_ldp_val[300] ;
   int   nobs_zl_val[300] ;
   float nqcd_ldp_val[300] ;
   float nlostlep_ldp_val[300] ;
   float nlostlep_ldp_err[300] ;
   float nhadtau_ldp_val[300] ;
   float nhadtau_ldp_err[300] ;
   float nznunu_ldp_val[300] ;
   float nznunu_ldp_err[300] ;

   float rqcd_val[300] ;
   float frac_val[300] ;
   float nqcd_ldp_sum ;
   float nonqcd_ldp_sum_val ;
   float nonqcd_ldp_sum_err ;
   int   nobs_ldp_sum ;
   int   nobs_zl_sum ;

   float Rprime_val, Rprime_err ;
   float qcdbg_val, qcdbg_err ;


   void fill_fb_arrays( int sb_mhtht_bi, int sb_njets_bi, int sb_nbjets_bi ) ;
   void get_nqcd_ldp( ifstream& ifs, int hti, int mhti, int nji, int nbi, float& val, float& err ) ;
   void get_nbg_ldp( ifstream& ifs, int hti, int mhti, int nji, int nbi, float& val, float& err ) ;
   void get_nobs( ifstream& ifs, int hti, int mhti, int nji, int nbi, int& val_ldp, int& val_zl ) ;
   void get_line_val_err( ifstream& ifs, const char* key, float& val, float& err ) ;
   void print_kpars( bool do_table, bool include_Rerr_in_table ) ;

  //----------

   void gen_combine_table5(
                            bool sum_njets = false, bool sum_nbjets = false, bool sum_mht = false, bool sum_ht = false,
                            int sb_mhtht_bi = -1, int sb_njets_bi = -1, int sb_nbjets_bi = -1,
                            bool include_Rerr_in_table = true,
                            const char* kpars_file   = "outputfiles/qcdlhfit-results-ws-kqcd-lhfit-perfect-qcd-closure-random-nobs-with-constraints/kqcd-parameter-fit-results.txt",
                            const char* nobs_file    = "outputfiles/finebin-input-fakedata-perfect-qcd-closure-random-nobs.txt",
                            const char* lostlep_file = "outputfiles/finebin-input-lostlep.txt",
                            const char* hadtau_file  = "outputfiles/finebin-input-hadtau.txt",
                            const char* znunu_file   = "outputfiles/finebin-input-znunu.txt"
                          ) {

      char njsumstr[100] ;
      char nbsumstr[100] ;
      char mhthtsumstr[100] ;
      char configstr[1000] ;
      if ( sum_njets ) { sprintf( njsumstr, "-njsum" ) ; } else { sprintf( njsumstr, "" ) ; }
      if ( sum_nbjets ) { sprintf( nbsumstr, "-nbsum" ) ; } else { sprintf( nbsumstr, "" ) ; }
      if ( (sum_mht && sum_ht) ) { sprintf( mhthtsumstr, "-mhthtsum" ) ; } else { sprintf( mhthtsumstr, "" ) ; }
      sprintf( configstr, "qcdbg-results%s%s%s", njsumstr, nbsumstr, mhthtsumstr ) ;

      printf("\n\n %s\n\n", configstr ) ;




      ifstream ifs ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin(0.35) ;
      gStyle -> SetMarkerSize( 0.7 ) ;

      ifs.open( kpars_file ) ;
      if ( !ifs.good() ) {
         printf("\n\n *** Problem with kpars file : %s\n\n", kpars_file ) ; return ;
      }

      get_line_val_err( ifs, "Kqcd_ht1", kqcd_ht_val[1], kqcd_ht_err[1] ) ;
      get_line_val_err( ifs, "Kqcd_ht2", kqcd_ht_val[2], kqcd_ht_err[2] ) ;
      get_line_val_err( ifs, "Kqcd_ht3", kqcd_ht_val[3], kqcd_ht_err[3] ) ;

      kqcd_mht_val[1] = 1. ;  kqcd_mht_err[1] = 0. ;
      get_line_val_err( ifs, "Kqcd_mht2", kqcd_mht_val[2], kqcd_mht_err[2] ) ;
      get_line_val_err( ifs, "Kqcd_mht3", kqcd_mht_val[3], kqcd_mht_err[3] ) ;
      get_line_val_err( ifs, "Kqcd_mht4", kqcd_mht_val[4], kqcd_mht_err[4] ) ;

      kqcd_njets_val[1] = 1. ;  kqcd_njets_err[1] = 0. ;
      get_line_val_err( ifs, "Kqcd_njet2", kqcd_njets_val[2], kqcd_njets_err[2] ) ;
      get_line_val_err( ifs, "Kqcd_njet3", kqcd_njets_val[3], kqcd_njets_err[3] ) ;
      get_line_val_err( ifs, "Kqcd_njet4", kqcd_njets_val[4], kqcd_njets_err[4] ) ;
      get_line_val_err( ifs, "Kqcd_njet5", kqcd_njets_val[5], kqcd_njets_err[5] ) ;

      printf("\n\n") ;

      printf( "  Kqcd_ht1   %.3f +/- %.3f\n", kqcd_ht_val[1], kqcd_ht_err[1] ) ;
      printf( "  Kqcd_ht2   %.3f +/- %.3f\n", kqcd_ht_val[2], kqcd_ht_err[2] ) ;
      printf( "  Kqcd_ht3   %.3f +/- %.3f\n", kqcd_ht_val[3], kqcd_ht_err[3] ) ;

      printf( "  Kqcd_mht1   %.3f +/- %.3f\n", kqcd_mht_val[1], kqcd_mht_err[1] ) ;
      printf( "  Kqcd_mht2   %.3f +/- %.3f\n", kqcd_mht_val[2], kqcd_mht_err[2] ) ;
      printf( "  Kqcd_mht3   %.3f +/- %.3f\n", kqcd_mht_val[3], kqcd_mht_err[3] ) ;
      printf( "  Kqcd_mht4   %.3f +/- %.3f\n", kqcd_mht_val[4], kqcd_mht_err[4] ) ;

      printf( "  Kqcd_njets1   %.3f +/- %.3f\n", kqcd_njets_val[1], kqcd_njets_err[1] ) ;
      printf( "  Kqcd_njets2   %.3f +/- %.3f\n", kqcd_njets_val[2], kqcd_njets_err[2] ) ;
      printf( "  Kqcd_njets3   %.3f +/- %.3f\n", kqcd_njets_val[3], kqcd_njets_err[3] ) ;
      printf( "  Kqcd_njets4   %.3f +/- %.3f\n", kqcd_njets_val[4], kqcd_njets_err[4] ) ;
      printf( "  Kqcd_njets5   %.3f +/- %.3f\n", kqcd_njets_val[5], kqcd_njets_err[5] ) ;

      printf("\n\n") ;

      ifs.close() ;

      //--------


      ifstream ifs_nobs ;
      ifs_nobs.open( nobs_file ) ;
      if ( !ifs_nobs.good() ) {
         printf("\n\n *** problem opening nobs file : %s\n\n", nobs_file ) ; return ;
      }

      ifstream ifs_lostlep ;
      ifs_lostlep.open( lostlep_file ) ;
      if ( !ifs_lostlep.good() ) {
         printf("\n\n *** problem opening lostlep file : %s\n\n", lostlep_file ) ; return ;
      }

      ifstream ifs_hadtau ;
      ifs_hadtau.open( hadtau_file ) ;
      if ( !ifs_hadtau.good() ) {
         printf("\n\n *** problem opening hadtau file : %s\n\n", hadtau_file ) ; return ;
      }

      ifstream ifs_znunu ;
      ifs_znunu.open( znunu_file ) ;
      if ( !ifs_znunu.good() ) {
         printf("\n\n *** problem opening znunu file : %s\n\n", znunu_file ) ; return ;
      }


      if ( sb_mhtht_bi >=1 && sb_njets_bi >=1 && sb_nbjets_bi >=0 ) {

         //-- do it for one search bin

         nfb = 0 ;
         fill_fb_arrays( sb_mhtht_bi, sb_njets_bi, sb_nbjets_bi ) ;

         nqcd_ldp_sum = 0. ;
         nobs_ldp_sum = 0 ;
         nobs_zl_sum = 0 ;
         nonqcd_ldp_sum_val = 0. ;
         float nonqcd_ldp_sum_err2 = 0. ;
         for ( int bi=0; bi<nfb; bi++ ) {
            get_nobs(     ifs_nobs   ,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nobs_ldp_val[bi],      nobs_zl_val[bi]       ) ;
            get_nbg_ldp(  ifs_lostlep,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nlostlep_ldp_val[bi],  nlostlep_ldp_err[bi]  ) ;
            get_nbg_ldp(  ifs_hadtau ,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nhadtau_ldp_val[bi] ,  nhadtau_ldp_err[bi]   ) ;
            get_nbg_ldp(  ifs_znunu  ,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nznunu_ldp_val[bi]  ,  nznunu_ldp_err[bi]    ) ;
            float this_fb_nonqcd = nlostlep_ldp_val[bi] + nhadtau_ldp_val[bi] + nznunu_ldp_val[bi] ;
            nqcd_ldp_val[bi] = nobs_ldp_val[bi] - this_fb_nonqcd ;
            if ( nqcd_ldp_val[bi] < 0. ) {
               printf(" *** Warning: Njet%d, Nb%d, MHT%d, HT%d has negative Nobs-NnonQCD: Nobs=%d, NnonQCD=%.1f\n",
                  fbji[bi], fbbi[bi], fbmi[bi], fbhi[bi], nobs_ldp_val[bi], this_fb_nonqcd ) ;
               nqcd_ldp_val[bi] = 0. ;
            }
            nqcd_ldp_sum += nqcd_ldp_val[bi] ;
            nobs_ldp_sum += nobs_ldp_val[bi] ;
            nobs_zl_sum += nobs_zl_val[bi] ;
            nonqcd_ldp_sum_val += this_fb_nonqcd ;
            nonqcd_ldp_sum_err2 += pow( nlostlep_ldp_err[bi], 2 ) + pow( nhadtau_ldp_err[bi], 2 ) + pow( nznunu_ldp_err[bi], 2 ) ; //*** wrong (placeholder)
         } // bi
         nonqcd_ldp_sum_err = sqrt( nonqcd_ldp_sum_err2 ) ;

         printf("\n\n") ;

         printf( "                   LDP fraction     QCD high/low ratio      Nqcd for low dPhi control\n") ;
         for ( int bi=0; bi<nfb; bi++ ) {
            float kht  = kqcd_ht_val[fbhi[bi]] ;
            float kmht = kqcd_mht_val[fbmi[bi]] ;
            float knj  = kqcd_njets_val[fbji[bi]] ;
            float kht_err  = kqcd_ht_err[fbhi[bi]] ;
            float kmht_err = kqcd_mht_err[fbmi[bi]] ;
            float knj_err  = kqcd_njets_err[fbji[bi]] ;
            float rval = kht * kmht * knj ;
            float rerr = rval * sqrt( pow( kht_err / kht, 2. ) + pow( kmht_err / kmht, 2. ) + pow( knj_err / knj, 2. ) ) ;
            rqcd_val[bi] = rval ;
            float frac(0.) ;
            if ( nqcd_ldp_sum > 0 ) {
               frac = nqcd_ldp_val[bi] / nqcd_ldp_sum ;
            } else {
               frac = 1. / (1.*nfb) ;
            }
            frac_val[bi] = frac ;
            printf(" %3d : HT%d  MHT%d  NJ%d  NB%d  fldp = %.3f    R = %.3f +/- %.3f    %9.2f\n",
                 bi, fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], frac, rval, rerr, nqcd_ldp_val[bi] ) ;
         } // bi


         printf("\n\n") ;

         print_kpars( true, include_Rerr_in_table ) ;

         float nqcd_ldp_bgsub_err = sqrt( nobs_ldp_sum ) ; // *** wrong, placeholder

         qcdbg_val = nqcd_ldp_sum * Rprime_val ;
         qcdbg_err = 0. ;
         if ( Rprime_val > 0 ) {
            qcdbg_err = sqrt( pow( nqcd_ldp_sum * Rprime_err, 2 ) + pow( Rprime_val * nqcd_ldp_bgsub_err, 2 ) ) ;
         }


         int sb72ind = 24*(sb_njets_bi-1) + 6*(sb_nbjets_bi) + sb_mhtht_bi ;

         printf(" SB %2d :  QCD BG prediction =  (%.3f +/- %.3f) * (%9.2f +/- %9.2f) = %9.2f +/- %5.2f\n\n",
             sb72ind, Rprime_val, Rprime_err, nqcd_ldp_sum, nqcd_ldp_bgsub_err, qcdbg_val, qcdbg_err ) ;

      } else {

         //-- do it for several bins.

         char texfile[10000] ;

         sprintf( texfile, "outputfiles/%s-tables.tex", configstr ) ;

         printf("\n  Tables file : %s\n\n", texfile ) ;

         FILE* ofp_tex(0x0) ;
         if ( (ofp_tex=fopen( texfile, "w" ))==NULL ) {
            printf("\n\n *** Problem opening %s\n\n", texfile ) ;
            return ;
         }


         fprintf( ofp_tex, "  \\documentclass[11pt]{article}\n" ) ;
         fprintf( ofp_tex, "  \\usepackage{graphicx}\n\n" ) ;

         fprintf( ofp_tex, "   \\textwidth 7.0in\n" ) ;
         fprintf( ofp_tex, "   \\textheight 9.0in\n" ) ;
         fprintf( ofp_tex, "   \\topmargin -1.0in\n" ) ;
         fprintf( ofp_tex, "   \\linewidth 7.0in\n" ) ;
         fprintf( ofp_tex, "   \\oddsidemargin -0.3in\n" ) ;
         fprintf( ofp_tex, "   \\evensidemargin -0.3in\n\n" ) ;

         fprintf( ofp_tex, "  \\begin{document}\n\n" ) ;

         fprintf( ofp_tex, "  %%%%==========================================\n\n" ) ;

         if ( !sum_njets && !sum_nbjets && !sum_mht && !sum_ht ) {
            fprintf( ofp_tex, "\n \\scriptsize \n\n") ;
         }

         fprintf( ofp_tex, "\\begin{tabular}{|l||r|r|r|c||r|}\n" ) ;
         fprintf( ofp_tex, "\\hline\n" ) ;
         fprintf( ofp_tex, "   Bin     &    $N_{obs}^{LDP}$    &  Non-QCD   &   $N_{obs}$ - Non-QCD  &    $R_{QCD}$    &    QCD BG  \\\\\n" ) ;
         fprintf( ofp_tex, "\\hline\n" ) ;
         fprintf( ofp_tex, "\\hline\n" ) ;




         vector<TH1F*> hist_list ;

         TH1F* h_rqcd = new TH1F( "h_rqcd", "Rqcd", 100, 0.5, 100.5 ) ; hist_list.push_back( h_rqcd ) ;
         TH1F* h_nqcd_zl = new TH1F( "h_nqcd_zl", "Nqcd, ZL", 100, 0.5, 100.5 ) ; hist_list.push_back( h_nqcd_zl ) ;
         TH1F* h_nobs_ldp = new TH1F( "h_nobs_ldp", "Nobs, LDP", 100, 0.5, 100.5 ) ; hist_list.push_back( h_nobs_ldp ) ;
         TH1F* h_nobs_zl = new TH1F( "h_nobs_zl", "Nobs, ZL", 100, 0.5, 100.5 ) ; hist_list.push_back( h_nobs_zl ) ;
         TH1F* h_nnonqcd_ldp = new TH1F( "h_nnonqcd_ldp", "NnonQCD, LDP", 100, 0.5, 100.5 ) ; hist_list.push_back( h_nnonqcd_ldp ) ;
         TH1F* h_nobs_minus_nnonqcd_ldp = new TH1F( "h_nobs_minus_nnonqcd_ldp", "Nobs - NnonQCD, LDP", 100, 0.5, 100.5 ) ; hist_list.push_back( h_nobs_minus_nnonqcd_ldp ) ;

         printf("\n\n") ;
         printf("       Search bin          N_ldp      Nnonqcd +/- err      Rqcd  Rqcd*Nqldp  Kht1   Kht2   Kht3   Kmht2  Kmht3  Kmht4   Knj2   Knj3   Knj4  Knj5          Rqcd      Model Nqcd_zl\n") ;

         int sb_njmin = 1 ; int sb_njmax = 3 ;
         int sb_nbmin = 0 ; int sb_nbmax = 3 ;
         int sb_mhthtmin = 1 ; int sb_mhthtmax = 6 ;

         if ( sum_njets ) { sb_njmin = -1 ; sb_njmax = -1 ; }
         if ( sum_nbjets ) { sb_nbmin = -1 ; sb_nbmax = -1 ; }
         if ( sum_mht && sum_ht ) { sb_mhthtmin = -1 ; sb_mhthtmax = -1 ; }

         int hist_bin(1) ;

         for ( sb_njets_bi = sb_njmin; sb_njets_bi<=sb_njmax; sb_njets_bi++ ) {
            for ( sb_nbjets_bi = sb_nbmin; sb_nbjets_bi<=sb_nbmax; sb_nbjets_bi++ ) {
               for ( sb_mhtht_bi = sb_mhthtmin; sb_mhtht_bi<=sb_mhthtmax; sb_mhtht_bi++ ) {

                  int il_sb_njmin = 1 ; int il_sb_njmax = 3 ;
                  int il_sb_nbmin = 0 ; int il_sb_nbmax = 3 ;
                  int il_sb_mhthtmin = 1 ; int il_sb_mhthtmax = 6 ;

                  if ( !sum_njets ) { il_sb_njmin = sb_njets_bi ; il_sb_njmax = sb_njets_bi ; }
                  if ( !sum_nbjets ) { il_sb_nbmin = sb_nbjets_bi ; il_sb_nbmax = sb_nbjets_bi ; }
                  if ( !(sum_mht && sum_ht) ) { il_sb_mhthtmin = sb_mhtht_bi ; il_sb_mhthtmax = sb_mhtht_bi ; }


                  nfb = 0 ;
                  for ( int il_sb_nj = il_sb_njmin; il_sb_nj<=il_sb_njmax ; il_sb_nj++ ) {
                     for ( int il_sb_nb = il_sb_nbmin; il_sb_nb<=il_sb_nbmax; il_sb_nb++ ) {
                        for ( int il_sb_mhtht = il_sb_mhthtmin; il_sb_mhtht<=il_sb_mhthtmax; il_sb_mhtht++ ) {
                           fill_fb_arrays( il_sb_mhtht, il_sb_nj, il_sb_nb ) ;
                        }
                     }
                  }

                  nqcd_ldp_sum = 0. ;
                  nobs_ldp_sum = 0 ;
                  nobs_zl_sum = 0 ;
                  nonqcd_ldp_sum_val = 0. ;
                  float nonqcd_ldp_sum_err2 = 0. ;
                  for ( int bi=0; bi<nfb; bi++ ) {
                     get_nobs(     ifs_nobs   ,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nobs_ldp_val[bi]    ,  nobs_zl_val[bi]       ) ;
                     get_nbg_ldp(  ifs_lostlep,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nlostlep_ldp_val[bi],  nlostlep_ldp_err[bi]  ) ;
                     get_nbg_ldp(  ifs_hadtau ,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nhadtau_ldp_val[bi] ,  nhadtau_ldp_err[bi]   ) ;
                     get_nbg_ldp(  ifs_znunu  ,  fbhi[bi], fbmi[bi], fbji[bi], fbbi[bi], nznunu_ldp_val[bi]  ,  nznunu_ldp_err[bi]    ) ;
                     float this_fb_nonqcd = nlostlep_ldp_val[bi] + nhadtau_ldp_val[bi] + nznunu_ldp_val[bi] ;
                     nqcd_ldp_val[bi] = nobs_ldp_val[bi] - this_fb_nonqcd ;
                     if ( nqcd_ldp_val[bi] < 0. ) {
                        //////printf(" *** Warning: Njet%d, Nb%d, MHT%d, HT%d has negative Nobs-NnonQCD: Nobs=%d, NnonQCD=%.1f\n",
                        //////   fbji[bi], fbbi[bi], fbmi[bi], fbhi[bi], nobs_ldp_val[bi], this_fb_nonqcd ) ;
                        nqcd_ldp_val[bi] = 0. ;
                     }
                     nqcd_ldp_sum += nqcd_ldp_val[bi] ;
                     nobs_ldp_sum += nobs_ldp_val[bi] ;
                     nobs_zl_sum  += nobs_zl_val[bi] ;
                     nonqcd_ldp_sum_val += this_fb_nonqcd ;
                     nonqcd_ldp_sum_err2 += pow( nlostlep_ldp_err[bi], 2 ) + pow( nhadtau_ldp_err[bi], 2 ) + pow( nznunu_ldp_err[bi], 2 ) ; //*** wrong (placeholder)
                     rqcd_val[bi] = kqcd_ht_val[fbhi[bi]] * kqcd_mht_val[fbmi[bi]] * kqcd_njets_val[fbji[bi]] ;
                  } // bi
                  nonqcd_ldp_sum_err = sqrt( nonqcd_ldp_sum_err2 ) ;

                  for ( int bi=0; bi<nfb; bi++ ) {
                     float frac(0.) ;
                     if ( nqcd_ldp_sum > 0 ) {
                        frac = nqcd_ldp_val[bi] / nqcd_ldp_sum ;
                     } else {
                        frac = 1. / (1.*nfb) ;
                     }
                     frac_val[bi] = frac ;
                  } // bi

                  int sb72ind = 24*(sb_njets_bi-1) + 6*(sb_nbjets_bi) + sb_mhtht_bi ;

                  char njbinlabel[100] ;
                  char nbbinlabel[100] ;
                  char mhthtbinlabel[100] ;
                  char binnumberlabel[100] ;
                  char binlabel[100] ;
                  if ( sum_njets ) { sprintf( njbinlabel,"" ) ; } else { sprintf( njbinlabel, "-Njet%d", sb_njets_bi ) ; }
                  if ( sum_nbjets ) { sprintf( nbbinlabel,"" ) ; } else { sprintf( nbbinlabel, "-Nb%d", sb_nbjets_bi ) ; }
                  if ( sum_mht && sum_ht ) { sprintf( mhthtbinlabel,"" ) ; } else { sprintf( mhthtbinlabel, "-MHTHT%d", sb_mhtht_bi ) ; }
                  if ( sum_njets || sum_nbjets || (sum_mht && sum_ht) ) { sprintf( binnumberlabel,"" ) ; } else { sprintf( binnumberlabel,"%2d", sb72ind ) ; }
                  sprintf( binlabel, "%s SB%s%s%s", binnumberlabel, njbinlabel, nbbinlabel, mhthtbinlabel ) ;

                  printf( " %22s ", binlabel ) ;

                  print_kpars( false, include_Rerr_in_table ) ;

                  h_rqcd -> SetBinContent( hist_bin, Rprime_val ) ;
                  h_rqcd -> SetBinError( hist_bin, Rprime_err ) ;

                  h_nqcd_zl -> SetBinContent( hist_bin, qcdbg_val ) ;
                  h_nqcd_zl -> SetBinError( hist_bin, qcdbg_err ) ;

                  h_nobs_ldp -> SetBinContent( hist_bin, nobs_ldp_sum ) ;

                  h_nobs_zl  -> SetBinContent( hist_bin, nobs_zl_sum ) ;

                  h_nnonqcd_ldp -> SetBinContent( hist_bin, nonqcd_ldp_sum_val ) ;
                  h_nnonqcd_ldp -> SetBinError( hist_bin, nonqcd_ldp_sum_err ) ;

                  h_nobs_minus_nnonqcd_ldp -> SetBinContent( hist_bin, nobs_ldp_sum - nonqcd_ldp_sum_val ) ;
                  h_nobs_minus_nnonqcd_ldp -> SetBinError( hist_bin, sqrt( nobs_ldp_sum + pow( nonqcd_ldp_sum_err, 2 ) ) ) ;

                  for ( unsigned long hi=0; hi<hist_list.size(); hi++ ) { hist_list.at(hi) -> GetXaxis() -> SetBinLabel( hist_bin, binlabel ) ; }

                  fprintf( ofp_tex, " %s  &  %8d  &  $%9.1f \\pm %7.1f$   &   $%9.1f \\pm %7.1f$  &   $%5.3f \\pm %5.3f$  &   $%9.1f \\pm %7.1f$   \\\\ \n",
                    binlabel, nobs_ldp_sum, nonqcd_ldp_sum_val, nonqcd_ldp_sum_err,
                    (nobs_ldp_sum - nonqcd_ldp_sum_val), sqrt( nobs_ldp_sum + pow( nonqcd_ldp_sum_err, 2 ) ),
                    Rprime_val, Rprime_err,
                    qcdbg_val, qcdbg_err
                    ) ;

                  hist_bin++ ;

               } // sb_mhtht_bi
               if ( !(sum_mht && sum_ht) ) fprintf( ofp_tex, "\\hline\n" ) ;
            } // sb_nbjets_bi
            fprintf( ofp_tex, "\\hline\n" ) ;
         } // sb_njets_bi

         fprintf( ofp_tex, "\\hline\n" ) ;
         fprintf( ofp_tex, "\\end{tabular}\n\n" ) ;

         fprintf( ofp_tex, "  %%%%==========================================\n\n" ) ;
         fprintf( ofp_tex, "\n\n \\end{document}\n" ) ;
         fclose( ofp_tex ) ;




         int n_hist_bins = hist_bin - 1 ;

         for ( unsigned long hi=0; hi<hist_list.size(); hi++ ) {
            hist_list.at(hi) -> GetXaxis() -> LabelsOption("v") ;
            hist_list.at(hi) -> GetXaxis() -> SetRange( 1, n_hist_bins ) ;
         } // hi


         h_nobs_minus_nnonqcd_ldp -> SetMarkerStyle(20) ;

         char pdffile[10000] ;

        //------

         TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
         if ( can1 == 0x0 ) {
            can1 = new TCanvas( "can1", "", 1000, 1100 ) ;
         }
         can1 -> Clear() ;
         can1 -> Divide(1,2) ;

         h_nnonqcd_ldp -> SetMarkerStyle(20) ;
         h_nobs_ldp -> SetMinimum(0.1) ;
         h_nnonqcd_ldp -> SetFillColor( kBlue-10 ) ;
         h_nqcd_zl -> SetMarkerStyle(20) ;
         h_nqcd_zl -> SetFillColor( kRed-10 ) ;
         h_nobs_zl -> SetMinimum(0.1) ;
         h_nobs_ldp -> SetFillColor( 18 ) ;
         h_nobs_zl -> SetFillColor( 18 ) ;

         TLegend* legend_ldp = new TLegend( 0.75, 0.80, 0.98, 0.98 ) ;
         legend_ldp -> AddEntry( h_nobs_ldp, "Nobs, LDP" ) ;
         legend_ldp -> AddEntry( h_nnonqcd_ldp, "non-QCD, LDP" ) ;

         TLegend* legend_zl = new TLegend( 0.80, 0.80, 0.98, 0.98 ) ;
         legend_zl -> AddEntry( h_nobs_zl, "Nobs, ZL" ) ;
         legend_zl -> AddEntry( h_nqcd_zl, "QCD, ZL" ) ;

         h_nobs_ldp -> SetTitle( "LDP event counts" ) ;
         h_nobs_zl  -> SetTitle( "ZL event counts" ) ;

         float scale_ldp = h_nobs_ldp -> GetMaximum() ;
         float scale_zl = h_nobs_zl -> GetMaximum() ;

         can1 -> cd(1) ;
         h_nobs_ldp -> Draw() ;
         h_nnonqcd_ldp -> Draw("hist same") ;
         h_nnonqcd_ldp -> Draw("same") ;
         h_nnonqcd_ldp -> Draw("axis same") ;
         h_nnonqcd_ldp -> Draw("axig same") ;
         legend_ldp -> Draw() ;
         gPad -> SetLogy(1) ;
         gPad -> SetGridy(1) ;

         can1 -> cd(2) ;
         h_nobs_zl -> Draw("same") ;
         h_nqcd_zl -> Draw("hist same") ;
         h_nqcd_zl -> Draw("same") ;
         h_nqcd_zl -> Draw("axis same") ;
         h_nqcd_zl -> Draw("axig same") ;
         legend_zl -> Draw() ;
         gPad -> SetGridy(1) ;
         gPad -> SetLogy(1) ;

         sprintf( pdffile, "outputfiles/%s-events-logy.pdf", configstr ) ;
         can1 -> SaveAs( pdffile ) ;


        //------

         can1 -> cd(1) ;
         gPad -> SetLogy(0) ;

         can1 -> cd(2) ;
         gPad -> SetLogy(0) ;

         sprintf( pdffile, "outputfiles/%s-events-liny.pdf", configstr ) ;
         can1 -> SaveAs( pdffile ) ;


        //------

         h_nobs_ldp -> SetMaximum( 0.1 * scale_ldp ) ;
         h_nobs_zl -> SetMaximum( 0.1 * scale_zl ) ;
         can1 -> Update() ; can1 -> Draw() ;

         sprintf( pdffile, "outputfiles/%s-events-zoom1.pdf", configstr ) ;
         can1 -> SaveAs( pdffile ) ;


        //------

         h_nobs_ldp -> SetMaximum( 0.03 * scale_ldp ) ;
         h_nobs_zl -> SetMaximum( 0.03 * scale_zl ) ;
         can1 -> Update() ; can1 -> Draw() ;

         sprintf( pdffile, "outputfiles/%s-events-zoom2.pdf", configstr ) ;
         can1 -> SaveAs( pdffile ) ;

        //------

         h_nobs_ldp -> SetMaximum( 0.006 * scale_ldp ) ;
         h_nobs_zl -> SetMaximum( 0.006 * scale_zl ) ;
         can1 -> Update() ; can1 -> Draw() ;

         sprintf( pdffile, "outputfiles/%s-events-zoom3.pdf", configstr ) ;
         can1 -> SaveAs( pdffile ) ;



        //--------------------------

         h_rqcd -> SetMarkerStyle(20) ;

         TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
         if ( can3 == 0x0 ) {
            can3 = new TCanvas( "can3", "", 800, 550 ) ;
         }
         can3 -> Clear() ;

         can3 -> cd() ;
         h_rqcd -> Draw() ;
         gPad -> SetGridy(1) ;

         sprintf( pdffile, "outputfiles/%s-rqcd.pdf", configstr ) ;
         can3 -> SaveAs( pdffile ) ;

      }



   } // gen_combine_table5

  //========================================================================================================

   void fill_fb_arrays( int sb_mhtht_bi, int sb_njets_bi, int sb_nbjets_bi ) {

      if ( sb_mhtht_bi < 1 || sb_mhtht_bi > 6 ) gSystem -> Exit(0) ;
      if ( sb_njets_bi < 1 || sb_njets_bi > 3 ) gSystem -> Exit(0) ;

      /////////nfb = 0 ;
      if ( sb_njets_bi == 1 ) {
         if ( sb_mhtht_bi == 1 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
         } else if ( sb_mhtht_bi == 2 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         } else if ( sb_mhtht_bi == 3 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         } else if ( sb_mhtht_bi == 4 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         } else if ( sb_mhtht_bi == 5 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 6 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         }
      } else if ( sb_njets_bi == 2 ) {
         if ( sb_mhtht_bi == 1 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 2 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 3 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 4 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 5 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 6 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         }
      } else {
         if ( sb_mhtht_bi == 1 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 2 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 3 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 4 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 5 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 6 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         }
      }

   } // fill_fb_arrays

  //========================================================================================================

   void print_kpars( bool do_table, bool include_Rerr_in_table ) {

      float prefactor ;
      float dk_over_k ;
      float datacard_err ;
      float kval, kerr ;
      char kname[100] ;
      int hi, mi, ji ;

      Rprime_val = 0. ;
      for ( int bi=0; bi<nfb; bi++ ) { Rprime_val += frac_val[bi] * rqcd_val[bi] ; }
      if ( Rprime_val <= 0 ) { printf("\n\n *** Zero R' val???\n\n") ; return ; }

      float sum_err2(0.) ;

      qcdbg_val = nqcd_ldp_sum * Rprime_val ;

      if ( !do_table ) {

         //printf( " %9.2f   %5.3f   %9.2f ", nqcd_ldp_sum, Rprime_val, nqcd_ldp_sum * Rprime_val ) ;

         printf( " %9d   %7.2f +/- %6.2f  %5.3f   %9.2f ",
            nobs_ldp_sum, nonqcd_ldp_sum_val, nonqcd_ldp_sum_err, Rprime_val, qcdbg_val ) ;
      }


      if ( do_table ) printf(" parameter       Value        pre factor   rel error on K    What goes in card line\n") ;

      for ( int hi = 1; hi<=3; hi++ ) {
         bool is_used = false ;
         for ( int bi=0; bi<nfb; bi++ ) { if ( fbhi[bi] == hi ) is_used = true ; }
         if ( is_used ) {

            float pfsum(0.) ;
            for ( int bi=0; bi<nfb; bi++ ) {
               if ( fbhi[bi] == hi ) {
                  pfsum += frac_val[bi] * rqcd_val[bi] ;
                  /////////printf("    debug:  bi=%d, hi=%d : fldp=%.3f             Rqcd=%.3f\n",
                  /////////   bi, hi, frac_val[bi], rqcd_val[bi] ) ;
               }
            }
            sprintf( kname, "Kht%d", hi ) ;
            kerr = kqcd_ht_err[hi] ;
            kval = kqcd_ht_val[hi] ;

            prefactor = pfsum / Rprime_val ;
               /////////printf("      debug:  pfsum=%.5f, R'=%.4f, prefactor=%.4f\n",
               /////////   pfsum, Rprime_val, prefactor ) ;
            dk_over_k = kerr / kval  ;
            datacard_err = 1. + prefactor * dk_over_k ;

            sum_err2 += pow( prefactor * dk_over_k, 2 ) ;

            if ( do_table ) printf( "%7s     %.3f +/- %.3f   pre = %.3f   dk/k = %.3f     datacard err = %.3f\n",
               kname, kval, kerr, prefactor, dk_over_k, datacard_err ) ;

            if ( !do_table ) {
               if ( kerr > 0 ) {
                  printf( " %5.3f ", datacard_err ) ;
               } else {
                  printf( "   -   " ) ;
               }
            }

         } else {
            if ( !do_table) printf( "   -   " ) ;
         }
      }


      for ( int mi = 1; mi<=4; mi++ ) {
         bool is_used = false ;
         for ( int bi=0; bi<nfb; bi++ ) { if ( fbmi[bi] == mi ) is_used = true ; }
         if ( is_used ) {

            float pfsum(0.) ;
            for ( int bi=0; bi<nfb; bi++ ) { if ( fbmi[bi] == mi ) pfsum += frac_val[bi] * rqcd_val[bi] ; }
            sprintf( kname, "Kmht%d", mi ) ;
            kerr = kqcd_mht_err[mi] ;
            kval = kqcd_mht_val[mi] ;

            prefactor = pfsum / Rprime_val ;
            dk_over_k = kerr / kval  ;
            datacard_err = 1. + prefactor * dk_over_k ;

            sum_err2 += pow( prefactor * dk_over_k, 2 ) ;

            if ( do_table ) printf( "%7s     %.3f +/- %.3f   pre = %.3f   dk/k = %.3f     datacard err = %.3f\n",
               kname, kval, kerr, prefactor, dk_over_k, datacard_err ) ;

            if ( !do_table && mi>1 ) {
               if ( kerr > 0 ) {
                  printf( " %5.3f ", datacard_err ) ;
               } else {
                  printf( "   -   " ) ;
               }
            }

         } else {
            if ( mi>1 && !do_table ) printf( "   -   " ) ;
         }
      }


      for ( int ji = 1; ji<=5; ji++ ) {
         bool is_used = false ;
         for ( int bi=0; bi<nfb; bi++ ) { if ( fbji[bi] == ji ) is_used = true ; }
         if ( is_used ) {

            float pfsum(0.) ;
            for ( int bi=0; bi<nfb; bi++ ) {
               if ( fbji[bi] == ji ) {
                  pfsum += frac_val[bi] * rqcd_val[bi] ;
                  /////////printf("    debug:  bi=%d, ji=%d : fldp=%.3f             Rqcd=%.3f\n",
                  /////////   bi, ji, frac_val[bi], rqcd_val[bi] ) ;
               }
            }
            sprintf( kname, "Knjet%d", ji ) ;
            kerr = kqcd_njets_err[ji] ;
            kval = kqcd_njets_val[ji] ;

            prefactor = pfsum / Rprime_val ;
               /////////printf("      debug:  pfsum=%.5f, R'=%.4f, prefactor=%.4f\n",
               /////////   pfsum, Rprime_val, prefactor ) ;
            dk_over_k = kerr / kval  ;
            datacard_err = 1. + prefactor * dk_over_k ;

            sum_err2 += pow( prefactor * dk_over_k, 2 ) ;

            if ( do_table ) printf( "%7s     %.3f +/- %.3f   pre = %.3f   dk/k = %.3f     datacard err = %.3f\n",
               kname, kval, kerr, prefactor, dk_over_k, datacard_err ) ;

            if ( !do_table && ji>1 ) {
               if ( kerr > 0 ) {
                  printf( " %5.3f ", datacard_err ) ;
               } else {
                  printf( "   -   " ) ;
               }
            }
         } else {
            if ( ji> 1  && !do_table) printf( "   -   " ) ;
         }
      }

      if ( do_table ) printf("\n\n") ;

      Rprime_err = 0. ;
      if ( sum_err2 > 0 ) Rprime_err = Rprime_val * sqrt( sum_err2 ) ;

      float nqcd_ldp_bgsub_err = sqrt( nqcd_ldp_sum + nonqcd_ldp_sum_val ) ; // simple sqrt(N) estimate.

      qcdbg_err = 0. ;
      qcdbg_err = sqrt( pow( nqcd_ldp_sum * Rprime_err, 2 ) + pow( Rprime_val * nqcd_ldp_bgsub_err, 2 ) ) ;

      if ( do_table ) printf("     Rprime = %.3f +/- %.3f\n", Rprime_val, Rprime_err ) ;

      if ( do_table ) printf("\n\n") ;

      if ( !do_table ) {
         if ( include_Rerr_in_table) {
            printf("   %.3f +/- %.3f , %8.2f +/- %5.2f\n", Rprime_val, Rprime_err, qcdbg_val, qcdbg_err ) ;
         } else {
            printf("\n") ;
         }
      }

   } // print_kpars

  //========================================================================================================

   void get_nqcd_ldp( ifstream& ifs, int hti, int mhti, int nji, int nbi, float& val, float& err ) {

      val = 0. ;
      err = 0. ;

      ifs.seekg(0) ;
      TString ts ;

      char fbname[1000] ;
      sprintf( fbname, "FB-Njet%d-Nb%d-MHT%d-HT%d", nji, nbi, mhti, hti ) ;

      while ( ifs.good() ) {
         TString ts ;
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         if ( strcmp( first->GetString().Data(), fbname ) == 0 ) {
            TObjString* val_str = (TObjString*) tokens->At( 3 ) ;
            TObjString* err_str = (TObjString*) tokens->At( 5 ) ;
            sscanf( val_str->GetString().Data(), "%f", &val ) ;
            sscanf( err_str->GetString().Data(), "%f", &err ) ;
            return ;
         }
      }

      printf("\n\n *** can't find bin %s in nqcd file.\n\n", fbname ) ;
      gSystem -> Exit(0) ;


   } // get_nqcd_ldp

  //========================================================================================================

   void get_nll_ldp( ifstream& ifs, int hti, int mhti, int nji, int nbi, float& val, float& err ) {

      val = 0. ;
      err = 0. ;

      ifs.seekg(0) ;
      TString ts ;

      char fbname[1000] ;
      sprintf( fbname, "FB-Njet%d-Nb%d-MHT%d-HT%d", nji, nbi, mhti, hti ) ;

      while ( ifs.good() ) {
         TString ts ;
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         if ( strcmp( first->GetString().Data(), fbname ) == 0 ) {
            TObjString* val_str = (TObjString*) tokens->At( 3 ) ;
            TObjString* err_str = (TObjString*) tokens->At( 5 ) ;
            sscanf( val_str->GetString().Data(), "%f", &val ) ;
            sscanf( err_str->GetString().Data(), "%f", &err ) ;
            return ;
         }
      }

      printf("\n\n *** can't find bin %s in nqcd file.\n\n", fbname ) ;
      gSystem -> Exit(0) ;


   } // get_nll_ldp

  //========================================================================================================

   void get_nbg_ldp( ifstream& ifs, int hti, int mhti, int nji, int nbi, float& val, float& err ) {

      val = 0. ;
      err = 0. ;

      ifs.seekg(0) ;
      TString ts ;

      char fbname[1000] ;
      sprintf( fbname, "FB-Njet%d-Nb%d-MHT%d-HT%d", nji, nbi, mhti, hti ) ;

      while ( ifs.good() ) {
         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( ts.Index( fbname ) >= 0 ) {
            int bi ;
            char bname[100] ;
            sscanf( ts.Data(), "%d %s %f +/- %f", &bi, bname, &val, &err ) ;
            ///////printf(" get_nbg_ldp : found %s.  val=%f  err=%f\n", fbname, val, err ) ;
            return ;
         }
      }

      printf("\n\n *** can't find bin %s in file.\n\n", fbname ) ;
      gSystem -> Exit(0) ;


   } // get_nbg_ldp

  //========================================================================================================

   void get_nobs( ifstream& ifs, int hti, int mhti, int nji, int nbi, int& val_ldp, int& val_zl ) {

      val_ldp = 0 ;
      val_zl = 0 ;

      ifs.seekg(0) ;
      TString ts ;

      char fbname[1000] ;
      sprintf( fbname, "FB-Njet%d-Nb%d-MHT%d-HT%d", nji, nbi, mhti, hti ) ;

      while ( ifs.good() ) {
         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( ts.Index( fbname ) >= 0 ) {
            int bi ;
            char bname[100] ;
            sscanf( ts.Data(), "%d %s %d %d", &bi, bname, &val_ldp, &val_zl ) ;
            ///////printf(" get_nobs : found %s.  val=%d\n", fbname, val ) ;
            return ;
         }
      }

      printf("\n\n *** can't find bin %s in file.\n\n", fbname ) ;
      gSystem -> Exit(0) ;


   } // get_nobs

  //========================================================================================================

   void get_line_val_err( ifstream& ifs, const char* key, float& val, float& err ) {

      val = 0. ;
      err = 0. ;

      ifs.seekg(0) ;
      TString ts ;

      int lc(0) ;
      while ( ifs.good() ) {
         ts.ReadLine( ifs ) ;
         if ( ts.Index( key ) >= 0 ) {
            //printf(" found key %s : line = %s\n", key, ts.Data() ) ;
            if ( ts.Index( "constrained" ) >= 0 ) {
               char pname[100] ;
               float mean ;
               sscanf( ts.Data(), "%s %f constrained by %f +/- %f", pname, &val, &mean, &err ) ;
               //printf("    constrained mean=%f, err=%f\n", mean, err ) ;
            } else {
               char pname[100] ;
               sscanf( ts.Data(), "%s %f +/- %f", pname, &val, &err ) ;
            }
            ///////printf("  get_line_val_err: found %s, val=%f, err = %f\n", key, val, err ) ;
            return ;
         }
        //-------
      }

      printf("\n\n *** can't fine line with key %s\n\n", key ) ;
      gSystem -> Exit(0) ;


   } // get_line_val_err

  //========================================================================================================




