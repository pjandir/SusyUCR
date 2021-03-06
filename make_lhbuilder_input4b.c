

#include "TString.h"
#include "TSystem.h"

#include "histio.c"
#include <stdio.h>

   TH1F* get_hist( const char* hname ) ;

  //---------------------------

   void make_lhbuilder_input4b( bool perfect_closure = true,
                               float true_sig_strength = 1.0,
                            const char* infile = "outputfiles/fill-bg-hists4-t1bbbbH-postdraw.root",
                            const char* signame = "t1bbbbH",
                            const char* outfilebase = "outputfiles/lhbuilder-input-v4b",
                            bool setup_qcdlhfit = false,
                            bool include_mht_ratios = false
                           ) {


      printf("\n\n") ;
      printf("  Input root file: %s\n", infile ) ;
      if ( setup_qcdlhfit ) { printf("   Setting up QCD LH fit.\n") ; } else { printf("  Not setting up QCD LH fit.\n") ; }

      printf("\n\n") ;

      gDirectory -> Delete( "*" ) ;


      float qcd_kht_val[10] ;
      float qcd_kht_err[10] ;
      float qcd_kmht_val[10] ;
      float qcd_kmht_err[10] ;
      float qcd_knjet_val[10] ;
      float qcd_knjet_err[10] ;
      float qcd_knb_val[10] ;
      float qcd_knb_err[10] ;

   //----------------------------
   //
   //   Pars for minDeltaPhiN > 6.0, no isotrk veto
   //
////  int n_qcd_kht_pars(3) ;
////  qcd_kht_val[1] = 0.156 ;  qcd_kht_err[1] = 0.02 ;
////  qcd_kht_val[2] = 0.087 ;  qcd_kht_err[2] = 0.01 ;
////  qcd_kht_val[3] = 0.043 ;  qcd_kht_err[3] = 0.01 ;

////  int n_qcd_kmht_pars ;
////  n_qcd_kmht_pars = 4 ;
////  qcd_kmht_val[1] = 1    ;  qcd_kmht_err[1] = 0    ;
////  qcd_kmht_val[2] = 1.34 ;  qcd_kmht_err[2] = 0.2  ;
////  qcd_kmht_val[3] = 1.93 ;  qcd_kmht_err[3] = 1.0  ;
////  qcd_kmht_val[4] = 2.2  ;  qcd_kmht_err[4] = 2.0  ;

////  int n_qcd_knjet_pars(5) ;
////  qcd_knjet_val[1] = 1    ;  qcd_knjet_err[1] = 0    ;
////  qcd_knjet_val[2] = 0.58 ;  qcd_knjet_err[2] = 0.1  ;
////  qcd_knjet_val[3] = 0.38 ;  qcd_knjet_err[3] = 0.1  ;
////  qcd_knjet_val[4] = 0.44 ;  qcd_knjet_err[4] = 0.2  ;
////  qcd_knjet_val[5] = 0.53 ;  qcd_knjet_err[5] = 0.5  ;

////  int n_qcd_knb_pars(1) ;
////  qcd_knb_val[1] = 1    ;  qcd_knb_err[1] = 0    ;

   //----------------------------
   //
   //  These are for minDeltaPhi_MHT > 0.5
   //

      int n_qcd_kht_pars(3) ;
      qcd_kht_val[1] = 0.062 ;  qcd_kht_err[1] = 0.010 ;
      qcd_kht_val[2] = 0.048 ;  qcd_kht_err[2] = 0.008 ;
      qcd_kht_val[3] = 0.036 ;  qcd_kht_err[3] = 0.004 ;

      int n_qcd_kmht_pars ;
      n_qcd_kmht_pars = 4 ;
      qcd_kmht_val[1] = 1    ;  qcd_kmht_err[1] = 0    ;
      qcd_kmht_val[2] = 0.472 ;  qcd_kmht_err[2] = 0.08 ;
      qcd_kmht_val[3] = 0.328 ;  qcd_kmht_err[3] = 0.16 ;
      qcd_kmht_val[4] = 0.308 ;  qcd_kmht_err[4] = 0.30 ;

      int n_qcd_knjet_pars(5) ;
      qcd_knjet_val[1] = 1    ;  qcd_knjet_err[1] = 0    ;
      qcd_knjet_val[2] = 1.45   ;  qcd_knjet_err[2] = 0.15 ;
      qcd_knjet_val[3] = 1.45   ;  qcd_knjet_err[3] = 0.15 ;
      qcd_knjet_val[4] = 2.11   ;  qcd_knjet_err[4] = 0.3  ;
      qcd_knjet_val[5] = 4.0    ;  qcd_knjet_err[5] = 2.0  ;

      int n_qcd_knb_pars(1) ;
      qcd_knb_val[1] = 1    ;  qcd_knb_err[1] = 0    ;

   //----------------------------
   //
   //  These are for minDeltaPhi30 > 0.4
   //

///   int n_qcd_kht_pars(3) ;
///   qcd_kht_val[1] = 0.080 ;  qcd_kht_err[1] = 0.024 ;
///   qcd_kht_val[2] = 0.064 ;  qcd_kht_err[2] = 0.009 ;
///   qcd_kht_val[3] = 0.053 ;  qcd_kht_err[3] = 0.004 ;

///   int n_qcd_kmht_pars ;
///   n_qcd_kmht_pars = 4 ;
///   qcd_kmht_val[1] = 1    ;  qcd_kmht_err[1] = 0    ;
///   qcd_kmht_val[2] = 0.443 ;  qcd_kmht_err[2] = 0.05 ;
///   qcd_kmht_val[3] = 0.305 ;  qcd_kmht_err[3] = 0.10 ;
///   qcd_kmht_val[4] = 0.426 ;  qcd_kmht_err[4] = 0.30 ;

///   int n_qcd_knjet_pars(5) ;
///   qcd_knjet_val[1] = 1    ;  qcd_knjet_err[1] = 0    ;
///   qcd_knjet_val[2] = 1.31   ;  qcd_knjet_err[2] = 0.1  ;
///   qcd_knjet_val[3] = 1.49   ;  qcd_knjet_err[3] = 0.2  ;
///   qcd_knjet_val[4] = 1.66   ;  qcd_knjet_err[4] = 0.3  ;
///   qcd_knjet_val[5] = 3.5    ;  qcd_knjet_err[5] = 1.0  ;

///   int n_qcd_knb_pars(1) ;
///   qcd_knb_val[1] = 1    ;  qcd_knb_err[1] = 0    ;

   //----------------------------



      int n_rmht(0) ;
      char rmht_name[500][100] ;

      if ( include_mht_ratios ) {

         sprintf( rmht_name[n_rmht], "Rmht_sbnj1_ht1_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj1_ht2_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj1_ht2_fbmht4over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj1_ht3_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj1_ht3_fbmht4over2" ) ; n_rmht ++ ;

         sprintf( rmht_name[n_rmht], "Rmht_sbnj2_ht1_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj2_ht2_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj2_ht2_fbmht4over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj2_ht3_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj2_ht3_fbmht4over2" ) ; n_rmht ++ ;

         sprintf( rmht_name[n_rmht], "Rmht_sbnj3_ht1_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj3_ht2_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj3_ht2_fbmht4over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj3_ht3_fbmht3over2" ) ; n_rmht ++ ;
         sprintf( rmht_name[n_rmht], "Rmht_sbnj3_ht3_fbmht4over2" ) ; n_rmht ++ ;

      }



      int fb_qcd_ht_par_ind[10] ;
      int fb_qcd_mht_par_ind[10] ;
      int fb_qcd_njet_par_ind[10] ;
      int fb_qcd_nb_par_ind[10] ;

      fb_qcd_ht_par_ind[1] = 1 ;
      fb_qcd_ht_par_ind[2] = 2 ;
      fb_qcd_ht_par_ind[3] = 3 ;

      fb_qcd_mht_par_ind[1] = 1 ;
      fb_qcd_mht_par_ind[2] = 2 ;
      fb_qcd_mht_par_ind[3] = 3 ;
      fb_qcd_mht_par_ind[4] = 4 ;

      fb_qcd_njet_par_ind[1] = 1 ;
      fb_qcd_njet_par_ind[2] = 2 ;
      fb_qcd_njet_par_ind[3] = 3 ;
      fb_qcd_njet_par_ind[4] = 4 ;
      fb_qcd_njet_par_ind[5] = 5 ;

      fb_qcd_nb_par_ind[0] = 1 ;
      fb_qcd_nb_par_ind[1] = 1 ;
      fb_qcd_nb_par_ind[2] = 1 ;
      fb_qcd_nb_par_ind[3] = 1 ;




      char hname[1000] ;
      TH1F* hp ;

      loadHist( infile ) ;

      char outfile[10000] ;
      if ( setup_qcdlhfit ) {
         sprintf( outfile, "%s-%s-qcdlhfit.txt", outfilebase, signame ) ;
      } else {
         sprintf( outfile, "%s-%s.txt", outfilebase, signame ) ;
      }
      FILE* outfp ;
      printf("\n\n Creating output file: %s\n\n", outfile ) ;
      if ( (outfp=fopen( outfile, "w" ))==NULL ) {
         printf("\n\n *** can't open output file.\n\n") ; return ;
      }

      FILE* outfp_qcdcounts ;
      sprintf( outfile, "%s-qcdcounts.txt", outfilebase ) ;
      if ( (outfp_qcdcounts=fopen( outfile, "w" ))==NULL ) { printf("\n\n *** can't open %s file.\n\n", outfile) ; return ; }

      FILE* outfp_llcounts ;
      sprintf( outfile, "%s-llcounts.txt", outfilebase ) ;
      if ( (outfp_llcounts=fopen( outfile, "w" ))==NULL ) { printf("\n\n *** can't open %s file.\n\n", outfile) ; return ; }


      hp = get_hist( "h_binning_njets_bg" ) ;
      int nnjetbins = hp->GetNbinsX() ;
      double njetbins[10] ;
      printf("FB-Njet-Binning  %d ", nnjetbins-1 ) ;
      fprintf( outfp, "FB-Njet-Binning  %d ", nnjetbins-1 ) ;
      for ( int bi=1; bi<=hp->GetNbinsX(); bi++ ) {
         printf(" %.1f", hp->GetXaxis()->GetBinUpEdge(bi) ) ;
         fprintf( outfp, " %.1f", hp->GetXaxis()->GetBinUpEdge(bi) ) ;
         njetbins[bi] = hp->GetXaxis()->GetBinUpEdge(bi) ;
      }
      printf( "\n" ) ;
      fprintf( outfp, "\n" ) ;

      hp = get_hist( "h_binning_nb_bg" ) ;
      int nnbbins = hp->GetNbinsX() ;
      double nbbins[10] ;
      printf("FB-Nb-Binning    %d ", nnbbins ) ;
      fprintf( outfp, "FB-Nb-Binning    %d ", nnbbins ) ;
      for ( int bi=1; bi<=hp->GetNbinsX(); bi++ ) {
         printf("  %.1f", hp->GetXaxis()->GetBinLowEdge(bi) ) ;
         fprintf( outfp, "  %.1f", hp->GetXaxis()->GetBinLowEdge(bi) ) ;
         nbbins[bi-1] = hp->GetXaxis()->GetBinLowEdge(bi) ;
      }
      printf( "  %.1f\n", hp->GetXaxis()->GetBinUpEdge(hp->GetNbinsX()) ) ;
      fprintf( outfp, "  %.1f\n", hp->GetXaxis()->GetBinUpEdge(hp->GetNbinsX()) ) ;
      nbbins[nnbbins] = hp->GetXaxis()->GetBinUpEdge(nnbbins) ;

      hp = get_hist( "h_binning_mht_bg" ) ;
      int nmhtbins = hp->GetNbinsX() ;
      double mhtbins[10] ;
      printf("FB-MHT-Binning   %d ", nmhtbins-1 ) ;
      fprintf( outfp, "FB-MHT-Binning   %d ", nmhtbins-1 ) ;
      for ( int bi=1; bi<=hp->GetNbinsX(); bi++ ) {
         printf("  %.0f", hp->GetXaxis()->GetBinUpEdge(bi) ) ;
         fprintf( outfp, "  %.0f", hp->GetXaxis()->GetBinUpEdge(bi) ) ;
         mhtbins[bi] = hp->GetXaxis()->GetBinUpEdge(bi) ;
      }
      printf( "\n" ) ;
      fprintf( outfp, "\n" ) ;

      hp = get_hist( "h_binning_ht_bg" ) ;
      int nhtbins = hp->GetNbinsX() ;
      double htbins[10] ;
      printf("FB-HT-Binning    %d ", nhtbins-1 ) ;
      fprintf( outfp, "FB-HT-Binning    %d ", nhtbins-1 ) ;
      for ( int bi=1; bi<=hp->GetNbinsX(); bi++ ) {
         printf("  %.0f", hp->GetXaxis()->GetBinUpEdge(bi) ) ;
         fprintf( outfp, "  %.0f", hp->GetXaxis()->GetBinUpEdge(bi) ) ;
         htbins[bi] = hp->GetXaxis()->GetBinUpEdge(bi) ;
      }
      printf( "\n" ) ;
      fprintf( outfp, "\n" ) ;



      printf( "N-fine-bins %d\n", (nnjetbins-1)*(nnbbins)*(nmhtbins-1)*(nhtbins-1) ) ;
      fprintf( outfp, "N-fine-bins %d\n", (nnjetbins-1)*(nnbbins)*(nmhtbins-1)*(nhtbins-1) ) ;


      hp = get_hist( "h_njvsmhtvsht_nb0_zl_ttbar_nj1_1d" ) ;
      int nhbins = hp -> GetNbinsX() ;



      float zl_ll_count[10][10][10] ; // htbi, mbi, nji
      float zl_ll_err2[10][10][10] ;

      for ( int nji=1; nji<nnjetbins; nji++ ) {
         for ( int mbi=1; mbi<nmhtbins; mbi++ ) {
            for ( int htbi=1; htbi<nhtbins; htbi++ ) {
               zl_ll_count[htbi][mbi][nji] = 0. ;
               zl_ll_err2[htbi][mbi][nji] = 0. ;
            } // htbi
         } // mbi
      } // nji


      float rslldp_ldp_nj1_nb0_val[10][10] ;
      /////////float rsl_ldp_nj1_nb0_val[10][10] ;

      const int N_FB_MAX(1000) ;
      int n_fine_bins(0) ;
      char fb_name[N_FB_MAX][100] ;


      for ( int nji=1; nji<nnjetbins; nji++ ) {

         for ( int nbi=0; nbi<nnbbins; nbi++ ) {

            for ( int hbi=1; hbi<nhbins; hbi++ ) {

               bool emptyBin(false) ;


               sprintf( hname, "h_njvsmhtvsht_nb%d_ldp_qcd_nj%d_1d", nbi, nji ) ;
               TH1F* hp_qcd_ldp = get_hist( hname ) ;
               sprintf( hname, "h_njvsmhtvsht_nb%d_zl_qcd_nj%d_1d", nbi, nji ) ;
               TH1F* hp_qcd_zl = get_hist( hname ) ;

               sprintf( hname, "h_njvsmhtvsht_nb%d_ldp_%s_nj%d_1d", nbi, signame, nji ) ;
               TH1F* hp_sig_ldp = get_hist( hname ) ;
               sprintf( hname, "h_njvsmhtvsht_nb%d_zl_%s_nj%d_1d", nbi, signame, nji ) ;
               TH1F* hp_sig_zl = get_hist( hname ) ;
               sprintf( hname, "h_njvsmhtvsht_nb%d_sl_%s_nj%d_1d", nbi, signame, nji ) ;
               TH1F* hp_sig_sl = get_hist( hname ) ;
               sprintf( hname, "h_njvsmhtvsht_nb%d_slldp_%s_nj%d_1d", nbi, signame, nji ) ;
               TH1F* hp_sig_slldp = get_hist( hname ) ;



               sprintf( hname, "h_bgsum_nb%d_ldp_nj%d", nbi, nji ) ;
               TH1F* hp_bgsum_ldp = get_hist( hname ) ;
               sprintf( hname, "h_bgsum_nb%d_zl_nj%d", nbi, nji ) ;
               TH1F* hp_bgsum_zl = get_hist( hname ) ;
               sprintf( hname, "h_bgsum_nb%d_sl_nj%d", nbi, nji ) ;
               TH1F* hp_bgsum_sl = get_hist( hname ) ;
               sprintf( hname, "h_bgsum_nb%d_slldp_nj%d", nbi, nji ) ;
               TH1F* hp_bgsum_slldp = get_hist( hname ) ;



               sprintf( hname, "h_llbgsum_nb%d_ldp_nj%d", nbi, nji ) ;
               TH1F* hp_llbgsum_ldp = get_hist( hname ) ;
               sprintf( hname, "h_llbgsum_nb%d_zl_nj%d", nbi, nji ) ;
               TH1F* hp_llbgsum_zl = get_hist( hname ) ;
               sprintf( hname, "h_llbgsum_nb%d_sl_nj%d", nbi, nji ) ;
               TH1F* hp_llbgsum_sl = get_hist( hname ) ;
               sprintf( hname, "h_llbgsum_nb%d_slldp_nj%d", nbi, nji ) ;
               TH1F* hp_llbgsum_slldp = get_hist( hname ) ;







               float nsl_val = hp_bgsum_sl -> GetBinContent( hbi ) ;
               float nsl_err = hp_bgsum_sl -> GetBinError( hbi ) ;

               float nldp_val = hp_bgsum_ldp -> GetBinContent( hbi ) ;
               float nldp_err = hp_bgsum_ldp -> GetBinError( hbi ) ;

               float nzl_val = hp_bgsum_zl -> GetBinContent( hbi ) ;
               float nzl_err = hp_bgsum_zl -> GetBinError( hbi ) ;

               float nslldp_val = hp_bgsum_slldp -> GetBinContent( hbi ) ;
               float nslldp_err = hp_bgsum_slldp -> GetBinError( hbi ) ;



               float nsig_sl_val = hp_sig_sl -> GetBinContent( hbi ) ;
               float nsig_sl_err = hp_sig_sl -> GetBinError( hbi ) ;

               float nsig_ldp_val = hp_sig_ldp -> GetBinContent( hbi ) ;
               float nsig_ldp_err = hp_sig_ldp -> GetBinError( hbi ) ;

               float nsig_zl_val = hp_sig_zl -> GetBinContent( hbi ) ;
               float nsig_zl_err = hp_sig_zl -> GetBinError( hbi ) ;

               float nsig_slldp_val = hp_sig_slldp -> GetBinContent( hbi ) ;
               float nsig_slldp_err = hp_sig_slldp -> GetBinError( hbi ) ;




               float nll_sl_val = hp_llbgsum_sl -> GetBinContent( hbi ) ;
               float nll_sl_err = hp_llbgsum_sl -> GetBinError( hbi ) ;

               float nll_zl_val = hp_llbgsum_zl -> GetBinContent( hbi ) ;
               float nll_zl_err = hp_llbgsum_zl -> GetBinError( hbi ) ;

               float nll_ldp_val = hp_llbgsum_ldp -> GetBinContent( hbi ) ;
               float nll_ldp_err = hp_llbgsum_ldp -> GetBinError( hbi ) ;

               float nll_slldp_val = hp_llbgsum_slldp -> GetBinContent( hbi ) ;
               float nll_slldp_err = hp_llbgsum_slldp -> GetBinError( hbi ) ;




               float nqcd_zl_val = hp_qcd_zl -> GetBinContent( hbi ) ;
               float nqcd_zl_err = hp_qcd_zl -> GetBinError( hbi ) ;

               float nqcd_ldp_val = hp_qcd_ldp -> GetBinContent( hbi ) ;
               float nqcd_ldp_err = hp_qcd_ldp -> GetBinError( hbi ) ;





               char binlabel[100] ;

               sprintf( binlabel, "%s", hp_llbgsum_zl -> GetXaxis() -> GetBinLabel( hbi ) ) ;
               if ( strlen( binlabel ) == 0 ) { emptyBin = true ; continue ; }
               int htbin(-1), mhtbin(-1) ;
               sscanf( binlabel, "MHT%d_HT%d", &mhtbin, &htbin ) ;
               TString blts( binlabel ) ;
               blts.ReplaceAll("_","-") ;



               float rsl_zl_val = 1. ;
               float rsl_zl_err = 1. ;
               if ( nll_sl_val > 0 ) {
                  rsl_zl_val = nll_zl_val / nll_sl_val ;
                  if ( nll_zl_val > 0 ) {
                     rsl_zl_err = rsl_zl_val * sqrt( pow( nll_sl_err/nll_sl_val, 2 ) + pow( nll_zl_err/nll_zl_val, 2 ) ) ;
                  }
               }

               float rslldp_ldp_val = 1. ;
               float rslldp_ldp_err = 1. ;
               if ( nll_slldp_val > 0 ) {
                  rslldp_ldp_val = nll_ldp_val / nll_slldp_val ;
                  if ( nll_ldp_val > 0 ) {
                     rslldp_ldp_err = rslldp_ldp_val * sqrt( pow( nll_slldp_err/nll_slldp_val, 2 ) + pow( nll_ldp_err/nll_ldp_val, 2 ) ) ;
                  }
                  if ( nll_ldp_val <= 0 ) {
                     rslldp_ldp_val = rslldp_ldp_nj1_nb0_val[mhtbin][htbin] ;
                     rslldp_ldp_err = 2. * rslldp_ldp_val ;
                     nll_ldp_val = rslldp_ldp_val * nll_slldp_val ;
                  }
               }

      ///////  float rsl_ldp_val = 1. ;
      ///////  float rsl_ldp_err = 1. ;
      ///////  if ( nll_sl_val > 0 ) {
      ///////     rsl_ldp_val = nll_ldp_val / nll_sl_val ;
      ///////     if ( nll_ldp_val > 0 ) {
      ///////        rsl_ldp_err = rsl_ldp_val * sqrt( pow( nll_sl_err/nll_sl_val, 2 ) + pow( nll_ldp_err/nll_ldp_val, 2 ) ) ;
      ///////     }
      ///////     if ( nll_ldp_val <= 0 ) {
      ///////        rsl_ldp_val = rsl_ldp_nj1_nb0_val[mhtbin][htbin] ;
      ///////        rsl_ldp_err = 2. * rsl_ldp_val ;
      ///////        nll_ldp_val = rsl_ldp_val * nll_sl_val ;
      ///////     }
      ///////  }

               if ( nji==1 && nbi==0 ) {
                  rslldp_ldp_nj1_nb0_val[mhtbin][htbin] = rslldp_ldp_val ;
                  ////////rsl_ldp_nj1_nb0_val[mhtbin][htbin] = rsl_ldp_val ;
               }

               float calc_nll_zl = rsl_zl_val * nll_sl_val ;
               printf("     lost lepton ZL  :  MC = %10.3f ,  calc = %10.3f\n", nll_zl_val, calc_nll_zl ) ;

               float calc_nll_ldp = rslldp_ldp_val * nll_slldp_val ;
               printf("     lost lepton LDP :  MC = %10.3f ,  calc = %10.3f\n", nll_ldp_val, calc_nll_ldp ) ;

               float qcd_ratio =    qcd_kht_val[     fb_qcd_ht_par_ind[htbin]    ]
                                  * qcd_kmht_val[    fb_qcd_mht_par_ind[mhtbin]  ]
                                  * qcd_knjet_val[   fb_qcd_njet_par_ind[nji]    ]
                                  * qcd_knb_val[     fb_qcd_nb_par_ind[nbi]      ]     ;

               float calc_nqcd_zl = qcd_ratio * nqcd_ldp_val ;

               printf("   QCD Ratio   ht,mht,njet,nb=(%d,%d,%d,%d) :  %5.3f * %5.3f * %5.3f * %5.3f  = %5.3f,  Nzl calc = %7.1f,  mc = %7.1f\n",
                  fb_qcd_ht_par_ind[htbin], fb_qcd_mht_par_ind[mhtbin], fb_qcd_njet_par_ind[nji], fb_qcd_nb_par_ind[nbi],
                  qcd_kht_val[     fb_qcd_ht_par_ind[htbin]    ],
                  qcd_kmht_val[    fb_qcd_mht_par_ind[mhtbin]  ],
                  qcd_knjet_val[   fb_qcd_njet_par_ind[nji]    ],
                  qcd_knb_val[     fb_qcd_nb_par_ind[nbi]      ],
                  qcd_ratio,
                  calc_nqcd_zl, nqcd_zl_val) ;

               float calc_nzl_bg = calc_nll_zl + calc_nqcd_zl ;
               float mc_nzl_bg   = nll_zl_val + nqcd_zl_val ;

               float calc_nldp_bg = calc_nll_ldp + nqcd_ldp_val ;
               float mc_nldp_bg   = nll_ldp_val + nqcd_ldp_val ;

               printf("   Nzl bg  calc = %10.3f (ll=%10.3f, qcd=%10.3f),  mc = %10.3f\n",  calc_nzl_bg, nll_zl_val, calc_nqcd_zl, mc_nzl_bg ) ;

               float nzl_output ;
               float nldp_output ;
               if ( perfect_closure ) {
                  nzl_output = calc_nzl_bg + true_sig_strength * nsig_zl_val ;
                  nldp_output = calc_nldp_bg + true_sig_strength * nsig_ldp_val ;
                  zl_ll_count[htbin][mhtbin][nji] += calc_nll_zl ;
                  zl_ll_err2[htbin][mhtbin][nji] += nll_zl_err * nll_zl_err ;
               } else {
                  nzl_output =   mc_nzl_bg + true_sig_strength * nsig_zl_val ;
                  nldp_output =  mc_nldp_bg + true_sig_strength * nsig_ldp_val ;
                  zl_ll_count[htbin][mhtbin][nji] += nll_zl_val ;
                  zl_ll_err2[htbin][mhtbin][nji] += nll_zl_err * nll_zl_err ;
               }

               printf("   Nzl total = %10.3f   (calc bg=%10.3f, sig=%10.3f)\n", nzl_output, calc_nzl_bg, true_sig_strength * nsig_zl_val  ) ;

               float nsl_output = nll_sl_val + true_sig_strength * nsig_sl_val ;
               ////////float nldp_output = nll_ldp_val + nqcd_ldp_val + true_sig_strength * nsig_ldp_val ;

               float nslldp_output = nll_slldp_val + true_sig_strength * nsig_slldp_val ;



               sprintf( fb_name[n_fine_bins], "FB-Njet%d-Nb%d-%s", nji, nbi, blts.Data() ) ;

               printf(         "%s  ", fb_name[n_fine_bins] ) ;
               fprintf( outfp, "%s  ", fb_name[n_fine_bins] ) ;

               printf(         " %6.3f %6.3f    %6.3f %6.3f ", rsl_zl_val, rsl_zl_err,  rslldp_ldp_val, rslldp_ldp_err ) ;
               fprintf( outfp, " %6.3f %6.3f    %6.3f %6.3f ", rsl_zl_val, rsl_zl_err,  rslldp_ldp_val, rslldp_ldp_err ) ;


               printf(         "   %9.3f  %9.3f  %9.3f  %9.3f  ", nzl_output, nsl_output, nldp_output, nslldp_output ) ;
               fprintf( outfp, "   %9.3f  %9.3f  %9.3f  %9.3f  ", nzl_output, nsl_output, nldp_output, nslldp_output ) ;

               printf(         "   %9.3f  %9.3f  %9.3f  %9.3f", nsig_zl_val, nsig_sl_val, nsig_ldp_val, nsig_slldp_val ) ;
               fprintf( outfp, "   %9.3f  %9.3f  %9.3f  %9.3f", nsig_zl_val, nsig_sl_val, nsig_ldp_val, nsig_slldp_val ) ;

               if ( perfect_closure ) {
                  fprintf( outfp_qcdcounts, "FB-Njet%d-Nb%d-%s   QCD LDP   %9.2f +/- %6.2f      QCD ZL   %9.2f (pefect closure)\n",
                      nji, nbi, blts.Data(), nqcd_ldp_val, nqcd_ldp_err,   calc_nqcd_zl ) ;
               } else {
                  fprintf( outfp_qcdcounts, "FB-Njet%d-Nb%d-%s   QCD LDP   %9.2f +/- %6.2f      QCD ZL   %9.2f +/- %6.2f\n",
                      nji, nbi, blts.Data(), nqcd_ldp_val, nqcd_ldp_err,   nqcd_zl_val, nqcd_zl_err ) ;
               }

               /////////fprintf( outfp_llcounts, "FB-Njet%d-Nb%d-%s   LL ZL   %9.2f +/- %5.2f\n",
              /////////     nji, nbi, blts.Data(), nll_zl_val, nll_zl_err ) ;
               fprintf( outfp_llcounts, "FB-Njet%d-Nb%d-%s     LL LDP   %9.2f +/- %6.2f     LL SLLDP  %9.2f +/- %6.2f    LL ZL  %9.2f +/- %6.2f\n",
                   nji, nbi, blts.Data(), nll_ldp_val, nll_ldp_err,  nll_slldp_val, nll_slldp_err, nll_zl_val, nll_zl_err ) ;

               printf("\n") ;
               fprintf( outfp, "\n") ;

               n_fine_bins++ ;
               if ( n_fine_bins > N_FB_MAX ) {
                  printf("\n\n *** too many fine bins.\n\n" ) ; return ;
               }


            } // hbi

         } // nbi

      } // nji


      int n_search_njet = 3 ;
      int n_search_nb = 4 ;
      int n_search_mhtht = 6 ;

      char search_njet_name[3][10] = { "Njet456", "Njet78", "Njet9+" } ;
      char search_nb_name[4][10] = { "Nb0", "Nb1", "Nb2", "Nb3+" } ;
      char search_mhtht_name[6][10] = { "MHT1-HT1", "MHT1-HT2", "MHT1-HT3", "MHT2-HT12", "MHT2-HT3", "MHT3-HT23" } ;

      float search_njet_bins[4] = { 3.5, 6.5, 8.5, 20. } ;
      float search_nb_bins[5] = { -0.5, 0.5, 1.5, 2.5, 20. } ;
      float search_mhtht_mhtlow[6]  = { 200., 200., 200., 500., 500.,   750. } ;
      float search_mhtht_mhthigh[6] = { 500., 500., 500., 750., 750., 20000. } ;
      float search_mhtht_htlow[6]  = { 500.,  800.,  1200.,  500.,  1200.,   800. } ;
      float search_mhtht_hthigh[6] = { 800., 1200., 20000., 1200., 20000., 20000. } ;

      int n_search_bins = n_search_njet * n_search_nb * n_search_mhtht ;

      if ( setup_qcdlhfit ) {

         printf( "N-search-bins  %d\n", n_fine_bins ) ;
         fprintf( outfp, "N-search-bins  %d\n", n_fine_bins ) ;

         for ( int bi=0; bi<n_fine_bins; bi++ ) {
            TString sbname( fb_name[bi] ) ;
            sbname.ReplaceAll( "FB", "SB" ) ;
            printf( "%s\n", sbname.Data() ) ;
            fprintf( outfp, "%s\n", sbname.Data() ) ;
         } // bi

         printf( "Fine-bin-search-bin-map\n" ) ;
         fprintf( outfp, "Fine-bin-search-bin-map\n" ) ;

         for ( int bi=0; bi<n_fine_bins; bi++ ) {
            TString sbname( fb_name[bi] ) ;
            sbname.ReplaceAll( "FB", "SB" ) ;
            printf( "%s   %s\n", fb_name[bi], sbname.Data() ) ;
            fprintf( outfp, "%s   %s\n", fb_name[bi], sbname.Data() ) ;
         } // bi


      } else {

         printf( "N-search-bins  %d\n", n_search_bins ) ;
         fprintf( outfp, "N-search-bins  %d\n", n_search_bins ) ;

         for ( int nji=1; nji<=n_search_njet; nji++ ) {
            for ( int nbi=0; nbi<n_search_nb; nbi++ ) {
               for ( int i=0; i<n_search_mhtht; i++ ) {
                  printf( "SB-%s-%s-%s\n", search_njet_name[nji-1], search_nb_name[nbi], search_mhtht_name[i] ) ;
                  fprintf( outfp, "SB-%s-%s-%s\n", search_njet_name[nji-1], search_nb_name[nbi], search_mhtht_name[i] ) ;
               } // i
            } // nbi
         } // nji

         printf( "Fine-bin-search-bin-map\n" ) ;
         fprintf( outfp, "Fine-bin-search-bin-map\n" ) ;

         for ( int nji=1; nji<nnjetbins; nji++ ) {
            for ( int nbi=0; nbi<nnbbins; nbi++ ) {
               for ( int mbi=1; mbi<nmhtbins; mbi++ ) {
                  for ( int htbi=1; htbi<nhtbins; htbi++ ) {
                     printf( "FB-Njet%d-Nb%d-MHT%d-HT%d  ", nji, nbi, mbi, htbi ) ;
                     fprintf( outfp, "FB-Njet%d-Nb%d-MHT%d-HT%d  ", nji, nbi, mbi, htbi ) ;
                     int sb_nji(-1) ;
                     int sb_nbi(-1) ;
                     int sb_mhthtbi(-1) ;
                     for ( int i=0; i<n_search_njet; i++ ) { if (  njetbins[nji]>= search_njet_bins[i] && njetbins[nji+1]<= search_njet_bins[i+1] ) sb_nji = i ; }
                     for ( int i=0; i<n_search_nb; i++ ) { if (  nbbins[nbi]>= search_nb_bins[i] && nbbins[nbi+1]<= search_nb_bins[i+1] ) sb_nbi = i ; }
                     for ( int i=0; i<n_search_mhtht; i++ ) {
                        if ( mhtbins[mbi]>= search_mhtht_mhtlow[i] && mhtbins[mbi+1] <= search_mhtht_mhthigh[i]
                          && htbins[htbi]>= search_mhtht_htlow[i] && htbins[htbi+1] <= search_mhtht_hthigh[i] ) sb_mhthtbi = i ;
                     }
                     if ( sb_nji>-1 && sb_nbi>-1 && sb_mhthtbi > -1 ) {
                        printf( " SB-%s-%s-%s", search_njet_name[sb_nji], search_nb_name[sb_nbi], search_mhtht_name[sb_mhthtbi] ) ;
                        fprintf( outfp, " SB-%s-%s-%s", search_njet_name[sb_nji], search_nb_name[sb_nbi], search_mhtht_name[sb_mhthtbi] ) ;
                     } else {
                        printf( " X" ) ;
                        fprintf( outfp, " X" ) ;
                     }
                     printf("\n") ;
                     fprintf( outfp, "\n") ;

                  } // htbi
               } // mbi
            } // nbi
         } // nji

      }


      fprintf( outfp, "N-QCD-Kht-pars    %d\n", n_qcd_kht_pars ) ;
      for ( int i=1; i<=n_qcd_kht_pars; i++ ) {
         float err = qcd_kht_err[i] ;
         if ( setup_qcdlhfit ) err = -1. ;
         fprintf( outfp, "QCD-Kht%d    %6.3f  %6.3f\n", i, qcd_kht_val[i], err ) ;
      }
      fprintf( outfp, "N-QCD-Kmht-pars    %d\n", n_qcd_kmht_pars ) ;
      int startind(1) ;
      int lastind = n_qcd_kmht_pars ;
      for ( int i=startind; i<=lastind; i++ ) {
         float err = qcd_kmht_err[i] ;
         if ( setup_qcdlhfit ) err = -1. ;
         if ( i==startind ) err = 0. ;
         fprintf( outfp, "QCD-Kmht%d    %6.3f  %6.3f\n", i, qcd_kmht_val[i], err ) ;
      }
      fprintf( outfp, "N-QCD-Knjet-pars    %d\n", n_qcd_knjet_pars ) ;
      for ( int i=1; i<=n_qcd_knjet_pars; i++ ) {
         float err = qcd_knjet_err[i] ;
         if ( setup_qcdlhfit ) err = -1. ;
         if ( i==1 ) err = 0. ;
         fprintf( outfp, "QCD-Knjet%d    %6.3f  %6.3f\n", i, qcd_knjet_val[i], err ) ;
      }
      fprintf( outfp, "N-QCD-Knb-pars    %d\n", n_qcd_knb_pars ) ;
      for ( int i=1; i<=n_qcd_knb_pars; i++ ) {
         float err = qcd_knb_err[i] ;
         if ( setup_qcdlhfit ) err = -1. ;
         if ( i==1 ) err = 0. ;
         fprintf( outfp, "QCD-Knb%d    %6.3f  %6.3f\n", i, qcd_knb_val[i], err ) ;
      }


      fprintf( outfp, "FB-Njet1-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[1] ) ;
      fprintf( outfp, "FB-Njet2-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[2] ) ;
      fprintf( outfp, "FB-Njet3-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[3] ) ;
      fprintf( outfp, "FB-Njet4-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[4] ) ;
      fprintf( outfp, "FB-Njet5-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[5] ) ;

      for ( int bi=0; bi<nnbbins; bi++ ) {
         fprintf( outfp, "FB-Nb%d-QCD-Nb-par-ind    %d\n", bi, fb_qcd_nb_par_ind[bi] ) ;
      }

      fprintf( outfp, "FB-MHT1-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[1] ) ;
      fprintf( outfp, "FB-MHT2-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[2] ) ;
      fprintf( outfp, "FB-MHT3-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[3] ) ;
      fprintf( outfp, "FB-MHT4-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[4] ) ;

      fprintf( outfp, "FB-HT1-QCD-HT-par-ind    %d\n", fb_qcd_ht_par_ind[1] ) ;
      fprintf( outfp, "FB-HT2-QCD-HT-par-ind    %d\n", fb_qcd_ht_par_ind[2] ) ;
      fprintf( outfp, "FB-HT3-QCD-HT-par-ind    %d\n", fb_qcd_ht_par_ind[3] ) ;


      if ( include_mht_ratios ) {

         fprintf( outfp, "N-Rmht  %d\n", n_rmht ) ;

         for ( int ri=0; ri<n_rmht; ri++ ) {

            int rnji(-1), rhti(-1), rmbin(-1), rmbid(-1) ;
            sscanf( rmht_name[ri], "Rmht_sbnj%d_ht%d_fbmht%dover%d", &rnji, &rhti, &rmbin, &rmbid ) ;
            printf("  %3d : %40s : nji=%d, hti=%d, mbin=%d, mbid=%d\n", ri, rmht_name[ri], rnji, rhti, rmbin, rmbid ) ;

            float r_numer(0.), r_denom(0.) ;
            float r_numer_err2(0.), r_denom_err2(0.) ;
            for ( int nji=1; nji<nnjetbins; nji++ ) {
               int sb_nji(-1) ;
               for ( int i=0; i<n_search_njet; i++ ) { if (  njetbins[nji]>= search_njet_bins[i] && njetbins[nji+1]<= search_njet_bins[i+1] ) sb_nji = i ; }
               if ( sb_nji <0 ) { printf( "\n\n *** wtf? sb_nji = %d\n", sb_nji ) ; continue ; }
               sb_nji += 1 ; // counting starts at 1 for this in the name.
               ////////printf("  nji=%d, sb_nji=%d, rnji=%d\n", nji, sb_nji, rnji ) ;
               if ( sb_nji != rnji ) continue ;
               for ( int htbi=1; htbi<nhtbins; htbi++ ) {
                  ////////printf("     htbi=%d, rhti=%d\n", htbi, rhti ) ;
                  if ( htbi != rhti ) continue ;
                  printf( "      Numerator   rhti=%d, Njet%d, MHT%d = %8.2f +/- %8.2f\n", rhti, nji, rmbin, zl_ll_count[rhti][rmbin][nji], sqrt(zl_ll_err2[rhti][rmbin][nji]) ) ;
                  printf( "      Denominator rhti=%d, Njet%d, MHT%d = %8.2f +/- %8.2f\n", rhti, nji, rmbid, zl_ll_count[rhti][rmbid][nji], sqrt(zl_ll_err2[rhti][rmbid][nji]) ) ;
                  r_numer += zl_ll_count[rhti][rmbin][nji] ;
                  r_denom += zl_ll_count[rhti][rmbid][nji] ;
                  r_numer_err2 += zl_ll_err2[rhti][rmbin][nji] ;
                  r_denom_err2 += zl_ll_err2[rhti][rmbid][nji] ;
               } // htbi
            } // nji
            float ratio(0.) ;
            float ratio_err(0.) ;
            if ( r_denom > 0 ) {
               ratio = r_numer / r_denom ;
               if ( r_numer > 0 ) {
                  ratio_err = ratio * sqrt( r_numer_err2 / (r_numer*r_numer) + r_denom_err2 / (r_denom*r_denom) ) ;
               }
            }
            printf("        %40s :  Ratio = %9.3f / %9.3f = %7.4f +/- %7.4f\n", rmht_name[ri], r_numer, r_denom, ratio, ratio_err ) ;
            fprintf( outfp, "%-35s   %7.4f  %7.4f\n", rmht_name[ri], ratio, ratio_err ) ;

         } // ri


         int n_rmht_map_lines = (n_rmht/n_search_njet) * n_search_nb * (nnjetbins-1) ;

         fprintf( outfp, "N-Rmht-map-lines   %d\n", n_rmht_map_lines ) ;

         for ( int nji=1; nji<nnjetbins; nji++ ) {
            for ( int ri=0; ri<(n_rmht/n_search_njet); ri++ ) {

               int rnji(-1), rhti(-1), rmbin(-1), rmbid(-1) ;
               sscanf( rmht_name[ri], "Rmht_sbnj%d_ht%d_fbmht%dover%d", &rnji, &rhti, &rmbin, &rmbid ) ;

               int sbnji(1) ;
               if ( nji>3 ) sbnji = nji-2 ;
               for ( int nbi=0; nbi<nnbbins; nbi++ ) {
                  fprintf( outfp, "FB-Njet%d-Nb%d-MHT%d-HT%d     Rmht_sbnj%d_ht%d_fbmht%dover%d   FB-Njet%d-Nb%d-MHT%d-HT%d\n",
                                      nji, nbi, rmbin, rhti,
                                      sbnji, rhti, rmbin, rmbid,
                                      nji, nbi, rmbid, rhti ) ;
               } // nbi

            } // ri
         } // nji

      } // include_mht_ratios ?



      fclose( outfp ) ;
      fclose( outfp_qcdcounts ) ;
      fclose( outfp_llcounts ) ;

   } // make_lhbuilder_input4b

  //==================================================================================================

   TH1F* get_hist( const char* hname ) {

      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) { printf("\n\n **** Missing histogram: %s\n\n", hname ) ; gSystem->Exit(-1) ; }

      return hp ;

   } // get_hist

  //==================================================================================================





