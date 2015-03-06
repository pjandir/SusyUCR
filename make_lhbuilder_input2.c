

#include "TString.h"
#include "TSystem.h"

#include "histio.c"
#include <stdio.h>

   TH1F* get_hist( const char* hname ) ;

  //---------------------------

   void make_lhbuilder_input2( bool perfect_closure = true,
                               float true_sig_strength = 1.0,
                            const char* infile = "outputfiles/fill-bg-hists2-t1bbbbH-postdraw.root",
                            const char* signame = "t1bbbbH",
                            const char* outfilebase = "outputfiles/lhbuilder-input"
                           ) {




      float qcd_kht_val[10] ;
      float qcd_kht_err[10] ;
      float qcd_kmht_val[10] ;
      float qcd_kmht_err[10] ;
      float qcd_knjet_val[10] ;
      float qcd_knjet_err[10] ;
      float qcd_knb_val[10] ;
      float qcd_knb_err[10] ;

      int n_qcd_kht_pars(3) ;
      qcd_kht_val[1] = 0.30 ;  qcd_kht_err[1] = 0.03 ;
      qcd_kht_val[2] = 0.23 ;  qcd_kht_err[2] = 0.03 ;
      qcd_kht_val[3] = 0.12 ;  qcd_kht_err[3] = 0.02 ;

      int n_qcd_kmht_pars(4) ;
      qcd_kmht_val[1] = 1    ;  qcd_kmht_err[1] = 0    ;
      qcd_kmht_val[2] = 1.2  ;  qcd_kmht_err[2] = 0.1  ;
      qcd_kmht_val[3] = 1.5  ;  qcd_kmht_err[3] = 0.1  ;
      qcd_kmht_val[4] = 1.6  ;  qcd_kmht_err[3] = 0.1  ;

      int n_qcd_knjet_pars(3) ;
      qcd_knjet_val[1] = 1    ;  qcd_knjet_err[1] = 0    ;
      qcd_knjet_val[2] = 0.66 ;  qcd_knjet_err[2] = 0.1  ;
      qcd_knjet_val[3] = 0.50 ;  qcd_knjet_err[3] = 0.1  ;

      int n_qcd_knb_pars(1) ;
      qcd_knb_val[1] = 1    ;  qcd_knb_err[1] = 0    ;



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
      fb_qcd_njet_par_ind[4] = 3 ;
      fb_qcd_njet_par_ind[5] = 3 ;

      fb_qcd_nb_par_ind[0] = 1 ;
      fb_qcd_nb_par_ind[1] = 1 ;
      fb_qcd_nb_par_ind[2] = 1 ;
      fb_qcd_nb_par_ind[3] = 1 ;




      char hname[1000] ;
      TH1F* hp ;

      loadHist( infile ) ;

      char outfile[10000] ;
      sprintf( outfile, "%s-%s.txt", outfilebase, signame ) ;
      FILE* outfp ;
      printf("\n\n Creating output file: %s\n\n", outfile ) ;
      if ( (outfp=fopen( outfile, "w" ))==NULL ) {
         printf("\n\n *** can't open output file.\n\n") ; return ;
      }


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


      hp = get_hist( "h_mhtvsht_nb0_nj1_zl_ttbar_1d" ) ;
      int nhbins = hp -> GetNbinsX() ;





      for ( int nji=1; nji<nnjetbins; nji++ ) {

         for ( int nbi=0; nbi<nnbbins; nbi++ ) {

            for ( int hbi=1; hbi<nhbins; hbi++ ) {

               bool emptyBin(false) ;


               sprintf( hname, "h_mhtvsht_nb%d_nj%d_ldp_qcd_1d", nbi, nji ) ;
               TH1F* hp_qcd_ldp = get_hist( hname ) ;
               sprintf( hname, "h_mhtvsht_nb%d_nj%d_zl_qcd_1d", nbi, nji ) ;
               TH1F* hp_qcd_zl = get_hist( hname ) ;
               sprintf( hname, "h_mhtvsht_nb%d_nj%d_ldp_%s_1d", nbi, nji, signame ) ;
               TH1F* hp_sig_ldp = get_hist( hname ) ;
               sprintf( hname, "h_mhtvsht_nb%d_nj%d_zl_%s_1d", nbi, nji, signame ) ;
               TH1F* hp_sig_zl = get_hist( hname ) ;
               sprintf( hname, "h_mhtvsht_nb%d_nj%d_sl_%s_1d", nbi, nji, signame ) ;
               TH1F* hp_sig_sl = get_hist( hname ) ;



               sprintf( hname, "h_bgsum_nb%d_nj%d_ldp", nbi, nji ) ;
               TH1F* hp_bgsum_ldp = get_hist( hname ) ;
               sprintf( hname, "h_bgsum_nb%d_nj%d_zl", nbi, nji ) ;
               TH1F* hp_bgsum_zl = get_hist( hname ) ;
               sprintf( hname, "h_bgsum_nb%d_nj%d_sl", nbi, nji ) ;
               TH1F* hp_bgsum_sl = get_hist( hname ) ;



               sprintf( hname, "h_llbgsum_nb%d_nj%d_ldp", nbi, nji ) ;
               TH1F* hp_llbgsum_ldp = get_hist( hname ) ;
               sprintf( hname, "h_llbgsum_nb%d_nj%d_zl", nbi, nji ) ;
               TH1F* hp_llbgsum_zl = get_hist( hname ) ;
               sprintf( hname, "h_llbgsum_nb%d_nj%d_sl", nbi, nji ) ;
               TH1F* hp_llbgsum_sl = get_hist( hname ) ;




               char binlabel[100] ;

               sprintf( binlabel, "%s", hp_llbgsum_zl -> GetXaxis() -> GetBinLabel( hbi ) ) ;
               if ( strlen( binlabel ) == 0 ) { emptyBin = true ; continue ; }
               int htbin(-1), mhtbin(-1) ;
               sscanf( binlabel, "MHT%d_HT%d", &mhtbin, &htbin ) ;
               //////////// printf("      hist bin label %s : mht bin = %d ,  ht bin = %d\n", binlabel, mhtbin, htbin ) ;
               TString blts( binlabel ) ;
               blts.ReplaceAll("_","-") ;




               float nsl_val = hp_bgsum_sl -> GetBinContent( hbi ) ;
               float nsl_err = hp_bgsum_sl -> GetBinError( hbi ) ;

               float nldp_val = hp_bgsum_ldp -> GetBinContent( hbi ) ;
               float nldp_err = hp_bgsum_ldp -> GetBinError( hbi ) ;

               float nzl_val = hp_bgsum_zl -> GetBinContent( hbi ) ;
               float nzl_err = hp_bgsum_zl -> GetBinError( hbi ) ;

               float nsig_sl_val = hp_sig_sl -> GetBinContent( hbi ) ;
               float nsig_sl_err = hp_sig_sl -> GetBinError( hbi ) ;

               float nsig_ldp_val = hp_sig_ldp -> GetBinContent( hbi ) ;
               float nsig_ldp_err = hp_sig_ldp -> GetBinError( hbi ) ;

               float nsig_zl_val = hp_sig_zl -> GetBinContent( hbi ) ;
               float nsig_zl_err = hp_sig_zl -> GetBinError( hbi ) ;


               float nll_sl_val = hp_llbgsum_sl -> GetBinContent( hbi ) ;
               float nll_sl_err = hp_llbgsum_sl -> GetBinError( hbi ) ;

               float nll_zl_val = hp_llbgsum_zl -> GetBinContent( hbi ) ;
               float nll_zl_err = hp_llbgsum_zl -> GetBinError( hbi ) ;

               float nll_ldp_val = hp_llbgsum_ldp -> GetBinContent( hbi ) ;
               float nll_ldp_err = hp_llbgsum_ldp -> GetBinError( hbi ) ;


               float nqcd_zl_val = hp_qcd_zl -> GetBinContent( hbi ) ;
               float nqcd_zl_err = hp_qcd_zl -> GetBinError( hbi ) ;

               float nqcd_ldp_val = hp_qcd_ldp -> GetBinContent( hbi ) ;
               float nqcd_ldp_err = hp_qcd_ldp -> GetBinError( hbi ) ;



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

               float calc_nzl_bg = nll_zl_val + calc_nqcd_zl ;
               float mc_nzl_bg   = nll_zl_val + nqcd_zl_val ;

               printf("   Nzl bg  calc = %7.1f,  mc = %7.1f\n",  calc_nzl_bg, mc_nzl_bg ) ;

               float nzl_output ;
               if ( perfect_closure ) {
                  nzl_output = calc_nzl_bg + true_sig_strength * nsig_zl_val ;
               } else {
                  nzl_output =   mc_nzl_bg + true_sig_strength * nsig_zl_val ;
               }

               float nsl_output = nll_sl_val + true_sig_strength * nsig_sl_val ;
               float nldp_output = nll_ldp_val + nqcd_ldp_val + true_sig_strength * nsig_ldp_val ;


               float rsl_zl_val = 0. ;
               float rsl_zl_err = 0. ;
               if ( nll_sl_val > 0 ) {
                  rsl_zl_val = nll_zl_val / nll_sl_val ;
                  if ( nll_zl_val > 0 ) {
                     rsl_zl_err = rsl_zl_val * sqrt( pow( nll_sl_err/nll_sl_val, 2 ) + pow( nll_zl_err/nll_zl_val, 2 ) ) ;
                  }
               }

               float rsl_ldp_val = 0. ;
               float rsl_ldp_err = 0. ;
               if ( nll_sl_val > 0 ) {
                  rsl_ldp_val = nll_ldp_val / nll_sl_val ;
                  if ( nll_ldp_val > 0 ) {
                     rsl_ldp_err = rsl_ldp_val * sqrt( pow( nll_sl_err/nll_sl_val, 2 ) + pow( nll_ldp_err/nll_ldp_val, 2 ) ) ;
                  }
               }

               float rqcd_val = 0. ;
               float rqcd_err = 0. ;
               if ( nqcd_ldp_val > 0 ) {
                  rqcd_val = nqcd_zl_val / nqcd_ldp_val ;
                  if ( nqcd_zl_val > 0 ) {
                     rqcd_err = rqcd_val * sqrt( pow( nqcd_zl_err/nqcd_zl_val, 2 ) + pow( nqcd_ldp_err/nqcd_ldp_val, 2 ) ) ;
                  }
               }




               printf(         "FB-Njet%d-Nb%d-%s  ", nji, nbi, blts.Data() ) ;
               fprintf( outfp, "FB-Njet%d-Nb%d-%s  ", nji, nbi, blts.Data() ) ;

               ////////////printf(         " %6.3f %6.3f    %6.3f %6.3f ", rsl_zl_val, rsl_zl_err,  rsl_ldp_val, rsl_ldp_err ) ;
               ////////////fprintf( outfp, " %6.3f %6.3f    %6.3f %6.3f ", rsl_zl_val, rsl_zl_err,  rsl_ldp_val, rsl_ldp_err ) ;

               printf(         " %6.3f %6.3f    %6.3f %6.3f ", rsl_zl_val, 0.,  rsl_ldp_val, 0. ) ;
               fprintf( outfp, " %6.3f %6.3f    %6.3f %6.3f ", rsl_zl_val, 0.,  rsl_ldp_val, 0. ) ;

               ////////////////// printf( "   %7.1f  %7.1f  %7.1f  ", nzl_val, nsl_val, nldp_val ) ;
               ////////////////// fprintf( outfp, "   %7.1f  %7.1f  %7.1f  ", nzl_val, nsl_val, nldp_val ) ;

               printf(         "   %9.3f  %9.3f  %9.3f  ", nzl_output, nsl_output, nldp_output ) ;
               fprintf( outfp, "   %9.3f  %9.3f  %9.3f  ", nzl_output, nsl_output, nldp_output ) ;

               printf(         "   %9.3f  %9.3f  %9.3f", nsig_zl_val, nsig_sl_val, nsig_ldp_val ) ;
               fprintf( outfp, "   %9.3f  %9.3f  %9.3f", nsig_zl_val, nsig_sl_val, nsig_ldp_val ) ;


               printf("\n") ;
               fprintf( outfp, "\n") ;


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
               for ( int hbi=1; hbi<nhtbins; hbi++ ) {
                  printf( "FB-Njet%d-Nb%d-MHT%d-HT%d  ", nji, nbi, mbi, hbi ) ;
                  fprintf( outfp, "FB-Njet%d-Nb%d-MHT%d-HT%d  ", nji, nbi, mbi, hbi ) ;
                  int sb_nji(-1) ;
                  int sb_nbi(-1) ;
                  int sb_mhthtbi(-1) ;
                  for ( int i=0; i<n_search_njet; i++ ) { if (  njetbins[nji]>= search_njet_bins[i] && njetbins[nji+1]<= search_njet_bins[i+1] ) sb_nji = i ; }
                  for ( int i=0; i<n_search_nb; i++ ) { if (  nbbins[nbi]>= search_nb_bins[i] && nbbins[nbi+1]<= search_nb_bins[i+1] ) sb_nbi = i ; }
                  for ( int i=0; i<n_search_mhtht; i++ ) {
                     if ( mhtbins[mbi]>= search_mhtht_mhtlow[i] && mhtbins[mbi+1] <= search_mhtht_mhthigh[i]
                       && htbins[hbi]>= search_mhtht_htlow[i] && htbins[hbi+1] <= search_mhtht_hthigh[i] ) sb_mhthtbi = i ;
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

               } // hbi
            } // mbi
         } // nbi
      } // nji


      fprintf( outfp, "N-QCD-Kht-pars    %d\n", n_qcd_kht_pars ) ;
      for ( int i=1; i<=n_qcd_kht_pars; i++ ) {
         fprintf( outfp, "QCD-Kht%d    %5.2f  %5.2f\n", i, qcd_kht_val[i], qcd_kht_err[i] ) ;
         //fprintf( outfp, "QCD-Kht%d    %5.2f  %5.2f\n", i, qcd_kht_val[i], 0. ) ;
      }
      fprintf( outfp, "N-QCD-Kmht-pars    %d\n", n_qcd_kmht_pars ) ;
      for ( int i=1; i<=n_qcd_kmht_pars; i++ ) {
         fprintf( outfp, "QCD-Kmht%d    %5.2f  %5.2f\n", i, qcd_kmht_val[i], qcd_kmht_err[i] ) ;
         //fprintf( outfp, "QCD-Kmht%d    %5.2f  %5.2f\n", i, qcd_kmht_val[i], 0. ) ;
      }
      fprintf( outfp, "N-QCD-Knjet-pars    %d\n", n_qcd_knjet_pars ) ;
      for ( int i=1; i<=n_qcd_knjet_pars; i++ ) {
         fprintf( outfp, "QCD-Knjet%d    %5.2f  %5.2f\n", i, qcd_knjet_val[i], qcd_knjet_err[i] ) ;
         //fprintf( outfp, "QCD-Knjet%d    %5.2f  %5.2f\n", i, qcd_knjet_val[i], 0. ) ;
      }
      fprintf( outfp, "N-QCD-Knb-pars    %d\n", n_qcd_knb_pars ) ;
      for ( int i=1; i<=n_qcd_knb_pars; i++ ) {
         fprintf( outfp, "QCD-Knb%d    %5.2f  %5.2f\n", i, qcd_knb_val[i], qcd_knb_err[i] ) ;
         //fprintf( outfp, "QCD-Knb%d    %5.2f  %5.2f\n", i, qcd_knb_val[i], 0. ) ;
      }


      fprintf( outfp, "FB-Njet1-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[1] ) ;
      fprintf( outfp, "FB-Njet2-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[2] ) ;
      fprintf( outfp, "FB-Njet3-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[3] ) ;
      fprintf( outfp, "FB-Njet4-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[4] ) ;
      fprintf( outfp, "FB-Njet5-QCD-Njet-par-ind    %d\n", fb_qcd_njet_par_ind[5] ) ;

      fprintf( outfp, "FB-Nb0-QCD-Nb-par-ind    %d\n", fb_qcd_nb_par_ind[0] ) ;
      fprintf( outfp, "FB-Nb1-QCD-Nb-par-ind    %d\n", fb_qcd_nb_par_ind[1] ) ;
      fprintf( outfp, "FB-Nb2-QCD-Nb-par-ind    %d\n", fb_qcd_nb_par_ind[2] ) ;
      fprintf( outfp, "FB-Nb3-QCD-Nb-par-ind    %d\n", fb_qcd_nb_par_ind[3] ) ;

      fprintf( outfp, "FB-MHT1-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[1] ) ;
      fprintf( outfp, "FB-MHT2-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[2] ) ;
      fprintf( outfp, "FB-MHT3-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[3] ) ;
      fprintf( outfp, "FB-MHT4-QCD-MHT-par-ind    %d\n", fb_qcd_mht_par_ind[4] ) ;

      fprintf( outfp, "FB-HT1-QCD-HT-par-ind    %d\n", fb_qcd_ht_par_ind[1] ) ;
      fprintf( outfp, "FB-HT2-QCD-HT-par-ind    %d\n", fb_qcd_ht_par_ind[2] ) ;
      fprintf( outfp, "FB-HT3-QCD-HT-par-ind    %d\n", fb_qcd_ht_par_ind[3] ) ;

      fclose( outfp ) ;

   } // make_lhbuilder_input2

  //==================================================================================================

   TH1F* get_hist( const char* hname ) {

      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) { printf("\n\n **** Missing histogram: %s\n\n", hname ) ; gSystem->Exit(-1) ; }

      return hp ;

   } // get_hist

  //==================================================================================================





