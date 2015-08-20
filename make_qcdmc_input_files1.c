
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

#include <fstream>

   int fbhi[300] ;
   int fbmi[300] ;
   int fbji[300] ;
   int fbbi[300] ;
   int fbmhtiti[300] ;
   int nfb ;


   int fb_mhti_from_mhthti( int mhthti ) ;
   int fb_hti_from_mhthti( int mhthti ) ;

   void fill_fb_arrays( int sb_mhtht_bi, int sb_njets_bi, int sb_nbjets_bi ) ;

  //-------

   void make_qcdmc_input_files1( const char* infile = "outputfiles/fill-bg-hists-tm4-allinone-postdraw.root",
                                  const char* outfile_kqcd_fit = "outputfiles/kqcd-input-qcdmc.txt",
                                  const char* outfile_combine = "outputfiles/combine-input-qcdmc.txt",
                                  const char* outfile_finebins = "outputfiles/finebin-input-qcdmc.txt"
                                  ) {

      TFile* tf_input = new TFile( infile, "READ" ) ;
      if ( tf_input == 0x0 ) { printf("\n\n *** Bad input file: %s\n\n", infile ) ; return ; }

      TH1* h_lowdphi = (TH1*) tf_input -> Get( "h_1d_allinone_ldp_qcd" ) ;
      if ( h_lowdphi == 0x0 ) { printf("\n\n *** Missing hist: h_1d_allinone_ldp_qcd\n\n" ) ; return ; }

      TH1* h_highdphi = (TH1*) tf_input -> Get( "h_1d_allinone_zl_qcd" ) ;
      if ( h_highdphi == 0x0 ) { printf("\n\n *** Missing hist: h_1d_allinone_zl_qcd\n\n" ) ; return ; }



      FILE* ofp_kqcd_fit ;
      if ( (ofp_kqcd_fit=fopen( outfile_kqcd_fit, "w" ))==NULL ) {
         printf("\n\n *** Problem opening output file %s\n\n", outfile_kqcd_fit ) ;
      }

      TH1F* h_kqcd_input_qcdmc_lowdphi = new TH1F( "h_kqcd_input_qcdmc_lowdphi", "Kqcd fit input, QCD MC, low DeltaPhi", 55, 0.5, 55.5 ) ;
      TH1F* h_kqcd_input_qcdmc_highdphi = new TH1F( "h_kqcd_input_qcdmc_highdphi", "Kqcd fit input, QCD MC, high DeltaPhi", 55, 0.5, 55.5 ) ;

      TH1F* h_kqcd_input_qcdmc_highlow_ratio = new TH1F( "h_kqcd_input_qcdmc_highlow_ratio", "Kqcd fit input, QCD MC, high/low ratio", 55, 0.5, 55.5 ) ;

      for ( int nji=1; nji<=5; nji++ ) {
         for ( int mhthti=1; mhthti<=11; mhthti++ ) {

            float nbsum_lowdphi_val(0.) ;
            float nbsum_lowdphi_err2(0.) ;

            float nbsum_highdphi_val(0.) ;
            float nbsum_highdphi_err2(0.) ;

            for ( int nbi=0; nbi<=3; nbi++ ) {

               int owen_bi = 44*(nji-1) + nbi*11 + mhthti ;

               nbsum_lowdphi_val += h_lowdphi -> GetBinContent( owen_bi ) ;
               nbsum_lowdphi_err2 += pow( ( h_lowdphi -> GetBinError( owen_bi ) ), 2. ) ;

               nbsum_highdphi_val += h_highdphi -> GetBinContent( owen_bi ) ;
               nbsum_highdphi_err2 += pow( ( h_highdphi -> GetBinError( owen_bi ) ), 2. ) ;

            } // nbi

            float ratio_val(0.) ;
            float ratio_err(0.) ;
            if ( nbsum_lowdphi_val > 0 ) {
               ratio_val = nbsum_highdphi_val / nbsum_lowdphi_val ;
               if ( nbsum_highdphi_val > 0 ) {
                  ratio_err = ratio_val * sqrt( nbsum_lowdphi_err2 / (nbsum_lowdphi_val*nbsum_lowdphi_val)
                                              + nbsum_highdphi_err2 / (nbsum_highdphi_val*nbsum_highdphi_val) ) ;
               } else {
                  ratio_err = ratio_val * sqrt( nbsum_lowdphi_err2 / (nbsum_lowdphi_val*nbsum_lowdphi_val) ) ;
               }
            }

            int bi = 11*(nji-1) + mhthti ;

            int mhti = fb_mhti_from_mhthti( mhthti ) ;
            int hti  = fb_hti_from_mhthti( mhthti ) ;

            char binlabel[100] ;
            sprintf( binlabel, "FB-Njet%d-Nbsum-MHT%d-HT%d  %3d", nji, mhti, hti, bi ) ;

            h_kqcd_input_qcdmc_lowdphi -> SetBinContent( bi, nbsum_lowdphi_val ) ;
            h_kqcd_input_qcdmc_lowdphi -> SetBinError( bi, sqrt(nbsum_lowdphi_err2) ) ;
            h_kqcd_input_qcdmc_lowdphi -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;

            h_kqcd_input_qcdmc_highdphi -> SetBinContent( bi, nbsum_highdphi_val ) ;
            h_kqcd_input_qcdmc_highdphi -> SetBinError( bi, sqrt(nbsum_highdphi_err2) ) ;
            h_kqcd_input_qcdmc_highdphi -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;

            h_kqcd_input_qcdmc_highlow_ratio -> SetBinContent( bi, ratio_val ) ;
            h_kqcd_input_qcdmc_highlow_ratio -> SetBinError( bi, ratio_err ) ;
            h_kqcd_input_qcdmc_highlow_ratio -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;

            printf( " %3d FB-Njet%d-Nbsum-MHT%d-HT%d    %7.1f +/- %5.1f     %7.1f +/- %5.1f\n", bi, nji, mhti, hti,
               nbsum_lowdphi_val, sqrt( nbsum_lowdphi_err2 ),
               nbsum_highdphi_val, sqrt( nbsum_highdphi_err2 )
               ) ;

            fprintf( ofp_kqcd_fit, " %3d FB-Njet%d-Nbsum-MHT%d-HT%d    %7.1f +/- %5.1f     %7.1f +/- %5.1f\n", bi, nji, mhti, hti,
               nbsum_lowdphi_val, sqrt( nbsum_lowdphi_err2 ),
               nbsum_highdphi_val, sqrt( nbsum_highdphi_err2 )
               ) ;


         } // mhthti
      } // nji

      fclose( ofp_kqcd_fit ) ;

      h_kqcd_input_qcdmc_lowdphi -> GetXaxis() -> LabelsOption( "v" ) ;
      h_kqcd_input_qcdmc_highdphi -> GetXaxis() -> LabelsOption( "v" ) ;

      h_kqcd_input_qcdmc_lowdphi -> SetMinimum(0.1) ;
      h_kqcd_input_qcdmc_highdphi -> SetMinimum(0.1) ;

      h_kqcd_input_qcdmc_lowdphi -> SetFillColor( kRed-9 ) ;
      h_kqcd_input_qcdmc_highdphi -> SetFillColor( kRed-9 ) ;

      gStyle -> SetPadBottomMargin( 0.38 ) ;
      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      TCanvas* can_kqcd = new TCanvas( "can_kqcd", "Kqcd fit inputs, QCD MC", 1100, 900 ) ;
      can_kqcd -> Divide(1,2) ;

      char pdffile[10000] ;

    //---
      can_kqcd -> cd(1) ;
      h_kqcd_input_qcdmc_lowdphi -> Draw() ;
      h_kqcd_input_qcdmc_lowdphi -> Draw( "hist same" ) ;
      h_kqcd_input_qcdmc_lowdphi -> Draw( "same" ) ;
      h_kqcd_input_qcdmc_lowdphi -> Draw( "axis same" ) ;
      h_kqcd_input_qcdmc_lowdphi -> Draw( "axig same" ) ;

      can_kqcd -> cd(2) ;
      h_kqcd_input_qcdmc_highdphi -> Draw() ;
      h_kqcd_input_qcdmc_highdphi -> Draw( "hist same" ) ;
      h_kqcd_input_qcdmc_highdphi -> Draw( "same" ) ;
      h_kqcd_input_qcdmc_highdphi -> Draw( "axis same" ) ;
      h_kqcd_input_qcdmc_highdphi -> Draw( "axig same" ) ;

      sprintf( pdffile, "outputfiles/kqcd-input-qcdmc-liny.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;

    //---
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(1) ;
      can_kqcd -> cd(2) ;
      gPad -> SetLogy(1) ;

      sprintf( pdffile, "outputfiles/kqcd-input-qcdmc-logy.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;


    //---
      h_kqcd_input_qcdmc_lowdphi -> SetMaximum(20) ;
      h_kqcd_input_qcdmc_highdphi -> SetMaximum(20) ;
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/kqcd-input-qcdmc-zoom1.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;


    //---
      h_kqcd_input_qcdmc_lowdphi -> SetMaximum(5) ;
      h_kqcd_input_qcdmc_highdphi -> SetMaximum(5) ;
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/kqcd-input-qcdmc-zoom2.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;



   //---------------------------------------------

      TCanvas* can_kqcd_ratio = new TCanvas( "can_kqcd_ratio", "Kqcd fit inputs, QCD MC, high/low ratio", 1100, 900 ) ;

      h_kqcd_input_qcdmc_highlow_ratio -> SetMarkerStyle(20) ;
      h_kqcd_input_qcdmc_highlow_ratio -> SetMaximum(0.55) ;
      h_kqcd_input_qcdmc_highlow_ratio -> SetMinimum(-0.15) ;

      h_kqcd_input_qcdmc_highlow_ratio -> Draw() ;

      sprintf( pdffile, "outputfiles/kqcd-input-qcdmc-highlow-ratio.pdf" ) ;
      can_kqcd_ratio -> Update() ; can_kqcd_ratio -> Draw() ; can_kqcd_ratio -> SaveAs( pdffile ) ;



   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FILE* ofp_combine ;
      if ( (ofp_combine=fopen( outfile_combine, "w" ))==NULL ) {
         printf("\n\n *** Problem opening output file %s\n\n", outfile_combine ) ;
      }

      TH1F* h_combine_input_qcdmc_lowdphi = new TH1F( "h_combine_input_qcdmc_lowdphi", "Combine input, QCD MC, low DeltaPhi", 72, 0.5, 72.5 ) ;
      TH1F* h_combine_input_qcdmc_highdphi = new TH1F( "h_combine_input_qcdmc_highdphi", "Combine input, QCD MC, high DeltaPhi", 72, 0.5, 72.5 ) ;

      for ( int sb_nji=1; sb_nji<=3; sb_nji++ ) {
         for ( int sb_nbi=0; sb_nbi<=3; sb_nbi++ ) {
            for ( int sb_mhthti=1; sb_mhthti<=6; sb_mhthti++ ) {

               int owen_sbi = 24*(sb_nji-1) + 6*sb_nbi + sb_mhthti ;
               char sb_mhthti_string[100] ;
               if ( sb_mhthti==1 ) sprintf( sb_mhthti_string, "MHT1-HT1" ) ;
               if ( sb_mhthti==2 ) sprintf( sb_mhthti_string, "MHT1-HT2" ) ;
               if ( sb_mhthti==3 ) sprintf( sb_mhthti_string, "MHT1-HT3" ) ;
               if ( sb_mhthti==4 ) sprintf( sb_mhthti_string, "MHT2-HT12" ) ;
               if ( sb_mhthti==5 ) sprintf( sb_mhthti_string, "MHT2-HT3" ) ;
               if ( sb_mhthti==6 ) sprintf( sb_mhthti_string, "MHT3-HT23" ) ;

               nfb = 0 ;
               fill_fb_arrays( sb_mhthti, sb_nji, sb_nbi ) ;

               float fbsum_lowdphi_val(0.) ;
               float fbsum_lowdphi_err2(0.) ;

               float fbsum_highdphi_val(0.) ;
               float fbsum_highdphi_err2(0.) ;

               for ( int fbi=0; fbi<nfb; fbi++ ) {


                  int owen_fbi = 44*(fbji[fbi]-1) + fbbi[fbi]*11 + fbmhtiti[fbi] ;

                  //printf("  SB %3d,  SB-Njet%d-Nb%d-%-9s  |  owen FB %3d :  FB-Njet%d-Nb%d-MHT%d-HT%d\n",
                  //  owen_sbi, sb_nji, sb_nbi, sb_mhthti_string,
                  //  owen_fbi, fbji[fbi], fbbi[fbi], fbmi[fbi], fbhi[fbi] ) ;


                  fbsum_lowdphi_val += h_lowdphi -> GetBinContent( owen_fbi ) ;
                  fbsum_lowdphi_err2 += pow( (h_lowdphi -> GetBinError( owen_fbi )), 2. ) ;

                  fbsum_highdphi_val += h_highdphi -> GetBinContent( owen_fbi ) ;
                  fbsum_highdphi_err2 += pow( (h_highdphi -> GetBinError( owen_fbi )), 2. ) ;


               } // fbi

               char binlabel[100] ;
               sprintf( binlabel, "SB-Njet%d-Nb%d-%-9s  %3d", sb_nji, sb_nbi, sb_mhthti_string, owen_sbi ) ;

               h_combine_input_qcdmc_lowdphi -> SetBinContent( owen_sbi, fbsum_lowdphi_val ) ;
               h_combine_input_qcdmc_lowdphi -> SetBinError( owen_sbi, sqrt(fbsum_lowdphi_err2) ) ;

               h_combine_input_qcdmc_highdphi -> SetBinContent( owen_sbi, fbsum_highdphi_val ) ;
               h_combine_input_qcdmc_highdphi -> SetBinError( owen_sbi, sqrt(fbsum_highdphi_err2) ) ;

               h_combine_input_qcdmc_lowdphi  -> GetXaxis() -> SetBinLabel( owen_sbi, binlabel ) ;
               h_combine_input_qcdmc_highdphi -> GetXaxis() -> SetBinLabel( owen_sbi, binlabel ) ;


               printf( " %3d SB-Njet%d-Nb%d-%-9s    %7.1f +/- %5.1f     %7.1f +/- %5.1f\n",
                  owen_sbi, sb_nji, sb_nbi, sb_mhthti_string,
                  fbsum_lowdphi_val, sqrt( fbsum_lowdphi_err2 ),
                  fbsum_highdphi_val, sqrt( fbsum_highdphi_err2 )
                  ) ;

               fprintf( ofp_combine, " %3d SB-Njet%d-Nb%d-%-9s    %7.1f +/- %5.1f     %7.1f +/- %5.1f\n",
                  owen_sbi, sb_nji, sb_nbi, sb_mhthti_string,
                  fbsum_lowdphi_val, sqrt( fbsum_lowdphi_err2 ),
                  fbsum_highdphi_val, sqrt( fbsum_highdphi_err2 )
                  ) ;



            } // sb_mhthti
         } // sb_nbi
      } // sb_nji

      fclose( ofp_combine ) ;

      h_combine_input_qcdmc_lowdphi  -> GetXaxis() -> LabelsOption( "v" ) ;
      h_combine_input_qcdmc_highdphi -> GetXaxis() -> LabelsOption( "v" ) ;

      h_combine_input_qcdmc_lowdphi -> SetMinimum(0.1) ;
      h_combine_input_qcdmc_highdphi -> SetMinimum(0.1) ;

      h_combine_input_qcdmc_lowdphi -> SetFillColor( kRed-9 ) ;
      h_combine_input_qcdmc_highdphi -> SetFillColor( kRed-9 ) ;

      gStyle -> SetPadBottomMargin( 0.35 ) ;
      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      TCanvas* can_combine = new TCanvas( "can_combine", "combine fit inputs, QCD MC", 1100, 900 ) ;
      can_combine -> Divide(1,2) ;


    //---
      can_combine -> cd(1) ;
      h_combine_input_qcdmc_lowdphi -> Draw() ;
      h_combine_input_qcdmc_lowdphi -> Draw( "hist same" ) ;
      h_combine_input_qcdmc_lowdphi -> Draw( "same" ) ;
      h_combine_input_qcdmc_lowdphi -> Draw( "axis same" ) ;
      h_combine_input_qcdmc_lowdphi -> Draw( "axig same" ) ;

      can_combine -> cd(2) ;
      h_combine_input_qcdmc_highdphi -> Draw() ;
      h_combine_input_qcdmc_highdphi -> Draw( "hist same" ) ;
      h_combine_input_qcdmc_highdphi -> Draw( "same" ) ;
      h_combine_input_qcdmc_highdphi -> Draw( "axis same" ) ;
      h_combine_input_qcdmc_highdphi -> Draw( "axig same" ) ;

      sprintf( pdffile, "outputfiles/combine-input-qcdmc-liny.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;

    //---
      can_combine -> cd(1) ;
      gPad -> SetLogy(1) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(1) ;

      sprintf( pdffile, "outputfiles/combine-input-qcdmc-logy.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


    //---
      h_combine_input_qcdmc_lowdphi -> SetMaximum(20) ;
      h_combine_input_qcdmc_highdphi -> SetMaximum(20) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/combine-input-qcdmc-zoom1.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


    //---
      h_combine_input_qcdmc_lowdphi -> SetMaximum(5) ;
      h_combine_input_qcdmc_highdphi -> SetMaximum(5) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/combine-input-qcdmc-zoom2.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FILE* ofp_finebins ;
      if ( (ofp_finebins=fopen( outfile_finebins, "w" ))==NULL ) {
         printf("\n\n *** Problem opening output file %s\n\n", outfile_finebins ) ;
      }

      TH1F* h_finebin_input_qcdmc_lowdphi = (TH1F*)  h_lowdphi -> Clone( "h_finebin_input_qcdmc_lowdphi" ) ;
      TH1F* h_finebin_input_qcdmc_highdphi = (TH1F*)  h_highdphi -> Clone( "h_finebin_input_qcdmc_highdphi" ) ;

      for ( int fb_nji=1; fb_nji<=5; fb_nji++ ) {
         for ( int fb_nbi=0; fb_nbi<=3; fb_nbi++ ) {
            for ( int fb_mhthti=1; fb_mhthti<=11; fb_mhthti++ ) {


               int owen_fbi = 44*(fb_nji-1) + fb_nbi*11 + fb_mhthti ;

               char  mhtht_bin_label[100] ;
               if ( fb_mhthti ==  1 ) { sprintf( mhtht_bin_label, "MHT1-HT1" ) ; }
               if ( fb_mhthti ==  2 ) { sprintf( mhtht_bin_label, "MHT1-HT2" ) ; }
               if ( fb_mhthti ==  3 ) { sprintf( mhtht_bin_label, "MHT1-HT3" ) ; }
               if ( fb_mhthti ==  4 ) { sprintf( mhtht_bin_label, "MHT2-HT1" ) ; }
               if ( fb_mhthti ==  5 ) { sprintf( mhtht_bin_label, "MHT2-HT2" ) ; }
               if ( fb_mhthti ==  6 ) { sprintf( mhtht_bin_label, "MHT2-HT3" ) ; }
               if ( fb_mhthti ==  7 ) { sprintf( mhtht_bin_label, "MHT3-HT1" ) ; }
               if ( fb_mhthti ==  8 ) { sprintf( mhtht_bin_label, "MHT3-HT2" ) ; }
               if ( fb_mhthti ==  9 ) { sprintf( mhtht_bin_label, "MHT3-HT3" ) ; }
               if ( fb_mhthti == 10 ) { sprintf( mhtht_bin_label, "MHT4-HT2" ) ; }
               if ( fb_mhthti == 11 ) { sprintf( mhtht_bin_label, "MHT4-HT3" ) ; }

               char binlabel[100] ;
               sprintf( binlabel, "FB-Njet%d-Nb%d-%s  %3d", fb_nji, fb_nbi, mhtht_bin_label, owen_fbi ) ;

               float lowdphi_val = h_finebin_input_qcdmc_lowdphi -> GetBinContent( owen_fbi ) ;
               float lowdphi_err = h_finebin_input_qcdmc_lowdphi -> GetBinError( owen_fbi ) ;

               float highdphi_val = h_finebin_input_qcdmc_highdphi -> GetBinContent( owen_fbi ) ;
               float highdphi_err = h_finebin_input_qcdmc_highdphi -> GetBinError( owen_fbi ) ;

               h_finebin_input_qcdmc_lowdphi  -> GetXaxis() -> SetBinLabel( owen_fbi, binlabel ) ;
               h_finebin_input_qcdmc_highdphi -> GetXaxis() -> SetBinLabel( owen_fbi, binlabel ) ;

               printf( " %3d FB-Njet%d-Nb%d-%s    %7.1f +/- %5.1f     %7.1f +/- %5.1f\n",
                  owen_fbi, fb_nji, fb_nbi, mhtht_bin_label,
                  lowdphi_val, lowdphi_err,
                  highdphi_val, highdphi_err
                  ) ;

               fprintf( ofp_finebins, " %3d FB-Njet%d-Nb%d-%s    %7.1f +/- %5.1f     %7.1f +/- %5.1f\n",
                  owen_fbi, fb_nji, fb_nbi, mhtht_bin_label,
                  lowdphi_val, lowdphi_err,
                  highdphi_val, highdphi_err
                  ) ;



            } // fb_mhthti
         } // fb_nbi
      } // fb_nji

      fclose( ofp_finebins ) ;

      h_finebin_input_qcdmc_lowdphi  -> GetXaxis() -> LabelsOption( "v" ) ;
      h_finebin_input_qcdmc_highdphi -> GetXaxis() -> LabelsOption( "v" ) ;

      h_finebin_input_qcdmc_lowdphi -> SetMinimum(0.1) ;
      h_finebin_input_qcdmc_highdphi -> SetMinimum(0.1) ;

      h_finebin_input_qcdmc_lowdphi -> SetFillColor( kRed-9 ) ;
      h_finebin_input_qcdmc_highdphi -> SetFillColor( kRed-9 ) ;

      gStyle -> SetPadBottomMargin( 0.35 ) ;
      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      can_combine -> Clear() ;
      can_combine -> Divide(1,2) ;


    //---
      can_combine -> cd(1) ;
      h_finebin_input_qcdmc_lowdphi -> Draw() ;
      h_finebin_input_qcdmc_lowdphi -> Draw( "hist same" ) ;
      h_finebin_input_qcdmc_lowdphi -> Draw( "same" ) ;
      h_finebin_input_qcdmc_lowdphi -> Draw( "axis same" ) ;
      h_finebin_input_qcdmc_lowdphi -> Draw( "axig same" ) ;

      can_combine -> cd(2) ;
      h_finebin_input_qcdmc_highdphi -> Draw() ;
      h_finebin_input_qcdmc_highdphi -> Draw( "hist same" ) ;
      h_finebin_input_qcdmc_highdphi -> Draw( "same" ) ;
      h_finebin_input_qcdmc_highdphi -> Draw( "axis same" ) ;
      h_finebin_input_qcdmc_highdphi -> Draw( "axig same" ) ;

      sprintf( pdffile, "outputfiles/finebin-input-qcdmc-liny.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;

    //---
      can_combine -> cd(1) ;
      gPad -> SetLogy(1) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(1) ;

      sprintf( pdffile, "outputfiles/finebin-input-qcdmc-logy.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


    //---
      h_finebin_input_qcdmc_lowdphi -> SetMaximum(20) ;
      h_finebin_input_qcdmc_highdphi -> SetMaximum(20) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/finebin-input-qcdmc-zoom1.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


    //---
      h_finebin_input_qcdmc_lowdphi -> SetMaximum(5) ;
      h_finebin_input_qcdmc_highdphi -> SetMaximum(5) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/finebin-input-qcdmc-zoom2.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      h_kqcd_input_qcdmc_lowdphi -> SetMaximum( 1.10 * ( h_kqcd_input_qcdmc_lowdphi -> GetBinContent( h_kqcd_input_qcdmc_lowdphi -> GetMaximumBin() ) ) ) ;
      h_kqcd_input_qcdmc_highdphi -> SetMaximum( 1.10 * ( h_kqcd_input_qcdmc_highdphi -> GetBinContent( h_kqcd_input_qcdmc_highdphi -> GetMaximumBin() ) ) ) ;
      h_combine_input_qcdmc_lowdphi -> SetMaximum( 1.10 * ( h_combine_input_qcdmc_lowdphi -> GetBinContent( h_combine_input_qcdmc_lowdphi -> GetMaximumBin() ) ) ) ;
      h_combine_input_qcdmc_highdphi -> SetMaximum( 1.10 * ( h_combine_input_qcdmc_highdphi -> GetBinContent( h_combine_input_qcdmc_highdphi -> GetMaximumBin() ) ) ) ;
      h_finebin_input_qcdmc_lowdphi -> SetMaximum( 1.10 * ( h_finebin_input_qcdmc_lowdphi -> GetBinContent( h_finebin_input_qcdmc_lowdphi -> GetMaximumBin() ) ) ) ;
      h_finebin_input_qcdmc_highdphi -> SetMaximum( 1.10 * ( h_finebin_input_qcdmc_highdphi -> GetBinContent( h_finebin_input_qcdmc_highdphi -> GetMaximumBin() ) ) ) ;

      printf("\n\n Saving histograms to outputfiles/qcdmc-input.root\n\n") ;
      TFile* tf_out = new TFile( "outputfiles/qcdmc-input.root", "RECREATE" ) ;
      h_kqcd_input_qcdmc_lowdphi -> Write() ;
      h_kqcd_input_qcdmc_highdphi -> Write() ;
      h_kqcd_input_qcdmc_highlow_ratio -> Write() ;
      h_combine_input_qcdmc_lowdphi -> Write() ;
      h_combine_input_qcdmc_highdphi -> Write() ;
      h_finebin_input_qcdmc_lowdphi -> Write() ;
      h_finebin_input_qcdmc_highdphi -> Write() ;
      tf_out -> Close() ;





   } // make_qcdmc_input_files1

  //================================================================

   int fb_hti_from_mhthti( int mhthti ) {

      if ( mhthti == 1 || mhthti==4 || mhthti==7 ) return 1 ;
      if ( mhthti == 2 || mhthti==5 || mhthti==8 || mhthti==10 ) return 2 ;
      if ( mhthti == 3 || mhthti==6 || mhthti==9 || mhthti==11 ) return 3 ;
      return -99 ;

   } // fb_hti_from_mhthti

  //==========================

   int fb_mhti_from_mhthti( int mhthti ) {

      if ( mhthti==1 || mhthti==2 || mhthti==3 ) return 1 ;
      if ( mhthti==4 || mhthti==5 || mhthti==6 ) return 2 ;
      if ( mhthti==7 || mhthti==8 || mhthti==9 ) return 3 ;
      if ( mhthti==10 || mhthti==11 ) return 4 ;
      return -99 ;

   } // fb_mhti_from_mhthti


  //========================================================================================================

   void fill_fb_arrays( int sb_mhtht_bi, int sb_njets_bi, int sb_nbjets_bi ) {

      if ( sb_mhtht_bi < 1 || sb_mhtht_bi > 6 ) gSystem -> Exit(0) ;
      if ( sb_njets_bi < 1 || sb_njets_bi > 3 ) gSystem -> Exit(0) ;

      /////////nfb = 0 ;
      if ( sb_njets_bi == 1 ) {
         if ( sb_mhtht_bi == 1 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 1 ; fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 1 ; fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 1 ; fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 4 ; fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 4 ; fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 4 ; fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;  nfb ++ ;
         } else if ( sb_mhtht_bi == 2 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 2 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 2 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 2 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 5 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 5 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 5 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         } else if ( sb_mhtht_bi == 3 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 3 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 3 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 3 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 6 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 6 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 6 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         } else if ( sb_mhtht_bi == 4 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 7 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 8 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 7 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 8 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 7 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 8 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         } else if ( sb_mhtht_bi == 5 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 9 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 9 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 9 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 6 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 10 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 11 ;  fbji[nfb] = 1 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 10 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 11 ;  fbji[nfb] = 2 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 10 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 11 ;  fbji[nfb] = 3 ;  fbbi[nfb] = sb_nbjets_bi ;   nfb ++ ;
         }
      } else if ( sb_njets_bi == 2 ) {
         if ( sb_mhtht_bi == 1 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 1 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 4 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 2 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 2 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 5 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 3 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 3 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 6 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 4 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 7 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 8 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         } else if ( sb_mhtht_bi == 5 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 9 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 6 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 10 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 11 ;  fbji[nfb] = 4 ;  fbbi[nfb] = sb_nbjets_bi ;    nfb ++ ;
         }
      } else {
         if ( sb_mhtht_bi == 1 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 1 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 1 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 4 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 2 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 2 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 5 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 3 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 1 ;  fbmhtiti[nfb] = 3 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 2 ;  fbmhtiti[nfb] = 6 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 4 ) {
            fbhi[nfb] = 1 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 7 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 2 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 8 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 5 ) {
            fbhi[nfb] = 3 ;  fbmi[nfb] = 3 ;  fbmhtiti[nfb] = 9 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         } else if ( sb_mhtht_bi == 6 ) {
            fbhi[nfb] = 2 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 10 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
            fbhi[nfb] = 3 ;  fbmi[nfb] = 4 ;  fbmhtiti[nfb] = 11 ;  fbji[nfb] = 5 ;  fbbi[nfb] = sb_nbjets_bi ;     nfb ++ ;
         }
      }

   } // fill_fb_arrays

  //========================================================================================================







