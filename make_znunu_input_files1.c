
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

   float find_line_val( ifstream& ifs, const char* key ) ;

   void fill_fb_arrays( int sb_mhtht_bi, int sb_njets_bi, int sb_nbjets_bi ) ;

  //-------

   void make_znunu_input_files1( const char* infile_lowdphi = "non-QCD-bg-inputs/znunu-LDP-v1.txt",
                                 const char* infile_highdphi = "non-QCD-bg-inputs/znunu-HDP-v1.txt",
                                  const char* outfile_kqcd_fit = "outputfiles/kqcd-input-znunu.txt",
                                  const char* outfile_combine = "outputfiles/combine-input-znunu.txt" ) {

      ifstream ifs_lowdphi ;
      ifstream ifs_highdphi ;

      ifs_lowdphi.open( infile_lowdphi ) ;
      if ( !ifs_lowdphi.good() ) {
         printf("\n\n *** Problem with input file : %s\n\n", infile_lowdphi ) ; return ;
      }
      ifs_highdphi.open( infile_highdphi ) ;
      if ( !ifs_highdphi.good() ) {
         printf("\n\n *** Problem with input file : %s\n\n", infile_highdphi ) ; return ;
      }

      FILE* ofp_kqcd_fit ;
      if ( (ofp_kqcd_fit=fopen( outfile_kqcd_fit, "w" ))==NULL ) {
         printf("\n\n *** Problem opening output file %s\n\n", outfile_kqcd_fit ) ;
      }

      TH1F* h_kqcd_input_znunu_lowdphi = new TH1F( "h_kqcd_input_znunu_lowdphi", "Kqcd fit input, Znunu, low DeltaPhi", 55, 0.5, 55.5 ) ;
      TH1F* h_kqcd_input_znunu_highdphi = new TH1F( "h_kqcd_input_znunu_highdphi", "Kqcd fit input, Znunu, high DeltaPhi", 55, 0.5, 55.5 ) ;

      for ( int nji=1; nji<=5; nji++ ) {
         for ( int mhthti=1; mhthti<=11; mhthti++ ) {

            float nbsum_lowdphi_val(0.) ;
            float nbsum_lowdphi_err2(0.) ;

            float nbsum_highdphi_val(0.) ;
            float nbsum_highdphi_err2(0.) ;

            for ( int nbi=0; nbi<=3; nbi++ ) {

               int owen_bi = 44*(nji-1) + nbi*11 + mhthti ;

               char znunu_file_key[100] ;

               sprintf( znunu_file_key, "Ynj%d_nb%d_kin%d", nji, nbi, mhthti ) ;
               float val_lowdphi  = find_line_val( ifs_lowdphi, znunu_file_key ) ;
               float val_highdphi = find_line_val( ifs_highdphi, znunu_file_key ) ;

               sprintf( znunu_file_key, "Enj%d_nb%d_kin%d", nji, nbi, mhthti ) ;
               float err_lowdphi  = find_line_val( ifs_lowdphi, znunu_file_key ) ;
               float err_highdphi = find_line_val( ifs_highdphi, znunu_file_key ) ;


               nbsum_lowdphi_val += val_lowdphi ;
               nbsum_lowdphi_err2 += pow( err_lowdphi, 2. ) ;

               nbsum_highdphi_val += val_highdphi ;
               nbsum_highdphi_err2 += pow( err_highdphi, 2. ) ;

            } // nbi

            int bi = 11*(nji-1) + mhthti ;

            int mhti = fb_mhti_from_mhthti( mhthti ) ;
            int hti  = fb_hti_from_mhthti( mhthti ) ;

            char binlabel[100] ;
            sprintf( binlabel, "FB-Njet%d-Nbsum-MHT%d-HT%d  %3d", nji, mhti, hti, bi ) ;

            h_kqcd_input_znunu_lowdphi -> SetBinContent( bi, nbsum_lowdphi_val ) ;
            h_kqcd_input_znunu_lowdphi -> SetBinError( bi, sqrt(nbsum_lowdphi_err2) ) ;
            h_kqcd_input_znunu_lowdphi -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;

            h_kqcd_input_znunu_highdphi -> SetBinContent( bi, nbsum_highdphi_val ) ;
            h_kqcd_input_znunu_highdphi -> SetBinError( bi, sqrt(nbsum_highdphi_err2) ) ;
            h_kqcd_input_znunu_highdphi -> GetXaxis() -> SetBinLabel( bi, binlabel ) ;

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

      h_kqcd_input_znunu_lowdphi -> GetXaxis() -> LabelsOption( "v" ) ;
      h_kqcd_input_znunu_highdphi -> GetXaxis() -> LabelsOption( "v" ) ;

      h_kqcd_input_znunu_lowdphi -> SetMinimum(0.1) ;
      h_kqcd_input_znunu_highdphi -> SetMinimum(0.1) ;

      h_kqcd_input_znunu_lowdphi -> SetFillColor( kGreen-7 ) ;
      h_kqcd_input_znunu_highdphi -> SetFillColor( kGreen-7 ) ;

      gStyle -> SetPadBottomMargin( 0.38 ) ;
      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      TCanvas* can_kqcd = new TCanvas( "can_kqcd", "Kqcd fit inputs, Znunu", 1100, 900 ) ;
      can_kqcd -> Divide(1,2) ;

      char pdffile[10000] ;

    //---
      can_kqcd -> cd(1) ;
      h_kqcd_input_znunu_lowdphi -> Draw() ;
      h_kqcd_input_znunu_lowdphi -> Draw( "hist same" ) ;
      h_kqcd_input_znunu_lowdphi -> Draw( "same" ) ;
      h_kqcd_input_znunu_lowdphi -> Draw( "axis same" ) ;
      h_kqcd_input_znunu_lowdphi -> Draw( "axig same" ) ;

      can_kqcd -> cd(2) ;
      h_kqcd_input_znunu_highdphi -> Draw() ;
      h_kqcd_input_znunu_highdphi -> Draw( "hist same" ) ;
      h_kqcd_input_znunu_highdphi -> Draw( "same" ) ;
      h_kqcd_input_znunu_highdphi -> Draw( "axis same" ) ;
      h_kqcd_input_znunu_highdphi -> Draw( "axig same" ) ;

      sprintf( pdffile, "outputfiles/kqcd-input-znunu-liny.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;

    //---
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(1) ;
      can_kqcd -> cd(2) ;
      gPad -> SetLogy(1) ;

      sprintf( pdffile, "outputfiles/kqcd-input-znunu-logy.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;


    //---
      h_kqcd_input_znunu_lowdphi -> SetMaximum(20) ;
      h_kqcd_input_znunu_highdphi -> SetMaximum(20) ;
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/kqcd-input-znunu-zoom1.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;


    //---
      h_kqcd_input_znunu_lowdphi -> SetMaximum(5) ;
      h_kqcd_input_znunu_highdphi -> SetMaximum(5) ;
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/kqcd-input-znunu-zoom2.pdf" ) ;
      can_kqcd -> Update() ; can_kqcd -> Draw() ; can_kqcd -> SaveAs( pdffile ) ;






   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FILE* ofp_combine ;
      if ( (ofp_combine=fopen( outfile_combine, "w" ))==NULL ) {
         printf("\n\n *** Problem opening output file %s\n\n", outfile_combine ) ;
      }

      TH1F* h_combine_input_znunu_lowdphi = new TH1F( "h_combine_input_znunu_lowdphi", "Combine input, Znunu, low DeltaPhi", 72, 0.5, 72.5 ) ;
      TH1F* h_combine_input_znunu_highdphi = new TH1F( "h_combine_input_znunu_highdphi", "Combine input, Znunu, high DeltaPhi", 72, 0.5, 72.5 ) ;

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

                  char znunu_file_key[100] ;

                  sprintf( znunu_file_key, "Ynj%d_nb%d_kin%d", fbji[fbi], fbbi[fbi], fbmhtiti[fbi] ) ;
                  float val_lowdphi  = find_line_val( ifs_lowdphi, znunu_file_key ) ;
                  float val_highdphi = find_line_val( ifs_highdphi, znunu_file_key ) ;

                  sprintf( znunu_file_key, "Enj%d_nb%d_kin%d", fbji[fbi], fbbi[fbi], fbmhtiti[fbi] ) ;
                  float err_lowdphi  = find_line_val( ifs_lowdphi, znunu_file_key ) ;
                  float err_highdphi = find_line_val( ifs_highdphi, znunu_file_key ) ;


                  fbsum_lowdphi_val += val_lowdphi ;
                  fbsum_lowdphi_err2 += pow( err_lowdphi, 2. ) ;

                  fbsum_highdphi_val += val_highdphi ;
                  fbsum_highdphi_err2 += pow( err_highdphi, 2. ) ;


               } // fbi

               char binlabel[100] ;
               sprintf( binlabel, "SB-Njet%d-Nb%d-%-9s  %3d", sb_nji, sb_nbi, sb_mhthti_string, owen_sbi ) ;

               h_combine_input_znunu_lowdphi -> SetBinContent( owen_sbi, fbsum_lowdphi_val ) ;
               h_combine_input_znunu_lowdphi -> SetBinError( owen_sbi, sqrt(fbsum_lowdphi_err2) ) ;

               h_combine_input_znunu_highdphi -> SetBinContent( owen_sbi, fbsum_highdphi_val ) ;
               h_combine_input_znunu_highdphi -> SetBinError( owen_sbi, sqrt(fbsum_highdphi_err2) ) ;

               h_combine_input_znunu_lowdphi  -> GetXaxis() -> SetBinLabel( owen_sbi, binlabel ) ;
               h_combine_input_znunu_highdphi -> GetXaxis() -> SetBinLabel( owen_sbi, binlabel ) ;


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

      h_combine_input_znunu_lowdphi  -> GetXaxis() -> LabelsOption( "v" ) ;
      h_combine_input_znunu_highdphi -> GetXaxis() -> LabelsOption( "v" ) ;

      h_combine_input_znunu_lowdphi -> SetMinimum(0.1) ;
      h_combine_input_znunu_highdphi -> SetMinimum(0.1) ;

      h_combine_input_znunu_lowdphi -> SetFillColor( kGreen-7 ) ;
      h_combine_input_znunu_highdphi -> SetFillColor( kGreen-7 ) ;

      gStyle -> SetPadBottomMargin( 0.35 ) ;
      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      TCanvas* can_combine = new TCanvas( "can_combine", "combine fit inputs, Znunu", 1100, 900 ) ;
      can_combine -> Divide(1,2) ;


    //---
      can_combine -> cd(1) ;
      h_combine_input_znunu_lowdphi -> Draw() ;
      h_combine_input_znunu_lowdphi -> Draw( "hist same" ) ;
      h_combine_input_znunu_lowdphi -> Draw( "same" ) ;
      h_combine_input_znunu_lowdphi -> Draw( "axis same" ) ;
      h_combine_input_znunu_lowdphi -> Draw( "axig same" ) ;

      can_combine -> cd(2) ;
      h_combine_input_znunu_highdphi -> Draw() ;
      h_combine_input_znunu_highdphi -> Draw( "hist same" ) ;
      h_combine_input_znunu_highdphi -> Draw( "same" ) ;
      h_combine_input_znunu_highdphi -> Draw( "axis same" ) ;
      h_combine_input_znunu_highdphi -> Draw( "axig same" ) ;

      sprintf( pdffile, "outputfiles/combine-input-znunu-liny.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;

    //---
      can_combine -> cd(1) ;
      gPad -> SetLogy(1) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(1) ;

      sprintf( pdffile, "outputfiles/combine-input-znunu-logy.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


    //---
      h_combine_input_znunu_lowdphi -> SetMaximum(20) ;
      h_combine_input_znunu_highdphi -> SetMaximum(20) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/combine-input-znunu-zoom1.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;


    //---
      h_combine_input_znunu_lowdphi -> SetMaximum(5) ;
      h_combine_input_znunu_highdphi -> SetMaximum(5) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;
      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/combine-input-znunu-zoom2.pdf" ) ;
      can_combine -> Update() ; can_combine -> Draw() ; can_combine -> SaveAs( pdffile ) ;



      h_kqcd_input_znunu_lowdphi -> SetMaximum( 1.10 * ( h_kqcd_input_znunu_lowdphi -> GetBinContent( h_kqcd_input_znunu_lowdphi -> GetMaximumBin() ) ) ) ;
      h_kqcd_input_znunu_highdphi -> SetMaximum( 1.10 * ( h_kqcd_input_znunu_highdphi -> GetBinContent( h_kqcd_input_znunu_highdphi -> GetMaximumBin() ) ) ) ;
      h_combine_input_znunu_lowdphi -> SetMaximum( 1.10 * ( h_combine_input_znunu_lowdphi -> GetBinContent( h_combine_input_znunu_lowdphi -> GetMaximumBin() ) ) ) ;
      h_combine_input_znunu_highdphi -> SetMaximum( 1.10 * ( h_combine_input_znunu_highdphi -> GetBinContent( h_combine_input_znunu_highdphi -> GetMaximumBin() ) ) ) ;

      printf("\n\n Saving histograms to outputfiles/znunu-input.root\n\n") ;
      TFile* tf_out = new TFile( "outputfiles/znunu-input.root", "RECREATE" ) ;
      h_kqcd_input_znunu_lowdphi -> Write() ;
      h_kqcd_input_znunu_highdphi -> Write() ;
      h_combine_input_znunu_lowdphi -> Write() ;
      h_combine_input_znunu_highdphi -> Write() ;
      tf_out -> Close() ;


   } // make_znunu_input_files1

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






