

#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"

#include "histio.c"


  //--------------------------------------------------------------------------------------------------------

   int nhtbins = 4 ;   // first bin is not used in analysis.
   double htbins[5] = { 0., 500., 800., 1200., 20000. } ;  // first bin is not used in analysis.

   /////////int nmhtbins = 6 ;  // first two bins are not used in analysis.
   /////////double mhtbins[7] = { 0., 150., 200., 300., 500., 750., 20000. } ;    // first two bins are not used in analysis.
   int nmhtbins ;
   double mhtbins[10] ;

   int nnbbins = 4 ;
   double nbbins[5] = { -0.5, 0.5, 1.5, 2.5, 20. } ;

   int nnjetbins = 6 ; // first bin is not used in analysis.
   double njetbins[7] = { 0., 3.5, 4.5, 5.5, 6.5, 8.5, 20. } ;  // first bin is not used in analysis.

  //--------------------------------------------------------------------------------------------------------

   int n_merged_bins ;
   int merged_bin_index[8][8][8][8] ; // nbi, nji, mbi, hbi
   char merged_bin_label[300][20] ;

  //--------------------------------------------------------------------------------------------------------

   int n_samples ;
   char sname[100][100] ;
   int scolor[100] ;
   bool isll[100] ;

   TCanvas* can ;

   void draw_stack( const char* configstr, int nbi, int nji ) ;
   TH1F* make_1d( TH3F* h3p, int nbi, int nji ) ;
   void make_legend(  ) ;
   void setup_merged_bins() ;

  //--------------------------------------------------------------------------------------------------------

   void draw_bg_hists4( const char* infile = "outputfiles/fill-bg-hists4.root" ) {

         nmhtbins = 5 ;
         mhtbins[0] = 0. ;
         mhtbins[1] = 200. ;
         mhtbins[2] = 300. ;
         mhtbins[3] = 500. ;
         mhtbins[4] = 750. ;
         mhtbins[5] = 20000. ;

      gStyle -> SetPadBottomMargin( 0.20 ) ;
      gStyle -> SetOptStat(0) ;

      setup_merged_bins() ;

      gDirectory -> Delete( "h*" ) ;
      printf("\n\n Reading in 2D MHT vs HT histograms from %s\n\n", infile ) ;
      loadHist( infile ) ;

      int n_selections = 4 ;
      char selname[4][100] ;
      sprintf( selname[0], "zl" ) ;
      sprintf( selname[1], "sl" ) ;
      sprintf( selname[2], "ldp" ) ;
      sprintf( selname[3], "slldp" ) ;

      int si(0) ;
      sprintf( sname[si], "qcd" )     ; isll[si] = false ; scolor[si] = kRed + 3 ; si++ ;
      sprintf( sname[si], "ttbar" )   ; isll[si] = true  ; scolor[si] = kBlue - 4 ; si++ ;
      sprintf( sname[si], "wjets" )   ; isll[si] = true  ; scolor[si] = kGreen + 2 ; si++ ;
      sprintf( sname[si], "sngltop" ) ; isll[si] = true  ; scolor[si] = kYellow - 7 ; si++ ;
      sprintf( sname[si], "znunu" )   ; isll[si] = false ; scolor[si] = kOrange - 3 ; si++ ;
      n_samples = si ;

      can = (TCanvas*) gDirectory -> FindObject( "can_draw_bg_hists" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_draw_bg_hists", "BG hists", 700, 600 ) ;
      can -> Clear() ;

      for ( int selind=0; selind<n_selections; selind++ ) {

         for ( int nbi=0; nbi<nnbbins; nbi++ ) {

            for ( int nji=1; nji<nnjetbins; nji++ ) {

               char configstr[1000] ;
               sprintf( configstr, "nb%d_%s", nbi, selname[selind] ) ;
               printf("  Drawing stack for %s\n", configstr ) ;

               draw_stack( configstr, nbi, nji ) ;

            } // nji

         } // nbi.

      } // selind.

      make_legend( ) ;

      TString savefilets( infile ) ;
      char newend[100] ;
      sprintf( newend, "-postdraw.root" ) ;
      savefilets.ReplaceAll( ".root", newend ) ;
      printf("\n\n Saving all histograms as %s\n", savefilets.Data() ) ;
      saveHist( savefilets.Data(), "h*" ) ;


   } // draw_bg_hists4

  //====================================================================================================================

   void draw_stack( const char* configstr, int nbi, int nji ) {

      char hname[1000] ;

      sprintf( hname, "h_stack_%s_nj%d", configstr, nji ) ;
      THStack* h_stack = new THStack( hname, configstr ) ;

      TH1F* h_bgsum(0x0) ;
      TH1F* h_llbgsum(0x0) ;

      TH1F* h1d_samples[100] ;

      for ( int si=0; si<n_samples; si++ ) {

         sprintf( hname, "h_njvsmhtvsht_%s_%s", configstr, sname[si] ) ;
         printf("   hname: %s\n", hname ) ;

         TH3F* h3d = (TH3F*) gDirectory -> FindObject( hname ) ;
         if ( h3d == 0x0 ) { printf("\n\n *** Missing histogram %s\n\n", hname ) ; return ; }

         TH1F* h1d = make_1d( h3d, nbi, nji ) ;
         h1d_samples[si] = h1d ;
         h1d -> SetFillColor( scolor[si] ) ;
         h_stack -> Add( h1d ) ;

         if ( si == 0 ) {
            sprintf( hname, "h_bgsum_%s_nj%d", configstr, nji ) ;
            h_bgsum = (TH1F*) h1d -> Clone( hname ) ;
            h_bgsum -> SetTitle( configstr ) ;
         } else {
            h_bgsum -> Add( h1d ) ;
         }

         if ( isll[si] ) {
            if ( h_llbgsum == 0x0 ) {
               sprintf( hname, "h_llbgsum_%s_nj%d", configstr, nji ) ;
               h_llbgsum = (TH1F*) h1d -> Clone( hname ) ;
               h_llbgsum -> SetTitle( configstr ) ;
            } else {
               h_llbgsum -> Add( h1d ) ;
            }
         }

      } // si

      sprintf( hname, "h_frac_stack_%s_nj%d", configstr, nji ) ;
      THStack* h_frac_stack = new THStack( hname, configstr ) ;

      TH1F* h_frac_sum(0x0) ;

      for ( int si=0; si<n_samples; si++ ) {
         sprintf( hname, "%s_fr", h1d_samples[si] -> GetName() ) ;
         TH1F* h1d_fr = (TH1F*) h1d_samples[si] -> Clone( hname ) ;
         char htitle[1000] ;
         sprintf( htitle, "%s fraction", h1d_samples[si] -> GetTitle() ) ;
         h1d_fr -> SetTitle( htitle ) ;
         for ( int bi=1; bi<=h1d_samples[0] -> GetNbinsX(); bi++ ) {
            float sum_val = h_bgsum -> GetBinContent( bi ) ;
            if ( sum_val <= 0 ) continue ;
            float sum_err = h_bgsum -> GetBinError( bi ) ;
            float comp_val = h1d_samples[si] -> GetBinContent( bi ) ;
            float comp_err = h1d_samples[si] -> GetBinError( bi ) ;
            float fr_val = comp_val / sum_val ;
            float rem_val = sum_val - comp_val ;
            float rem_err = sqrt( pow(sum_err,2) - pow(comp_err,2) ) ;
            float fr_err = 0. ;
            if ( rem_val > 0. && comp_val > 0. ) fr_err = sqrt( pow( (rem_val*comp_err) / (sum_val*sum_val), 2 ) + pow( (comp_val*rem_err) / (sum_val*sum_val), 2 ) ) ;
            h1d_fr -> SetBinContent( bi, fr_val ) ;
            h1d_fr -> SetBinError( bi, fr_err ) ;
         } // bi
         h_frac_stack -> Add( h1d_fr ) ;
         if ( si == 0 ) {
            sprintf( hname, "h_bgfracsum_%s_nj%d", configstr, nji ) ;
            h_frac_sum = (TH1F*) h1d_fr -> Clone( hname ) ;
            h_frac_sum -> SetTitle( configstr ) ;
         } else {
            h_frac_sum -> Add( h1d_fr ) ;
         }
      } // si

      char fname[10000] ;

      h_bgsum -> SetMinimum( 0.1 ) ;
      h_frac_sum -> SetMaximum( 1.3 ) ;

      h_bgsum -> Draw( ) ;
      h_stack -> Draw( "hist same" ) ;
      h_bgsum -> Draw( "same" ) ;

      gPad -> SetLogy(0) ;
      sprintf( fname, "outputfiles/bgplot-nj%d-%s-liny.pdf", nji, configstr ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

      gPad -> SetLogy(1) ;
      sprintf( fname, "outputfiles/bgplot-nj%d-%s-logy.pdf", nji, configstr ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

      h_bgsum -> SetMaximum(25.) ;
      h_bgsum -> Draw( ) ;
      h_stack -> Draw( "hist same" ) ;
      h_bgsum -> Draw( "same" ) ;

      gPad -> SetLogy(0) ;
      sprintf( fname, "outputfiles/bgplot-nj%d-%s-zoom.pdf", nji, configstr ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;


      h_bgsum -> SetMaximum(5.) ;
      h_bgsum -> Draw( ) ;
      h_stack -> Draw( "hist same" ) ;
      h_bgsum -> Draw( "same" ) ;

      gPad -> SetLogy(0) ;
      sprintf( fname, "outputfiles/bgplot-nj%d-%s-zoom2.pdf", nji, configstr ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

      h_frac_sum -> Draw( "hist" ) ;
      h_frac_stack -> Draw( "hist same" ) ;
      h_frac_stack -> Draw( "same" ) ;
      sprintf( fname, "outputfiles/bgplot-nj%d-%s-frac.pdf", nji, configstr ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;


   } // draw_stack


  //====================================================================================================================

   TH1F* make_1d( TH3F* h3p, int nbi, int nji ) {

      if ( h3p == 0x0 ) return 0x0 ;

      char hname1d[1000] ;
      sprintf( hname1d, "%s_nj%d_1d", h3p->GetName(), nji ) ;
      TH1F* h1p = new TH1F( hname1d, h3p->GetTitle(), n_merged_bins, 0.5, n_merged_bins+0.5 ) ;

      for ( int bi=1; bi<=n_merged_bins; bi++ ) {

         float val(0.) ;
         float err2(0.) ;

         for ( int mbi=1; mbi<nmhtbins; mbi++ ) {

            for ( int hbi=1; hbi<nhtbins; hbi++ ) {

               if ( merged_bin_index[nbi][nji][mbi][hbi] == bi ) {

                  int h3xbi = hbi+1 ;
                  int h3ybi = mbi+1 ;
                  int h3zbi = nji+1 ;
                  val += h3p -> GetBinContent( h3xbi, h3ybi, h3zbi ) ;
                  err2 += pow( h3p -> GetBinError( h3xbi, h3ybi, h3zbi ), 2. ) ;

               }

            } // hbi
         } // mbi

         h1p -> GetXaxis() -> SetBinLabel( bi, merged_bin_label[bi] ) ;

         h1p -> SetBinContent( bi, val ) ;
         h1p -> SetBinError( bi, sqrt(err2) ) ;
      }

      h1p -> GetXaxis() -> LabelsOption( "v" ) ; ;

      return h1p ;

   } // make_1d

  //====================================================================================================================

   void make_legend( ) {

      char hname[1000] ;

      TH1F* h_legend_dummy = new TH1F( "h_legend_dummy", "", 2, 0., 1. ) ;

      TLegend* legend = new TLegend( 0.1, 0.2, 0.9, 0.9 ) ;

      for ( int si=n_samples-1; si>=0; si-- ) {

         sprintf( hname, "h_njvsmhtvsht_nb0_zl_%s_nj1_1d", sname[si] ) ;
         printf("   hname: %s\n", hname ) ;

         TH1F* h1d = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( h1d == 0x0 ) { printf("\n\n *** Missing histogram %s\n\n", hname ) ; return ; }

         legend -> AddEntry( h1d, sname[si] ) ;

      } // si

      h_legend_dummy -> Draw( "AH" ) ;
      legend -> Draw() ;

      char fname[10000] ;
      sprintf( fname, "outputfiles/bgplot-legend.pdf" ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

   } // draw_stack

  //====================================================================================================================

   void setup_merged_bins() {

      printf("\n\n Setting up merged bins.\n\n") ;


      for ( int nbi=0; nbi<nnbbins; nbi++ ) {
         for ( int nji=1; nji<nnjetbins; nji++ ) {

            int mbi, hbi ;

            int bi = 1 ;


            sprintf( merged_bin_label[bi], "" ) ; bi ++ ; // blank bin in histogram

            sprintf( merged_bin_label[bi], "MHT1_HT1" ) ;
            mbi = 1 ; hbi = 1 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT1_HT2" ) ;
            mbi = 1 ; hbi = 2 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT1_HT3" ) ;
            mbi = 1 ; hbi = 3 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;


            sprintf( merged_bin_label[bi], "" ) ; bi ++ ; // blank bin in histogram

            sprintf( merged_bin_label[bi], "MHT2_HT1" ) ;
            mbi = 2 ; hbi = 1 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT2_HT2" ) ;
            mbi = 2 ; hbi = 2 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT2_HT3" ) ;
            mbi = 2 ; hbi = 3 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;


            sprintf( merged_bin_label[bi], "" ) ; bi ++ ; // blank bin in histogram

            sprintf( merged_bin_label[bi], "MHT3_HT1" ) ;
            mbi = 3 ; hbi = 1 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT3_HT2" ) ;
            mbi = 3 ; hbi = 2 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT3_HT3" ) ;
            mbi = 3 ; hbi = 3 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;


            sprintf( merged_bin_label[bi], "" ) ; bi ++ ; // blank bin in histogram

            sprintf( merged_bin_label[bi], "MHT4_HT1" ) ;
            mbi = 4 ; hbi = 1 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT4_HT2" ) ;
            mbi = 4 ; hbi = 2 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "MHT4_HT3" ) ;
            mbi = 4 ; hbi = 3 ;
            merged_bin_index[nbi][nji][mbi][hbi] = bi ;
            printf("  %3d : %s\n", bi, merged_bin_label[bi] ) ;
            bi ++ ;

            sprintf( merged_bin_label[bi], "" ) ; bi ++ ; // blank bin in histogram


            n_merged_bins = bi -1 ;



           //-- unused bins
            mbi = 4 ; hbi = 1 ;
            merged_bin_index[nbi][nji][mbi][hbi] = -1 ;


         } // nji
      } // nbi


      printf("\n\n") ;

   } // setup_merged_bins.

  //====================================================================================================================







