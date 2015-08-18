
#include "TFile.h"
#include "THStack.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"


   void draw_boundaries_kqcd( float hmin, float hmax ) ;
   void draw_boundaries_combine( float hmin, float hmax ) ;

   void draw_inputs1(
          const char* data_file = "outputfiles/fakedata-input.root",
          const char* sigmc_file = "outputfiles/sigmc-input.root",
          const char* lostlep_file = "outputfiles/lostlep-input.root",
          const char* hadtau_file = "outputfiles/hadtau-input.root",
          const char* znunu_file = "outputfiles/znunu-input.root",
          const char* qcdmc_file = "outputfiles/qcdmc-input.root",
          bool draw_legend = true
   ) {

      char pdf_file[1000] ;

      if ( draw_legend) {
         gStyle -> SetPadRightMargin(0.15) ;
      }

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin(0.37) ;

      gDirectory -> Delete( "h*" ) ;

      TFile* tf_lostlep = new TFile( lostlep_file, "READ" ) ;
      if ( tf_lostlep == 0x0 ) { printf("\n\n *** bad file lostlep_file\n\n" ) ; return ; }

      TFile* tf_hadtau = new TFile( hadtau_file, "READ" ) ;
      if ( tf_hadtau == 0x0 ) { printf("\n\n *** bad file hadtau_file\n\n" ) ; return ; }

      TFile* tf_znunu = new TFile( znunu_file, "READ" ) ;
      if ( tf_znunu == 0x0 ) { printf("\n\n *** bad file znunu_file\n\n" ) ; return ; }

      TFile* tf_qcdmc = new TFile( qcdmc_file, "READ" ) ;
      if ( tf_qcdmc == 0x0 ) { printf("\n\n *** bad file qcdmc_file\n\n" ) ; return ; }

      TFile* tf_data = new TFile( data_file, "READ" ) ;
      if ( tf_data == 0x0 ) { printf("\n\n *** Missing or bad data file.  Will not plot it.\n\n" ) ; }

      TFile* tf_sigmc = new TFile( sigmc_file, "READ" ) ;
      if ( tf_sigmc == 0x0 ) { printf("\n\n *** Missing or bad sigmc file.  Will not plot it.\n\n" ) ; }

      fflush(stdout) ;



   //------

      TH1* h_kqcd_input_lostlep_lowdphi = (TH1*) tf_lostlep -> Get( "h_kqcd_input_lostlep_lowdphi" ) ;
      TH1* h_kqcd_input_hadtau_lowdphi  = (TH1*) tf_hadtau  -> Get( "h_kqcd_input_hadtau_lowdphi" ) ;
      TH1* h_kqcd_input_znunu_lowdphi = (TH1*) tf_znunu -> Get( "h_kqcd_input_znunu_lowdphi" ) ;
      TH1* h_kqcd_input_qcdmc_lowdphi   = (TH1*) tf_qcdmc   -> Get( "h_kqcd_input_qcdmc_lowdphi" ) ;
      TH1* h_kqcd_input_data_lowdphi(0x0) ;
      if ( tf_data != 0x0 ) h_kqcd_input_data_lowdphi = (TH1*) tf_data -> Get( "h_kqcd_input_data_lowdphi" ) ;
      TH1* h_kqcd_input_sigmc_lowdphi(0x0) ;
      if ( tf_sigmc != 0x0 ) {
         h_kqcd_input_sigmc_lowdphi = (TH1*) tf_sigmc -> Get( "h_kqcd_input_sigmc_lowdphi" ) ;
         h_kqcd_input_sigmc_lowdphi -> SetFillColor( kMagenta-3 ) ;
         h_kqcd_input_sigmc_lowdphi -> SetLineColor( kMagenta-3 ) ;
         h_kqcd_input_sigmc_lowdphi -> SetLineWidth( 3 ) ;
         h_kqcd_input_sigmc_lowdphi -> SetFillStyle( 3005 ) ;
      }

      THStack* hs_kqcd_input_lowdphi = new THStack( "hs_kqcd_input_lowdphi", "KQCD fit input, low delta phi" ) ;
      hs_kqcd_input_lowdphi -> Add( h_kqcd_input_qcdmc_lowdphi ) ;
      hs_kqcd_input_lowdphi -> Add( h_kqcd_input_znunu_lowdphi ) ;
      hs_kqcd_input_lowdphi -> Add( h_kqcd_input_hadtau_lowdphi ) ;
      hs_kqcd_input_lowdphi -> Add( h_kqcd_input_lostlep_lowdphi ) ;

      TH1* h_sum_kqcd_input_lowdphi = (TH1*) h_kqcd_input_qcdmc_lowdphi -> Clone( "h_sum_kqcd_input_lowdphi" ) ;
      h_sum_kqcd_input_lowdphi -> Add( h_kqcd_input_znunu_lowdphi ) ;
      h_sum_kqcd_input_lowdphi -> Add( h_kqcd_input_hadtau_lowdphi ) ;
      h_sum_kqcd_input_lowdphi -> Add( h_kqcd_input_lostlep_lowdphi ) ;

      h_sum_kqcd_input_lowdphi -> SetMinimum( 0.1 ) ;


   //------

      TH1* h_kqcd_input_lostlep_highdphi = (TH1*) tf_lostlep -> Get( "h_kqcd_input_lostlep_highdphi" ) ;
      TH1* h_kqcd_input_hadtau_highdphi  = (TH1*) tf_hadtau  -> Get( "h_kqcd_input_hadtau_highdphi" ) ;
      TH1* h_kqcd_input_znunu_highdphi = (TH1*) tf_znunu -> Get( "h_kqcd_input_znunu_highdphi" ) ;
      TH1* h_kqcd_input_qcdmc_highdphi   = (TH1*) tf_qcdmc   -> Get( "h_kqcd_input_qcdmc_highdphi" ) ;
      TH1* h_kqcd_input_data_highdphi(0x0) ;
      if ( tf_data != 0x0 ) h_kqcd_input_data_highdphi = (TH1*) tf_data -> Get( "h_kqcd_input_data_highdphi" ) ;
      TH1* h_kqcd_input_sigmc_highdphi(0x0) ;
      if ( tf_sigmc != 0x0 ) {
         h_kqcd_input_sigmc_highdphi = (TH1*) tf_sigmc -> Get( "h_kqcd_input_sigmc_highdphi" ) ;
         h_kqcd_input_sigmc_highdphi -> SetFillColor( kMagenta-3 ) ;
         h_kqcd_input_sigmc_highdphi -> SetLineColor( kMagenta-3 ) ;
         h_kqcd_input_sigmc_highdphi -> SetLineWidth( 3 ) ;
         h_kqcd_input_sigmc_highdphi -> SetFillStyle( 3005 ) ;
      }

      THStack* hs_kqcd_input_highdphi = new THStack( "hs_kqcd_input_highdphi", "KQCD fit input, high delta phi" ) ;
      hs_kqcd_input_highdphi -> Add( h_kqcd_input_qcdmc_highdphi ) ;
      hs_kqcd_input_highdphi -> Add( h_kqcd_input_znunu_highdphi ) ;
      hs_kqcd_input_highdphi -> Add( h_kqcd_input_hadtau_highdphi ) ;
      hs_kqcd_input_highdphi -> Add( h_kqcd_input_lostlep_highdphi ) ;

      TH1* h_sum_kqcd_input_highdphi = (TH1*) h_kqcd_input_qcdmc_highdphi -> Clone( "h_sum_kqcd_input_highdphi" ) ;
      h_sum_kqcd_input_highdphi -> Add( h_kqcd_input_znunu_highdphi ) ;
      h_sum_kqcd_input_highdphi -> Add( h_kqcd_input_hadtau_highdphi ) ;
      h_sum_kqcd_input_highdphi -> Add( h_kqcd_input_lostlep_highdphi ) ;

      h_sum_kqcd_input_highdphi -> SetMinimum( 0.1 ) ;


      TLegend* legend = new TLegend( 0.87, 0.35, 0.98, 0.90 ) ;
      legend -> AddEntry( h_kqcd_input_lostlep_highdphi, "Lost Lep" ) ;
      legend -> AddEntry( h_kqcd_input_hadtau_highdphi, "Had Tau" ) ;
      legend -> AddEntry( h_kqcd_input_znunu_highdphi, "Znunu" ) ;
      legend -> AddEntry( h_kqcd_input_qcdmc_highdphi, "QCD" ) ;
      if ( h_kqcd_input_data_lowdphi != 0x0 ) legend -> AddEntry( h_kqcd_input_data_highdphi, "data (fake)" ) ;
      if ( h_kqcd_input_sigmc_lowdphi != 0x0 ) legend -> AddEntry( h_kqcd_input_sigmc_highdphi, "Signal MC" ) ;



      TCanvas* can_kqcd = new TCanvas( "can_kqcd", "Kqcd fit inputs", 1300, 1300 ) ;
      can_kqcd -> Divide(1,2) ;


      can_kqcd -> cd(1) ;
      h_sum_kqcd_input_lowdphi -> Draw() ;
      hs_kqcd_input_lowdphi -> Draw( "hist same" ) ;
      hs_kqcd_input_lowdphi -> Draw( "axis same" ) ;
      hs_kqcd_input_lowdphi -> Draw( "axig same" ) ;
      hs_kqcd_input_lowdphi -> Draw( "same" ) ;
      if ( h_kqcd_input_data_lowdphi != 0x0 ) h_kqcd_input_data_lowdphi -> Draw( "same e" ) ;
      if ( h_kqcd_input_sigmc_lowdphi != 0x0 ) h_kqcd_input_sigmc_lowdphi -> Draw( "same hist" ) ;
      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;
      draw_boundaries_kqcd(0.1, 1.1*(hs_kqcd_input_lowdphi->GetMaximum()) ) ;
      if ( draw_legend ) legend -> Draw() ;

      can_kqcd -> cd(2) ;
      h_sum_kqcd_input_highdphi -> Draw() ;
      hs_kqcd_input_highdphi -> Draw( "hist same" ) ;
      hs_kqcd_input_highdphi -> Draw( "axis same" ) ;
      hs_kqcd_input_highdphi -> Draw( "axig same" ) ;
      hs_kqcd_input_highdphi -> Draw( "same" ) ;
      if ( h_kqcd_input_data_highdphi != 0x0 ) h_kqcd_input_data_highdphi -> Draw( "same e" ) ;
      if ( h_kqcd_input_sigmc_highdphi != 0x0 ) h_kqcd_input_sigmc_highdphi -> Draw( "same hist" ) ;
      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;
      draw_boundaries_kqcd(0.1, 1.1*(hs_kqcd_input_highdphi->GetMaximum()) ) ;
      if ( draw_legend ) legend -> Draw() ;

      sprintf( pdf_file, "outputfiles/inputs-kqcd-logy.pdf" ) ;
      can_kqcd -> SaveAs( pdf_file ) ;

    //---

      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;

      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdf_file, "outputfiles/inputs-kqcd-liny.pdf" ) ;
      can_kqcd -> SaveAs( pdf_file ) ;

    //---

      h_sum_kqcd_input_lowdphi -> SetMaximum(20) ;
      h_sum_kqcd_input_highdphi -> SetMaximum(20) ;
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;

      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdf_file, "outputfiles/inputs-kqcd-zoom1.pdf" ) ;
      can_kqcd -> SaveAs( pdf_file ) ;


    //---

      h_sum_kqcd_input_lowdphi -> SetMaximum(5) ;
      h_sum_kqcd_input_highdphi -> SetMaximum(5) ;
      can_kqcd -> cd(1) ;
      gPad -> SetLogy(0) ;

      can_kqcd -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdf_file, "outputfiles/inputs-kqcd-zoom2.pdf" ) ;
      can_kqcd -> SaveAs( pdf_file ) ;


   //--------------------------------------------------------------------------------------------


      TH1* h_combine_input_lostlep_lowdphi = (TH1*) tf_lostlep -> Get( "h_combine_input_lostlep_lowdphi" ) ;
      TH1* h_combine_input_hadtau_lowdphi  = (TH1*) tf_hadtau  -> Get( "h_combine_input_hadtau_lowdphi" ) ;
      TH1* h_combine_input_znunu_lowdphi = (TH1*) tf_znunu -> Get( "h_combine_input_znunu_lowdphi" ) ;
      TH1* h_combine_input_qcdmc_lowdphi   = (TH1*) tf_qcdmc   -> Get( "h_combine_input_qcdmc_lowdphi" ) ;
      TH1* h_combine_input_data_lowdphi(0x0) ;
      if ( tf_data != 0x0 ) h_combine_input_data_lowdphi = (TH1*) tf_data -> Get( "h_combine_input_data_lowdphi" ) ;
      TH1* h_combine_input_sigmc_lowdphi(0x0) ;
      if ( tf_sigmc != 0x0 ) {
         h_combine_input_sigmc_lowdphi = (TH1*) tf_sigmc -> Get( "h_combine_input_sigmc_lowdphi" ) ;
         h_combine_input_sigmc_lowdphi -> SetFillColor( kMagenta-3 ) ;
         h_combine_input_sigmc_lowdphi -> SetLineColor( kMagenta-3 ) ;
         h_combine_input_sigmc_lowdphi -> SetLineWidth( 3 ) ;
         h_combine_input_sigmc_lowdphi -> SetFillStyle(3005) ;
      }

      THStack* hs_combine_input_lowdphi = new THStack( "hs_combine_input_lowdphi", "combine fit input, low delta phi" ) ;
      hs_combine_input_lowdphi -> Add( h_combine_input_qcdmc_lowdphi ) ;
      hs_combine_input_lowdphi -> Add( h_combine_input_znunu_lowdphi ) ;
      hs_combine_input_lowdphi -> Add( h_combine_input_hadtau_lowdphi ) ;
      hs_combine_input_lowdphi -> Add( h_combine_input_lostlep_lowdphi ) ;

      TH1* h_sum_combine_input_lowdphi = (TH1*) h_combine_input_qcdmc_lowdphi -> Clone( "h_sum_combine_input_lowdphi" ) ;
      h_sum_combine_input_lowdphi -> Add( h_combine_input_znunu_lowdphi ) ;
      h_sum_combine_input_lowdphi -> Add( h_combine_input_hadtau_lowdphi ) ;
      h_sum_combine_input_lowdphi -> Add( h_combine_input_lostlep_lowdphi ) ;

      h_sum_combine_input_lowdphi -> SetMinimum( 0.1 ) ;


   //------

      TH1* h_combine_input_lostlep_highdphi = (TH1*) tf_lostlep -> Get( "h_combine_input_lostlep_highdphi" ) ;
      TH1* h_combine_input_hadtau_highdphi  = (TH1*) tf_hadtau  -> Get( "h_combine_input_hadtau_highdphi" ) ;
      TH1* h_combine_input_znunu_highdphi = (TH1*) tf_znunu -> Get( "h_combine_input_znunu_highdphi" ) ;
      TH1* h_combine_input_qcdmc_highdphi   = (TH1*) tf_qcdmc   -> Get( "h_combine_input_qcdmc_highdphi" ) ;
      TH1* h_combine_input_data_highdphi(0x0) ;
      if ( tf_data != 0x0 ) h_combine_input_data_highdphi = (TH1*) tf_data -> Get( "h_combine_input_data_highdphi" ) ;
      TH1* h_combine_input_sigmc_highdphi(0x0) ;
      if ( tf_sigmc != 0x0 ) {
         h_combine_input_sigmc_highdphi = (TH1*) tf_sigmc -> Get( "h_combine_input_sigmc_highdphi" ) ;
         h_combine_input_sigmc_highdphi -> SetFillColor( kMagenta-3 ) ;
         h_combine_input_sigmc_highdphi -> SetLineColor( kMagenta-3 ) ;
         h_combine_input_sigmc_highdphi -> SetLineWidth( 3 ) ;
         h_combine_input_sigmc_highdphi -> SetFillStyle( 3005 ) ;
      }

      THStack* hs_combine_input_highdphi = new THStack( "hs_combine_input_highdphi", "combine fit input, high delta phi" ) ;
      hs_combine_input_highdphi -> Add( h_combine_input_qcdmc_highdphi ) ;
      hs_combine_input_highdphi -> Add( h_combine_input_znunu_highdphi ) ;
      hs_combine_input_highdphi -> Add( h_combine_input_hadtau_highdphi ) ;
      hs_combine_input_highdphi -> Add( h_combine_input_lostlep_highdphi ) ;

      TH1* h_sum_combine_input_highdphi = (TH1*) h_combine_input_qcdmc_highdphi -> Clone( "h_sum_combine_input_highdphi" ) ;
      h_sum_combine_input_highdphi -> Add( h_combine_input_znunu_highdphi ) ;
      h_sum_combine_input_highdphi -> Add( h_combine_input_hadtau_highdphi ) ;
      h_sum_combine_input_highdphi -> Add( h_combine_input_lostlep_highdphi ) ;

      h_sum_combine_input_highdphi -> SetMinimum( 0.1 ) ;


      TLegend* legend_combine = new TLegend( 0.87, 0.35, 0.98, 0.90 ) ;
      legend_combine -> AddEntry( h_combine_input_lostlep_highdphi, "Lost Lep" ) ;
      legend_combine -> AddEntry( h_combine_input_hadtau_highdphi, "Had Tau" ) ;
      legend_combine -> AddEntry( h_combine_input_znunu_highdphi, "Znunu" ) ;
      legend_combine -> AddEntry( h_combine_input_qcdmc_highdphi, "QCD" ) ;
      if ( h_combine_input_data_lowdphi != 0x0 ) legend_combine -> AddEntry( h_combine_input_data_highdphi, "data (fake)" ) ;
      if ( h_combine_input_sigmc_lowdphi != 0x0 ) legend_combine -> AddEntry( h_combine_input_sigmc_highdphi, "Signal MC" ) ;



      TCanvas* can_combine = new TCanvas( "can_combine", "combine fit inputs", 1300, 1300 ) ;
      can_combine -> Divide(1,2) ;


      can_combine -> cd(1) ;
      h_sum_combine_input_lowdphi -> Draw() ;
      hs_combine_input_lowdphi -> Draw( "hist same" ) ;
      hs_combine_input_lowdphi -> Draw( "axis same" ) ;
      hs_combine_input_lowdphi -> Draw( "axig same" ) ;
      hs_combine_input_lowdphi -> Draw( "same" ) ;
      if ( h_combine_input_data_lowdphi != 0x0 ) h_combine_input_data_lowdphi -> Draw( "same e" ) ;
      if ( h_combine_input_sigmc_lowdphi != 0x0 ) h_combine_input_sigmc_lowdphi -> Draw( "same hist" ) ;
      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;
      draw_boundaries_combine(0.1, 1.1*(hs_combine_input_lowdphi->GetMaximum()) ) ;
      if ( draw_legend ) legend_combine -> Draw() ;

      can_combine -> cd(2) ;
      h_sum_combine_input_highdphi -> Draw() ;
      hs_combine_input_highdphi -> Draw( "hist same" ) ;
      hs_combine_input_highdphi -> Draw( "axis same" ) ;
      hs_combine_input_highdphi -> Draw( "axig same" ) ;
      hs_combine_input_highdphi -> Draw( "same" ) ;
      if ( h_combine_input_data_highdphi != 0x0 ) h_combine_input_data_highdphi -> Draw( "same e" ) ;
      if ( h_combine_input_sigmc_highdphi != 0x0 ) h_combine_input_sigmc_highdphi -> Draw( "same hist" ) ;
      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;
      draw_boundaries_combine(0.1, 1.1*(hs_combine_input_highdphi->GetMaximum()) ) ;
      if ( draw_legend ) legend_combine -> Draw() ;

      sprintf( pdf_file, "outputfiles/inputs-combine-logy.pdf" ) ;
      can_combine -> SaveAs( pdf_file ) ;

    //---

      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;

      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdf_file, "outputfiles/inputs-combine-liny.pdf" ) ;
      can_combine -> SaveAs( pdf_file ) ;

    //---

      h_sum_combine_input_lowdphi -> SetMaximum(20) ;
      h_sum_combine_input_highdphi -> SetMaximum(20) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;

      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdf_file, "outputfiles/inputs-combine-zoom1.pdf" ) ;
      can_combine -> SaveAs( pdf_file ) ;


    //---

      h_sum_combine_input_lowdphi -> SetMaximum(5) ;
      h_sum_combine_input_highdphi -> SetMaximum(5) ;
      can_combine -> cd(1) ;
      gPad -> SetLogy(0) ;

      can_combine -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdf_file, "outputfiles/inputs-combine-zoom2.pdf" ) ;
      can_combine -> SaveAs( pdf_file ) ;


   } // draw_inputs1

//==========================

 void draw_boundaries_kqcd( float hmin, float hmax ) {

    TLine* l1 = new TLine() ;

    for ( int i=1; i<=4; i++ ) {
       float x = 11*i+0.5 ;
       l1 -> SetLineWidth(2) ;
       l1 -> SetLineStyle(1) ;
       l1 -> DrawLine( x, hmin, x, hmax ) ;
    }


 } // draw_boundaries_kqcd

//==========================

 void draw_boundaries_combine( float hmin, float hmax ) {

    TLine* l1 = new TLine() ;

    for ( int i=1; i<=2; i++ ) {
       float x = 24*i+0.5 ;
       l1 -> SetLineWidth(3) ;
       l1 -> SetLineStyle(1) ;
       l1 -> DrawLine( x, hmin, x, hmax ) ;
    }

    for ( int i=1; i<=11; i++ ) {
       float x = 6*i+0.5 ;
       l1 -> SetLineWidth(1) ;
       l1 -> SetLineStyle(2) ;
       l1 -> SetLineColor(4) ;
       l1 -> DrawLine( x, hmin, x, hmax ) ;
    }



 } // draw_boundaries_combine

//==========================




