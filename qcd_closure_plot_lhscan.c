#include "TSystem.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TText.h"
#include "THStack.h"
#include "TLine.h"
#include "TStyle.h"
#include "TPad.h"
#include "TPave.h"
#include "TBox.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TEllipse.h"


   void qcd_closure_plot_lhscan( bool do_logy = true, int zoom_option = 0,
                                 bool show_events = true, bool show_norm = true,
                                 const char* infile_mchist = "outputfiles/fill-bg-hists4-sb-allinone-postdraw.root",
                                 const char* infile_lhscan = "outputfiles/per-bin-lh-analysis-v4b-t1bbbbH-nonperfect-closure-ss1-fullfit.root" ) {

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin(0.30) ;
      gStyle -> SetPadRightMargin(0.20) ;

      char dirname[1000] ;

      gDirectory -> Delete( "h*" ) ;

      TFile* fp_mchist(0x0) ;
      TFile* fp_lhscan(0x0) ;


      sprintf( dirname, "%s:/", infile_mchist ) ;
      printf("\n\n Checking %s\n", dirname ) ;
      bool cdstat(false) ;
      cdstat = gDirectory -> cd( dirname ) ;
      if ( !cdstat ) {
         printf("\n\n === Opening file : %s\n\n", infile_mchist ) ;
         fp_mchist = new TFile( infile_mchist, "READ" ) ;
         gDirectory -> pwd() ;
      } else {
         printf("\n\n File %s already open.\n\n", infile_mchist ) ;
         fp_mchist = gFile ;
         gDirectory -> pwd() ;
         gDirectory -> ls() ;
         if ( fp_mchist == 0x0 ) { printf("\n\n *** null file pointer.\n\n" ) ; return; }
      }

      TH1F* h_qcd_mc = (TH1F*) fp_mchist -> Get( "h_1d_allinone_zl_qcd" ) ;
      if ( h_qcd_mc == 0x0 ) { printf("\n\n *** can't find h_1d_allinone_zl_qcd.\n\n") ; return ; }

      TH1F* h_allbg_mc = (TH1F*) fp_mchist -> Get( "h_1d_allinone_zl_bgsum" ) ;
      if ( h_allbg_mc == 0x0 ) { printf("\n\n *** can't find h_1d_allinone_zl_allbg.\n\n") ; return ; }




      sprintf( dirname, "%s:/", infile_lhscan ) ;
      cdstat = false ;
      printf("\n\n Checking %s\n", dirname ) ;
      cdstat = gDirectory -> cd( dirname ) ;
      if ( !cdstat ) {
         printf("\n\n === Opening file : %s\n\n", infile_lhscan ) ;
         fp_lhscan = new TFile( infile_lhscan, "READ" ) ;
         gDirectory -> pwd() ;
      } else {
         printf("\n\n File %s already open.\n\n", infile_lhscan ) ;
         fp_lhscan = gFile ;
         gDirectory -> pwd() ;
         gDirectory -> ls() ;
         if ( fp_lhscan == 0x0 ) { printf("\n\n *** null file pointer.\n\n" ) ; return; }
      }

      TGraphAsymmErrors* g_qcdbg = (TGraphAsymmErrors*) fp_lhscan -> Get( "g_qcdbg" ) ;
      if ( g_qcdbg == 0x0 ) { printf("\n\n *** can't find g_qcdbg.\n\n") ; return ; }




      gDirectory -> cd("Rint:/") ;
      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_draw_bg_dphi_1dsplit_hists" ) ;
      if ( can == 0x0 ) {
         can = new TCanvas( "can_draw_bg_dphi_1dsplit_hists", "QCD 1D split", 1500, 1300 ) ;
      }
      can -> Clear() ;
      can -> Divide(1,2) ;










      TH1F* h_qcd_pred = new TH1F( "h_qcd_pred", "QCD BG model", 72, 0.5, 72.5 ) ;
      Double_t* yvals = g_qcdbg -> GetY() ;
      for ( int i=0; i<72; i++ ) {  h_qcd_pred -> SetBinContent( i+1, yvals[i] ) ; } ;
      TH1F* h_qcd_pred2 = (TH1F*) h_qcd_pred -> Clone( "h_qcd_pred2" ) ;

      TH1F* h_allbg_mc2 = new TH1F( "h_allbg_mc2", "All BG, MC", 72, 0.5, 72.5 ) ;
      for ( int i=0; i<72; i++ ) {  h_allbg_mc2 -> SetBinContent( i+1, h_allbg_mc -> GetBinContent( i+1 ) ) ; } ;



      g_qcdbg -> SetFillColor( kRed-9 ) ;
      h_qcd_pred -> SetFillColor( kRed-10 ) ;
      h_qcd_mc -> SetMarkerStyle(20) ;
      h_allbg_mc2 -> SetLineColor( 0 ) ;
      h_allbg_mc2 -> SetFillColor( kBlue-10 ) ;


      if ( do_logy ) {
         h_qcd_mc -> SetMaximum( 5.0 * ( h_qcd_mc -> GetMaximum() ) ) ;
         h_qcd_mc -> SetMinimum(0.01 ) ;
      } else {
         h_qcd_mc -> SetMinimum( -1 ) ;
         if ( zoom_option == 1 ) {
            h_qcd_mc -> SetMaximum( 500. ) ;
         } else if ( zoom_option == 2 ) {
            h_qcd_mc -> SetMaximum( 100. ) ;
         } else if ( zoom_option == 3 ) {
            h_qcd_mc -> SetMaximum( 25. ) ;
         } else if ( zoom_option == 4 ) {
            h_qcd_mc -> SetMaximum( 10. ) ;
         } else if ( zoom_option == 5 ) {
            h_qcd_mc -> SetMaximum( 3. ) ;
         }
      }






      can -> cd(1) ;

      h_qcd_mc -> DrawCopy( ) ;
      if ( do_logy ) gPad -> SetLogy() ;
      h_allbg_mc2 -> DrawCopy( "hist same" ) ;
      h_qcd_pred -> DrawCopy( "hist same" ) ;
      g_qcdbg -> Draw( "2 same" ) ;
      h_qcd_mc -> DrawCopy( "same" ) ;
      h_qcd_pred2 -> DrawCopy( "hist same" ) ;
      h_qcd_mc -> DrawCopy("axis same") ;

         { // events legend
            can -> cd(0) ;
            float padyl(0.) ;
            if ( show_events && show_norm ) { padyl = 0.5 ; }
            TPad* evts_pad = new TPad( "evts_pad", "evts_pad", 0.80, padyl, 1.0, 1.0 ) ;
            evts_pad -> Range( 0., 0., 1., 1. ) ;
            evts_pad -> Draw() ;
            evts_pad -> cd() ;
            TBox* tb = new TBox( 0.10, 0.80, 0.20, 0.90 ) ;
            tb -> SetFillStyle(1001) ;
            tb -> SetFillColor( kBlue-10 ) ;
            tb -> DrawBox( 0.10, 0.80, 0.20, 0.90 ) ;
            TText* ltx = new TText() ;
            ltx -> SetTextAlign( 12 ) ;
            ltx -> SetTextSize( 0.06 ) ;
            ltx -> DrawText( 0.25, 0.85, "All BG, MC" ) ;
            tb -> SetFillColor( kRed-10 ) ;
            tb -> DrawBox( 0.10, 0.50, 0.20, 0.70 ) ;
            tb -> SetFillColor( kRed-9 ) ;
            tb -> DrawBox( 0.10, 0.63, 0.20, 0.70 ) ;
            TLine* ll = new TLine() ;
            ll -> DrawLine( 0.10, 0.665, 0.20, 0.665 ) ;
            ltx -> DrawText( 0.25, 0.60, "QCD BG prediction" ) ;
            TEllipse* lte = new TEllipse() ;
            lte -> SetFillStyle(1001) ;
            lte -> SetFillColor(1) ;
            lte -> DrawEllipse( 0.15, 0.40, 0.02, 0.02*0.2/0.5, 0, 360, 0 ) ;
            ll -> DrawLine( 0.10, 0.40, 0.20, 0.40 ) ;
            ll -> DrawLine( 0.15, 0.33, 0.15, 0.47 ) ;
            ltx -> DrawText( 0.25, 0.40, "QCD MC" ) ;
         }



     fp_mchist -> Close() ;
     fp_lhscan -> Close() ;


   } // qcd_closure_plot_lhscan






