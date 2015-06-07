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

   void label_plot( float ymax, bool logscale ) ;

  //---------

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
      Double_t* qcd_pred_yvals    = g_qcdbg -> GetY() ;
      Double_t* qcd_pred_xvals    = g_qcdbg -> GetX() ;
      Double_t* qcd_pred_yerrlow  = g_qcdbg -> GetEYlow() ;
      Double_t* qcd_pred_yerrhigh = g_qcdbg -> GetEYhigh() ;
      Double_t* qcd_pred_xerrlow  = g_qcdbg -> GetEXlow() ;
      Double_t* qcd_pred_xerrhigh = g_qcdbg -> GetEXhigh() ;



      TH1F* h_qcd_mc_diff_norm = new TH1F( "h_qcd_mc_diff_norm", "(value-pred)/allmc, QCD MC", 72, 0.5, 72.5 ) ;
      TH1F* h_qcd_mc_diff_norm_zeromc = new TH1F( "h_qcd_mc_diff_norm_zeromc", "(value-pred)/allmc, QCD MC", 72, 0.5, 72.5 ) ;

      h_qcd_mc_diff_norm -> SetMarkerStyle( 20 ) ;
      h_qcd_mc_diff_norm_zeromc -> SetMarkerStyle( 24 ) ;

      double qcd_pred_yvals_zero[100] ;
      double qcd_pred_yerrhigh_norm[100] ;
      double qcd_pred_yerrlow_norm[100] ;
      for ( int i=0; i<72; i++ ) {
         qcd_pred_yvals_zero[i] = 0. ;
         float allbgmc = h_allbg_mc -> GetBinContent( i+1 ) ;
         if ( allbgmc > 0 ) {
            qcd_pred_yerrhigh_norm[i] = qcd_pred_yerrhigh[i] / allbgmc ;
            qcd_pred_yerrlow_norm[i]  = qcd_pred_yerrlow[i] / allbgmc ;
            float qcdmc_val = h_qcd_mc -> GetBinContent( i+1 ) ;
            float qcdmc_err = h_qcd_mc -> GetBinError( i+1 ) ;
            float diff_norm_val = (qcdmc_val - qcd_pred_yvals[i]) / allbgmc ;
            float diff_norm_err = qcdmc_err / allbgmc ;
            if ( qcdmc_val > 0 ) {
               h_qcd_mc_diff_norm -> SetBinContent( i+1, diff_norm_val ) ;
               h_qcd_mc_diff_norm -> SetBinError( i+1, diff_norm_err ) ;
               h_qcd_mc_diff_norm_zeromc -> SetBinContent( i+1, +9 ) ;
            } else {
               h_qcd_mc_diff_norm_zeromc -> SetBinContent( i+1, diff_norm_val ) ;
               h_qcd_mc_diff_norm_zeromc -> SetBinError( i+1, diff_norm_err ) ;
            }
         } else {
            qcd_pred_yerrhigh_norm[i] = 0. ;
            qcd_pred_yerrlow_norm[i]  = 0. ;
            h_qcd_mc_diff_norm_zeromc -> SetBinContent( i+1, +9 ) ;
         }
      } // i


      TGraphAsymmErrors* g_pred_error_norm = new TGraphAsymmErrors( 72, qcd_pred_xvals, qcd_pred_yvals_zero,
                                                                        qcd_pred_xerrlow, qcd_pred_xerrhigh,
                                                                        qcd_pred_yerrlow_norm, qcd_pred_yerrhigh_norm ) ;

      g_pred_error_norm -> SetTitle( "g_pred_error_norm" ) ;
      g_pred_error_norm -> SetFillColor( kRed-9 ) ;










      gDirectory -> cd("Rint:/") ;
      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_draw_bg_dphi_1dsplit_hists" ) ;
      if ( can == 0x0 ) {
         can = new TCanvas( "can_draw_bg_dphi_1dsplit_hists", "QCD 1D split", 1500, 1300 ) ;
      }
      can -> Clear() ;
      can -> Divide(1,2) ;










      TH1F* h_qcd_pred = new TH1F( "h_qcd_pred", "QCD BG model", 72, 0.5, 72.5 ) ;
      for ( int i=0; i<72; i++ ) {  h_qcd_pred -> SetBinContent( i+1, qcd_pred_yvals[i] ) ; } ;
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

         float ymax = h_qcd_mc -> GetMaximum() ;
         label_plot( ymax, do_logy ) ;

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

    //-----------

         gPad -> SetLogy(0) ;

         float yrange = 0.6 ;
         h_qcd_mc_diff_norm -> SetMaximum( yrange ) ;
         h_qcd_mc_diff_norm -> SetMinimum(-1.*yrange ) ;
         h_qcd_mc_diff_norm -> SetMarkerStyle(20) ;
         h_qcd_mc_diff_norm_zeromc -> SetMarkerStyle(24) ;
         h_qcd_mc_diff_norm -> SetLineWidth(1) ;

         h_qcd_mc_diff_norm -> SetYTitle( "( value - prediction )/ (all BG, MC)" ) ;

         can -> cd(2) ;
         h_qcd_mc_diff_norm -> DrawCopy( ) ;
         TLine l ;
         l.DrawLine(0.5, 0., 72.5, 0.) ;
         g_pred_error_norm -> Draw( "2" ) ;
         h_qcd_mc_diff_norm_zeromc -> DrawCopy("same P0") ;
         h_qcd_mc_diff_norm -> DrawCopy( "same" ) ;
         h_qcd_mc_diff_norm -> DrawCopy("axis same") ;
         h_qcd_mc_diff_norm -> DrawCopy("axig same") ;

         label_plot( yrange, false ) ;
         gPad -> SetGridy(1) ;

         { // norm legend
            can -> cd(0) ;
            float padyh(1.) ;
            if ( show_events && show_norm ) { padyh = 0.5 ; }
            TPad* norm_pad = new TPad( "norm_pad", "norm_pad", 0.80, 0.00, 1.0, padyh ) ;
            norm_pad -> Range( 0., 0., 1., 1. ) ;
            norm_pad -> Draw() ;
            norm_pad -> cd() ;
            TBox* tb = new TBox() ;
            tb -> SetFillStyle(1001) ;
            tb -> SetFillColor( kRed-9 ) ;
            tb -> DrawBox( 0.10, 0.80, 0.20, 0.90 ) ;
            TText* ltx = new TText() ;
            ltx -> SetTextAlign( 12 ) ;
            ltx -> SetTextSize( 0.06 ) ;
            ltx -> DrawText( 0.25, 0.87, "QCD BG prediction error" ) ;
            ltx -> DrawText( 0.25, 0.82, "divided by all BG, MC" ) ;
            TEllipse* lte = new TEllipse() ;
            lte -> SetFillStyle(1001) ;
            lte -> SetFillColor(1) ;
            lte -> DrawEllipse( 0.15, 0.40, 0.02, 0.02*0.2/0.5, 0, 360, 0 ) ;
            TLine* ll = new TLine() ;
            ll -> DrawLine( 0.10, 0.40, 0.20, 0.40 ) ;
            ll -> DrawLine( 0.15, 0.33, 0.15, 0.47 ) ;
            ltx -> DrawText( 0.25, 0.42, "QCD MC" ) ;
            ltx -> DrawText( 0.25, 0.37, "divided by all BG, MC" ) ;
         }





        fp_mchist -> Close() ;
        fp_lhscan -> Close() ;


   } // qcd_closure_plot_lhscan


  //===============================================================
   void label_plot( float ymax, bool logscale ) {

      TText* ttext = new TText() ;
      TLine* line = new TLine() ;

      float txty ;

      line -> SetLineStyle(1) ;
      line -> SetLineWidth(1) ;
      line -> SetLineColor( kBlue ) ;
      ttext -> SetTextColor( kBlue ) ;
      ttext -> SetTextAlign(22) ;
      ttext -> SetTextSize( 0.030 ) ;
      for ( int i=1; i<=12; i++ ) {
         line -> DrawLine( 6*i+0.5, -0.85*ymax, 6*i+0.5, 0.87*ymax ) ;
         char text[100] ;
         sprintf( text, "Nb%d", (i-1)%4 ) ;
         if ( logscale ) { txty = 0.5*0.84*ymax ; } else { txty = 0.84*ymax ; }
         ttext -> DrawText( 6*i-2.5, txty, text ) ;
      }
      line -> SetLineStyle(1) ;
      line -> SetLineColor( 1 ) ;
      line -> DrawLine( 24.5, -ymax, 24.5, ymax ) ;
      line -> DrawLine( 48.5, -ymax, 48.5, ymax ) ;
      ttext -> SetTextColor( kRed ) ;
      ttext -> SetTextSize( 0.038 ) ;

      if ( logscale ) { txty = 2.02*ymax ; } else { txty = 1.10*ymax ; }
      ttext -> DrawText( 12.5, txty, "Njet1" ) ;
      ttext -> DrawText( 36.5, txty, "Njet2" ) ;
      ttext -> DrawText( 60.5, txty, "Njet3" ) ;


   }

  //===============================================================




