

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
#include "TEllipse.h"

#include "histio.c"


   //float mdp_cut[4] = { 0.35, 0.25, 0.20, 0.20 } ;
   ///float mdp_cut[4] = { 0.30, 0.30, 0.30, 0.30 } ;

   float mdp_cut ;

   TFile* p_hist_file ;

   TH1F* h_dphi_zl_qcd ;
   TH1F* h_dphi_zl_lostlep ;
   TH1F* h_dphi_zl_znunu ;
   TH1F* h_dphi_zl_sig ;

   TH1F* h_dphi_sl_qcd ;
   TH1F* h_dphi_sl_lostlep ;
   TH1F* h_dphi_sl_znunu ;
   TH1F* h_dphi_sl_sig ;

   TPad* text_pad ;
     char text_line[1000] ;
   TText* ttext ;

     float text_y_ndc(0.95) ;
     float text_x_ndc(0.0) ;
     float text_size(0.015) ;
     float text_dy_ndc(0.017) ;

   int sb_mht_ht_ind ;
   int sb_nj_ind ;

   int fb_mht_ind[10] ;
   int fb_ht_ind[10] ;
   int fb_nj_ind[10] ;
   int fb_nb_ind ;
   int n_fb_in_sb ;

   double n_zl_qcd_val[10] ;
   double n_zl_qcd_err[10] ;
   double n_zl_lostlep_val[10] ;
   double n_zl_lostlep_err[10] ;
   double n_zl_znunu_val[10] ;
   double n_zl_znunu_err[10] ;
   double n_zl_sig_val[10] ;
   double n_zl_sig_err[10] ;
   double n_zl_val[10] ;

   double n_ldp_qcd_val[10] ;
   double n_ldp_qcd_err[10] ;
   double n_ldp_lostlep_val[10] ;
   double n_ldp_lostlep_err[10] ;
   double n_ldp_znunu_val[10] ;
   double n_ldp_znunu_err[10] ;
   double n_ldp_val[10] ;

   double n_slldp_qcd_val[10] ;
   double n_slldp_qcd_err[10] ;
   double n_slldp_lostlep_val[10] ;
   double n_slldp_lostlep_err[10] ;
   double n_slldp_znunu_val[10] ;
   double n_slldp_znunu_err[10] ;
   double n_slldp_val[10] ;

   double n_sl_qcd_val[10] ;
   double n_sl_qcd_err[10] ;
   double n_sl_lostlep_val[10] ;
   double n_sl_lostlep_err[10] ;
   double n_sl_znunu_val[10] ;
   double n_sl_znunu_err[10] ;
   double n_sl_val[10] ;

   char sb_name[100] ;

   char signal_name[100] ;

   char varname[100] ;

  //---------
   TH1F* get_hist( const char* hname ) ;
   void fill_fb_index_arrays( ) ;
   void collect_dphi_hists() ;
   void calc_qcd( float& val, float& err ) ;
   float calc_qcd_err() ;
   void fill_kqcd_arrays() ;
   void label_plot( float ymax, bool logscale ) ;
   void print_table( TH1* hp_qcd_pred, TH1* hp_qcd_mc, TH1* hp_allbg, TGraph* gr_qcdbg_lhs ) ;
  //---------

   float Kqcd_ht_val[3] ;
   float Kqcd_ht_err[3] ;

   float Kqcd_mht_val[4] ;
   float Kqcd_mht_err[4] ;

   float Kqcd_nj_val[5] ;
   float Kqcd_nj_err[5] ;


  //===============================================================

   void qcd_closure_plot( bool do_logy = true,
                          int zoom_option = 0,
                          bool show_events = true,
                          bool show_norm = true,
                          bool norm_by_allmc = true,
                          float arg_mdp_cut = 6.0,
                          const char* arg_varname = "mdpn",
                          const char* infile =  "outputfiles/fill-bg-dphi-hists2-t1bbbbH-postdraw.root",
                          const char* signame = "t1bbbbH" ) {

      if ( !show_events && !show_norm ) return ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadBottomMargin(0.30) ;
      gStyle -> SetPadRightMargin(0.20) ;

      gDirectory -> cd("Rint:/") ;
      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_qcd_closure_plot" ) ;
      int canypx(650) ;
      if ( show_events && show_norm ) canypx = 2*canypx ;
      if ( can == 0x0 ) {
         can = new TCanvas( "can_qcd_closure_plot", "QCD closure plot", 1500, canypx ) ;
      }
      can -> SetCanvasSize( 1500, canypx ) ;
      can -> Clear() ;
      if ( show_events && show_norm ) {
         can -> Divide(1,2) ;
      }

      gDirectory -> Delete( "h*") ;
      char dirname[1000] ;
      sprintf( dirname, "%s:/", infile ) ;
      bool cdstat(false) ;
      cdstat = gDirectory -> cd( dirname ) ;
      if ( !cdstat ) {
         printf("\n\n === Opening file : %s\n\n", infile ) ;
         p_hist_file = new TFile( infile, "READ" ) ;
         gDirectory -> pwd() ;
      } else {
         printf("\n\n File %s already open.\n\n", infile ) ;
         p_hist_file = gFile ;
         if ( p_hist_file == 0x0 ) { printf("\n\n *** null file pointer.\n\n" ) ; return; }
      }


      sprintf( varname, "%s", arg_varname ) ;
      sprintf( signal_name, "%s", signame ) ;

      float shift = 0.06 ;
      TH1F* h_qcd_pred = new TH1F( "h_qcd_pred","", 72, 0.5-shift, 72.5-shift ) ;
      TH1F* h_qcd_mc   = new TH1F( "h_qcd_mc","", 72, 0.5+shift, 72.5+shift ) ;
      TH1F* h_qcd_pred_err_norm = new TH1F( "h_qcd_pred_err_norm","", 72, 0.5-shift, 72.5-shift ) ;
      TH1F* h_qcd_mc_diff_norm   = new TH1F( "h_qcd_mc_diff_norm","", 72, 0.5+shift, 72.5+shift ) ;
      TH1F* h_qcd_mc_diff_norm_zeromc   = new TH1F( "h_qcd_mc_diff_norm_zeromc","", 72, 0.5+shift, 72.5+shift ) ;
      TH1F* h_qcd_hl_ratio_ave = new TH1F( "h_qcd_hl_ratio_ave", "", 72, 0.5, 72.5 ) ;

      TH1F* h_allbg_mc   = new TH1F( "h_allbg_mc","", 72, 0.5+shift, 72.5+shift ) ;


      int hbi(0) ;

      mdp_cut = arg_mdp_cut ;

      for ( int nji=1; nji<=3; nji++ ) {
         for ( int nbi=0; nbi<=3; nbi++ ) {
            for ( int mhthti=1; mhthti<=6; mhthti++ ) {

               hbi++ ;

               sb_mht_ht_ind = mhthti ;
               sb_nj_ind = nji ;

               fb_nb_ind = nbi ;

               fill_kqcd_arrays() ;
               fill_fb_index_arrays( ) ;
               collect_dphi_hists() ;
               printf(" ======= %s\n", sb_name ) ;
               float pred_zl_qcd_val, pred_zl_qcd_err ;
               calc_qcd( pred_zl_qcd_val, pred_zl_qcd_err ) ;

               int nbins = h_dphi_zl_qcd -> GetNbinsX() ;
               int cutbin = h_dphi_zl_qcd -> FindBin( mdp_cut+1e-5 ) ;



               double sb_n_zl_qcd_err(0.) ;
               float sb_n_zl_qcd_val = h_dphi_zl_qcd -> IntegralAndError( cutbin, nbins, sb_n_zl_qcd_err ) ;
               double sb_n_zl_lostlep_err(0.) ;
               float sb_n_zl_lostlep_val = h_dphi_zl_lostlep -> IntegralAndError( cutbin, nbins, sb_n_zl_lostlep_err ) ;
               double sb_n_zl_znunu_err(0.) ;
               float sb_n_zl_znunu_val = h_dphi_zl_znunu -> IntegralAndError( cutbin, nbins, sb_n_zl_znunu_err ) ;

               float sb_n_zl_val = sb_n_zl_qcd_val + sb_n_zl_lostlep_val + sb_n_zl_znunu_val ;




               double sb_n_ldp_qcd_err(0.) ;
               float sb_n_ldp_qcd_val = h_dphi_zl_qcd -> IntegralAndError( 1, cutbin-1, sb_n_ldp_qcd_err ) ;
               double sb_n_ldp_lostlep_err(0.) ;
               float sb_n_ldp_lostlep_val = h_dphi_zl_lostlep -> IntegralAndError( 1, cutbin-1, sb_n_ldp_lostlep_err ) ;
               double sb_n_ldp_znunu_err(0.) ;
               float sb_n_ldp_znunu_val = h_dphi_zl_znunu -> IntegralAndError( 1, cutbin-1, sb_n_ldp_znunu_err ) ;

               float sb_n_ldp_val = sb_n_ldp_qcd_val + sb_n_ldp_lostlep_val + sb_n_ldp_znunu_val ;


               float ave_qcd_hl_ratio_val(0.) ;
               float ave_qcd_hl_ratio_err(0.) ;
               if ( sb_n_ldp_qcd_val > 0 ) {
                  ave_qcd_hl_ratio_val = pred_zl_qcd_val / sb_n_ldp_qcd_val ;
                  ave_qcd_hl_ratio_err = pred_zl_qcd_err / sb_n_ldp_qcd_val ;
               }

               h_qcd_hl_ratio_ave -> SetBinContent( hbi, ave_qcd_hl_ratio_val ) ;
               h_qcd_hl_ratio_ave -> SetBinError( hbi, ave_qcd_hl_ratio_err ) ;


               h_qcd_pred -> SetBinContent( hbi, pred_zl_qcd_val ) ;
               h_qcd_pred -> SetBinError( hbi, pred_zl_qcd_err ) ;

               h_qcd_pred_err_norm -> SetBinContent( hbi, 0. ) ;
               float err_norm(0.) ;
               if ( norm_by_allmc ) {
                  if ( sb_n_zl_val > 0 ) { err_norm = pred_zl_qcd_err / sb_n_zl_val ; }
               } else {
                  if ( pred_zl_qcd_val > 0 ) { err_norm = pred_zl_qcd_err / pred_zl_qcd_val ; }
               }
               h_qcd_pred_err_norm -> SetBinError( hbi, err_norm ) ;


               h_qcd_mc -> SetBinContent( hbi, sb_n_zl_qcd_val ) ;
               h_qcd_mc -> SetBinError( hbi, sb_n_zl_qcd_err ) ;

               float mc_diff_norm_val = 0. ;
               float mc_diff_norm_err = 0. ;
               if ( norm_by_allmc ) {
                  if ( sb_n_zl_val > 0 ) {
                     mc_diff_norm_val = ( sb_n_zl_qcd_val - pred_zl_qcd_val ) / sb_n_zl_val ;
                     mc_diff_norm_err = sb_n_zl_qcd_err / sb_n_zl_val ;
                  }
               } else {
                  if ( pred_zl_qcd_val > 0 && sb_n_zl_qcd_val > 0. ) {
                     mc_diff_norm_val = ( sb_n_zl_qcd_val - pred_zl_qcd_val ) / pred_zl_qcd_val ;
                     mc_diff_norm_err = sb_n_zl_qcd_err / pred_zl_qcd_val ;
                  }
               }
               if ( sb_n_zl_qcd_val > 0.001 ) {
                  h_qcd_mc_diff_norm -> SetBinContent( hbi, mc_diff_norm_val ) ;
                  h_qcd_mc_diff_norm -> SetBinError( hbi, mc_diff_norm_err ) ;
                  h_qcd_mc_diff_norm_zeromc -> SetBinContent( hbi, -99 ) ;
               } else {
                  printf("zero MC.  diff: %f +/- %f\n", mc_diff_norm_val, mc_diff_norm_err ) ;
                  h_qcd_mc_diff_norm_zeromc -> SetBinContent( hbi, mc_diff_norm_val+1e-6 ) ;
                  h_qcd_mc_diff_norm_zeromc -> SetBinError( hbi, mc_diff_norm_err ) ;
               }


               printf(" %3d : %s pred err norm %.3f\n", hbi, sb_name, err_norm ) ;

               h_qcd_mc -> GetXaxis() -> SetBinLabel( hbi, sb_name ) ;
               h_qcd_pred_err_norm -> GetXaxis() -> SetBinLabel( hbi, sb_name ) ;
               h_qcd_hl_ratio_ave -> GetXaxis() -> SetBinLabel( hbi, sb_name ) ;



               h_allbg_mc -> SetBinContent( hbi, sb_n_zl_val ) ;




               printf("           MC ZL QCD val : %.3f +/- %.3f\n", sb_n_zl_qcd_val, sb_n_zl_qcd_err ) ;


            } // mhthti
         } // nbi
      } // nji


      if ( show_events ) {

         gPad -> SetLogy(0) ;

         TH1F* h_qcd_pred2 = (TH1F*) h_qcd_pred -> Clone( "h_qcd_pred2" ) ;
         TH1F* h_qcd_pred3 = (TH1F*) h_qcd_pred -> Clone( "h_qcd_pred3" ) ;

         h_qcd_pred  -> SetFillColor( kRed-9 ) ;
         h_qcd_pred2 -> SetFillColor( kRed-10 ) ;
         h_qcd_pred  -> SetLineColor(0) ;
         h_qcd_pred2 -> SetLineColor(0) ;
         h_qcd_pred3 -> SetLineColor(1) ;
         h_qcd_mc -> SetLineWidth(1) ;
         h_qcd_mc -> SetMarkerStyle(20) ;
         //h_allbg_mc -> SetLineColor( 4 ) ;
         h_allbg_mc -> SetLineColor( 0 ) ;
         h_allbg_mc -> SetFillColor( kBlue-10 ) ;


         h_qcd_mc -> GetXaxis() -> LabelsOption("v") ;
         h_qcd_pred_err_norm -> GetXaxis() -> LabelsOption("v") ;

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

         h_qcd_mc -> SetYTitle( "Events" ) ;

         if ( show_events && show_norm ) can -> cd(1) ;
         h_qcd_mc -> Draw() ;
         if ( do_logy ) gPad -> SetLogy() ;
         h_allbg_mc -> Draw( "same" ) ;
         h_qcd_pred2 -> Draw( "hist same" ) ;
         h_qcd_pred -> Draw( "e2 same" ) ;
         h_qcd_pred3 -> Draw( "hist same" ) ;
         h_qcd_mc -> Draw( "same" ) ;
         h_qcd_mc -> Draw( "axis same" ) ;

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


      } // show_events ?

    //-----

      if ( show_norm ) {

         gPad -> SetLogy(0) ;

         h_qcd_pred_err_norm -> SetFillColor( kRed-9 ) ;
         float yrange ;
         if ( norm_by_allmc ) { yrange = 0.6 ; } else { yrange = 4 ; }
         h_qcd_pred_err_norm -> SetMaximum( yrange ) ;
         h_qcd_pred_err_norm -> SetMinimum(-1.*yrange ) ;
         h_qcd_mc_diff_norm -> SetMarkerStyle(20) ;
         h_qcd_mc_diff_norm_zeromc -> SetMarkerStyle(24) ;
         h_qcd_mc_diff_norm -> SetLineWidth(1) ;

         if ( norm_by_allmc ) {
           h_qcd_pred_err_norm -> SetYTitle( "( value - prediction )/ (all BG, MC)" ) ;
         } else {
           h_qcd_pred_err_norm -> SetYTitle( "( value - prediction )/ (qcd BG, MC)" ) ;
         }

         if ( show_events && show_norm ) can -> cd(2) ;
         h_qcd_pred_err_norm -> Draw( "e2" ) ;
         TLine l ;
         l.DrawLine(0.5, 0., 72.5, 0.) ;
         h_qcd_mc_diff_norm -> Draw("same") ;
         h_qcd_mc_diff_norm_zeromc -> Draw("same P0") ;
         h_qcd_mc_diff_norm -> Draw("axis same") ;
         h_qcd_mc_diff_norm -> Draw("axig same") ;

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

      } // show_norm?

      char logzoomstr[100] ;
      if ( do_logy ) { sprintf( logzoomstr, "logy" ) ; } else { sprintf( logzoomstr, "zoom%d", zoom_option ) ; }
      char evtsnormstr[100] ;
      if ( show_events && show_norm ) {
         sprintf( evtsnormstr, "both-evts-and-norm" ) ;
      } else {
         if ( show_events ) sprintf( evtsnormstr, "evts" ) ;
         if ( show_norm ) sprintf( evtsnormstr, "norm" ) ;
      }
      char savefile[1000] ;
      if ( show_events ) {
         sprintf( savefile, "outputfiles/qcd-closure-%s-%.2f-%s-%s.pdf", varname, mdp_cut, logzoomstr, evtsnormstr ) ;
      } else {
         sprintf( savefile, "outputfiles/qcd-closure-%s-%.2f-%s.pdf", varname, mdp_cut, evtsnormstr ) ;
      }
      can -> SaveAs( savefile ) ;



      TFile* lhscan_file = new TFile( "outputfiles/per-bin-lh-analysis-v4b-t1bbbbH-nonperfect-closure-ss1-fullfit.root", "READ" ) ;
      TGraph* g_qcdbg(0x0) ;
      if ( lhscan_file != 0 ) {
         g_qcdbg = (TGraph*) lhscan_file -> Get( "g_qcdbg" ) ;
         if ( g_qcdbg != 0 ) {
            double* y ;
            double* yel ;
            double* yeh ;
            y = g_qcdbg -> GetY() ;
            yel = g_qcdbg -> GetEYlow() ;
            yeh = g_qcdbg -> GetEYhigh() ;
            for ( int bi = 1; bi<=72; bi++ ) {
               printf(" %20s : %7.2f +/- %5.2f |  %7.2f -%.2f +%.2f\n",
                  h_qcd_mc -> GetXaxis() -> GetBinLabel( bi ) ,
                  h_qcd_pred -> GetBinContent( bi ) , h_qcd_pred -> GetBinError( bi ) ,
                  y[bi-1], yel[bi-1], yeh[bi-1] ) ;

            }
         }
      }

      print_table( h_qcd_pred, h_qcd_mc, h_allbg_mc, g_qcdbg ) ;

      //////h_qcd_mc -> Print("all") ;


   } // qcd_closure_plot


  //===============================================================

   TH1F* get_hist( const char* hname ) {

      TH1F* hp = (TH1F*) p_hist_file -> Get( hname ) ;
      if ( hp == 0x0 ) { printf("\n\n **** Missing histogram: %s\n\n", hname ) ; gSystem->Exit(-1) ; }

      return hp ;

   } // get_hist

  //===============================================================

   void fill_fb_index_arrays( ) {

      if ( sb_mht_ht_ind == 1 ) {

         sprintf( sb_name, "HT1_MHT1_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 2 ) {

         sprintf( sb_name, "HT2_MHT1_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 3 ) {

         sprintf( sb_name, "HT3_MHT1_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 4 ) {

         sprintf( sb_name, "HT12_MHT2_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 5 ) {

         sprintf( sb_name, "HT3_MHT2_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 6 ) {

         sprintf( sb_name, "HT23_MHT3_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      }

   } // fill_fb_index_arrays

  //===============================================================

   void collect_dphi_hists() {

     char hname[1000] ;
     char newhname[1000] ;

     for ( int fbi=0; fbi<n_fb_in_sb; fbi++ ) {

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_qcd_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_qcd = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_qcd_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_qcd = (TH1F*) h_zl_qcd -> Clone( newhname ) ; } else { h_dphi_zl_qcd -> Add( h_zl_qcd ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_lostlep_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_lostlep = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_lostlep_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_lostlep = (TH1F*) h_zl_lostlep -> Clone( newhname ) ; } else { h_dphi_zl_lostlep -> Add( h_zl_lostlep ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_znunu_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_znunu = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_znunu_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_znunu = (TH1F*) h_zl_znunu -> Clone( newhname ) ; } else { h_dphi_zl_znunu -> Add( h_zl_znunu ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_%s_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], signal_name, fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_sig = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_sig_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_sig = (TH1F*) h_zl_sig -> Clone( newhname ) ; } else { h_dphi_zl_sig -> Add( h_zl_sig ) ; }


        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_sl_qcd_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_sl_qcd = get_hist( hname ) ;
        sprintf( newhname, "h_%s_sl_qcd_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_sl_qcd = (TH1F*) h_sl_qcd -> Clone( newhname ) ; } else { h_dphi_sl_qcd -> Add( h_sl_qcd ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_sl_lostlep_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_sl_lostlep = get_hist( hname ) ;
        sprintf( newhname, "h_%s_sl_lostlep_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_sl_lostlep = (TH1F*) h_sl_lostlep -> Clone( newhname ) ; } else { h_dphi_sl_lostlep -> Add( h_sl_lostlep ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_sl_znunu_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_sl_znunu = get_hist( hname ) ;
        sprintf( newhname, "h_%s_sl_znunu_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_sl_znunu = (TH1F*) h_sl_znunu -> Clone( newhname ) ; } else { h_dphi_sl_znunu -> Add( h_sl_znunu ) ; }


      //-------

        int nbins = h_zl_qcd -> GetNbinsX() ;
        int cutbin = h_zl_qcd -> FindBin( mdp_cut+1e-5 ) ;

        n_zl_qcd_val[fbi] = h_zl_qcd -> IntegralAndError( cutbin, nbins, n_zl_qcd_err[fbi] ) ;
        n_zl_lostlep_val[fbi] = h_zl_lostlep -> IntegralAndError( cutbin, nbins, n_zl_lostlep_err[fbi] ) ;
        n_zl_znunu_val[fbi] = h_zl_znunu -> IntegralAndError( cutbin, nbins, n_zl_znunu_err[fbi] ) ;

        n_ldp_qcd_val[fbi] = h_zl_qcd -> IntegralAndError( 1, cutbin-1, n_ldp_qcd_err[fbi] ) ;
        n_ldp_lostlep_val[fbi] = h_zl_lostlep -> IntegralAndError( 1, cutbin-1, n_ldp_lostlep_err[fbi] ) ;
        n_ldp_znunu_val[fbi] = h_zl_znunu -> IntegralAndError( 1, cutbin-1, n_ldp_znunu_err[fbi] ) ;

        n_sl_qcd_val[fbi] = h_sl_qcd -> IntegralAndError( cutbin, nbins, n_sl_qcd_err[fbi] ) ;
        n_sl_lostlep_val[fbi] = h_sl_lostlep -> IntegralAndError( cutbin, nbins, n_sl_lostlep_err[fbi] ) ;
        n_sl_znunu_val[fbi] = h_sl_znunu -> IntegralAndError( cutbin, nbins, n_sl_znunu_err[fbi] ) ;

        n_slldp_qcd_val[fbi] = h_sl_qcd -> IntegralAndError( 1, cutbin-1, n_slldp_qcd_err[fbi] ) ;
        n_slldp_lostlep_val[fbi] = h_sl_lostlep -> IntegralAndError( 1, cutbin-1, n_slldp_lostlep_err[fbi] ) ;
        n_slldp_znunu_val[fbi] = h_sl_znunu -> IntegralAndError( 1, cutbin-1, n_slldp_znunu_err[fbi] ) ;

        n_sl_val[fbi]  = n_sl_lostlep_val[fbi] + n_sl_qcd_val[fbi] + n_sl_znunu_val[fbi] ;
        n_ldp_val[fbi] = n_ldp_lostlep_val[fbi] + n_ldp_qcd_val[fbi] + n_ldp_znunu_val[fbi] ;
        n_zl_val[fbi]  = n_zl_lostlep_val[fbi] + n_zl_qcd_val[fbi] + n_zl_znunu_val[fbi] ;

     } // fbi

     char htitle[1000] ;

     sprintf( htitle, "ZL QCD %s", sb_name ) ;
     h_dphi_zl_qcd -> SetTitle( htitle ) ;
     sprintf( htitle, "ZL LostLep %s", sb_name ) ;
     h_dphi_zl_lostlep -> SetTitle( htitle ) ;
     sprintf( htitle, "ZL Znunu %s", sb_name ) ;
     h_dphi_zl_znunu -> SetTitle( htitle ) ;
     sprintf( htitle, "ZL %s %s", signal_name, sb_name ) ;
     h_dphi_zl_sig -> SetTitle( htitle ) ;

     sprintf( htitle, "SL QCD %s", sb_name ) ;
     h_dphi_sl_qcd -> SetTitle( htitle ) ;
     sprintf( htitle, "SL LostLep %s", sb_name ) ;
     h_dphi_sl_lostlep -> SetTitle( htitle ) ;
     sprintf( htitle, "SL Znunu %s", sb_name ) ;
     h_dphi_sl_znunu -> SetTitle( htitle ) ;


   } // collect_dphi_hists

  //===============================================================

   void calc_qcd( float& pred_zl_qcd_val, float& pred_zl_qcd_err ) {

     pred_zl_qcd_val = 0. ;
     pred_zl_qcd_err = 0. ;

     bool verb(false) ;

     float total_pred_qcd_zl_val(0.) ;

     for ( int fbi=0; fbi<n_fb_in_sb; fbi++ ) {

        int mht_bin = fb_mht_ind[fbi] ;
        int ht_bin = fb_ht_ind[fbi] ;
        int nj_bin = fb_nj_ind[fbi] ;

        if (verb) printf("\n\n") ;
        sprintf( text_line, "====== Fine Bin: HT%d MHT%d Njet%d Nb%d", ht_bin, mht_bin, nj_bin, fb_nb_ind ) ;
        if (verb) printf( "%s\n", text_line ) ;
        text_y_ndc -= text_dy_ndc ;
        //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

        sprintf( text_line, "   lostlep : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f",
              n_sl_lostlep_val[fbi], n_sl_lostlep_err[fbi],
              n_slldp_lostlep_val[fbi], n_slldp_lostlep_err[fbi],
              n_ldp_lostlep_val[fbi], n_ldp_lostlep_err[fbi],
              n_zl_lostlep_val[fbi], n_zl_lostlep_err[fbi] ) ;
        text_y_ndc -= text_dy_ndc ;
        //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
        if (verb) printf( "%s\n", text_line ) ;

        sprintf( text_line, "   QCD     : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f",
              n_sl_qcd_val[fbi], n_sl_qcd_err[fbi],
              n_slldp_qcd_val[fbi], n_slldp_qcd_err[fbi],
              n_ldp_qcd_val[fbi], n_ldp_qcd_err[fbi],
              n_zl_qcd_val[fbi], n_zl_qcd_err[fbi] ) ;
        text_y_ndc -= text_dy_ndc ;
        //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
        if (verb) printf( "%s\n", text_line ) ;

        sprintf( text_line, "   Znunu   : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f",
              n_sl_znunu_val[fbi], n_sl_znunu_err[fbi],
              n_slldp_znunu_val[fbi], n_slldp_znunu_err[fbi],
              n_ldp_znunu_val[fbi], n_ldp_znunu_err[fbi],
              n_zl_znunu_val[fbi], n_zl_znunu_err[fbi] ) ;
        text_y_ndc -= text_dy_ndc ;
        //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
        if (verb) printf( "%s\n", text_line ) ;


        float Rqcd_val = Kqcd_ht_val[ht_bin-1] * Kqcd_mht_val[mht_bin-1] * Kqcd_nj_val[nj_bin-1] ;
        float Rqcd_err = Rqcd_val * sqrt( pow( Kqcd_ht_err[ht_bin-1]/Kqcd_ht_val[ht_bin-1], 2 )
                                        + pow( Kqcd_mht_err[mht_bin-1]/Kqcd_mht_val[mht_bin-1], 2 )
                                        + pow( Kqcd_nj_err[nj_bin-1]/Kqcd_nj_val[nj_bin-1], 2 ) ) ;

        if (verb) printf("\n") ;
        sprintf( text_line, "    Rqcd = Kqcd_ht * Kqcd_mht * Kqcd_nj = (%.3f +/- %.3f) * (%.3f +/- %.3f) * (%.3f +/- %.3f) = %.3f +/- %.3f",
              Kqcd_ht_val[ht_bin-1], Kqcd_ht_err[ht_bin-1],
              Kqcd_mht_val[mht_bin-1], Kqcd_mht_err[mht_bin-1],
              Kqcd_nj_val[nj_bin-1], Kqcd_nj_err[nj_bin-1],
              Rqcd_val, Rqcd_err ) ;
        text_y_ndc -= text_dy_ndc ;
        //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
        if (verb) printf( "%s\n", text_line ) ;


        float pred_qcd_zl_val = n_ldp_qcd_val[fbi] * Rqcd_val ;

        float pred_qcd_zl_err = 0. ;
        if ( n_ldp_qcd_val[fbi] > 0. ) {
           pred_qcd_zl_err = pred_qcd_zl_val * sqrt( pow( Rqcd_err / Rqcd_val, 2 ) + pow( sqrt(n_ldp_val[fbi])/n_ldp_qcd_val[fbi], 2 ) ) ;
        } else {
           pred_qcd_zl_err = pred_qcd_zl_val * sqrt( pow( Rqcd_err / Rqcd_val, 2 ) + n_ldp_val[fbi] ) ;
        }

        sprintf( text_line, "    Predicted ZL QCD : QCD_LDP * Rqcd_val = (%.2f +/- %.2f) * (%.3f +/- %.3f) =  %.2f +/- %.2f",
                  n_ldp_qcd_val[fbi], sqrt(n_ldp_val[fbi]),
                  Rqcd_val, Rqcd_err,
                  pred_qcd_zl_val, pred_qcd_zl_err ) ;
        text_y_ndc -= text_dy_ndc ;
        //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
        if (verb) printf( "%s\n", text_line ) ;

        float error_ratio(1.) ;
        if ( n_zl_val[fbi] > 0 ) error_ratio = pred_qcd_zl_err / sqrt( n_zl_val[fbi] ) ;

        if (verb) printf("   Error ratio ZL QCD / sqrt(N_ZL):  %.2f / %.2f = %.2f\n", pred_qcd_zl_err, sqrt( n_zl_val[fbi] ), error_ratio ) ;

        total_pred_qcd_zl_val += pred_qcd_zl_val ;

     } // fbi

    //-- Calculate the uncertainty, taking into account correlations between estimates for the fine bins.

     float total_pred_qcd_zl_err = calc_qcd_err() ;

     sprintf(text_line, "  Total predicted ZL QCD : %.3f +/- %.3f", total_pred_qcd_zl_val, total_pred_qcd_zl_err ) ;
     pred_zl_qcd_val = total_pred_qcd_zl_val ;
     pred_zl_qcd_err = total_pred_qcd_zl_err ;

     text_y_ndc -= text_dy_ndc ;
     //ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
     printf( "%s\n", text_line ) ;

   } // calc_qcd


  //===============================================================

   float calc_qcd_err() {

   //++++++++ The order matters for n_ldp array!  Check for consistency with fill_fb_index_arrays.
     float err2(0.) ;

     if ( sb_mht_ht_ind == 1 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[0] * Kqcd_nj_val[2] ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[0] * Kqcd_nj_val[2] ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_nj_val[2] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[0] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[3] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[0] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[4] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 2 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[1] * Kqcd_nj_val[2] ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[1] * Kqcd_nj_val[2] ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_nj_val[2] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[3] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[4] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 3 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[2] * Kqcd_nj_val[2] ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[2] * Kqcd_nj_val[2] ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_nj_val[2] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 4 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[1] * Kqcd_nj_val[2] ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[2] * Kqcd_nj_val[2] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[2] * Kqcd_nj_val[2] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[3] ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[3] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;


           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[4] ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[4] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;


           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 5 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[2] * Kqcd_nj_val[2]  ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_nj_val[2] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[2]   ,  2 ) ;


        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[3]   ,  2 ) ;


        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[4]   ,  2 ) ;


        }

     } else if ( sb_mht_ht_ind == 6 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[1] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[2] * Kqcd_nj_val[2]  ) * Kqcd_mht_err[3]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[3] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[3] * Kqcd_nj_val[2] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[3] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[3] * Kqcd_nj_val[2] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[2]   ,  2 ) ;



        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[3]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_nj_val[3] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[3]   ,  2 ) ;


        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[3]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_nj_val[4] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     }

     //-- LDP stats
     for ( int fbi=0; fbi<n_fb_in_sb; fbi++ ) {

      //-- stat error on N_ldp
        if ( n_ldp_val[fbi] > 0 ) {
           err2 += n_ldp_val[fbi] * pow( Kqcd_mht_val[fbi] * Kqcd_ht_val[fbi] * Kqcd_nj_val[fbi] , 2 ) ;
        } else {
           err2 += pow( Kqcd_mht_val[fbi] * Kqcd_ht_val[fbi] * Kqcd_nj_val[fbi] , 2 ) ;
        }

      //-- stat error on lost-lep subtraction.
        err2 += n_ldp_lostlep_val[fbi] * pow( Kqcd_mht_val[fbi] * Kqcd_ht_val[fbi] * Kqcd_nj_val[fbi] , 2 ) ;

     } // fbi

     return sqrt(err2) ;

   } // calc_qcd_err

  //===============================================================

   void fill_kqcd_arrays() {

      if ( strcmp( varname, "mdpn" ) == 0 ) {

      //-- these are for mdpN > 4.0
   ///// Kqcd_ht_val[0] = 0.268 ;  Kqcd_ht_err[0] = 0.015 ;
   ///// Kqcd_ht_val[1] = 0.200 ;  Kqcd_ht_err[1] = 0.009 ;
   ///// Kqcd_ht_val[2] = 0.096 ;  Kqcd_ht_err[2] = 0.007 ;

   ///// Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
   ///// Kqcd_mht_val[1] = 1.34 ; Kqcd_mht_err[1] = 0.1 ;
   ///// Kqcd_mht_val[2] = 2.54 ; Kqcd_mht_err[2] = 0.9 ;
   ///// Kqcd_mht_val[3] = 3.5  ; Kqcd_mht_err[3] = 1.5 ;

   ///// Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
   ///// Kqcd_nj_val[1] = 0.92 ; Kqcd_nj_err[1] = 0.06 ;
   ///// Kqcd_nj_val[2] = 0.62 ; Kqcd_nj_err[2] = 0.07 ;
   ///// Kqcd_nj_val[3] = 0.53 ; Kqcd_nj_err[3] = 0.11 ;
   ///// Kqcd_nj_val[4] = 0.64 ; Kqcd_nj_err[4] = 0.41 ;

      //-- these are for mdpN > 6.0, no isotrk veto
         Kqcd_ht_val[0] = 0.140 ;  Kqcd_ht_err[0] = 0.012 ;
         Kqcd_ht_val[1] = 0.096 ;  Kqcd_ht_err[1] = 0.007 ;
         Kqcd_ht_val[2] = 0.038 ;  Kqcd_ht_err[2] = 0.005 ;

         Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
         Kqcd_mht_val[1] = 1.29 ; Kqcd_mht_err[1] = 0.2 ;
         Kqcd_mht_val[2] = 2.09 ; Kqcd_mht_err[2] = 1.1 ;
         Kqcd_mht_val[3] = 3.3  ; Kqcd_mht_err[3] = 1.5 ;

         Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
         Kqcd_nj_val[1] = 0.69 ; Kqcd_nj_err[1] = 0.08 ;
         Kqcd_nj_val[2] = 0.69 ; Kqcd_nj_err[2] = 0.13 ;
         Kqcd_nj_val[3] = 0.50 ; Kqcd_nj_err[3] = 0.23 ;
         Kqcd_nj_val[4] = 0.50 ; Kqcd_nj_err[4] = 0.50 ;

  //  //-- these are for mdpN > 6.0, with aggressive isotrk veto
  //     Kqcd_ht_val[0] = 0.128 ;  Kqcd_ht_err[0] = 0.012 ;
  //     Kqcd_ht_val[1] = 0.089 ;  Kqcd_ht_err[1] = 0.007 ;
  //     Kqcd_ht_val[2] = 0.040 ;  Kqcd_ht_err[2] = 0.005 ;

  //     Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
  //     Kqcd_mht_val[1] = 1.00 ; Kqcd_mht_err[1] = 0.3 ;
  //     Kqcd_mht_val[2] = 2.10 ; Kqcd_mht_err[2] = 2.0 ;
  //     Kqcd_mht_val[3] = 3.0  ; Kqcd_mht_err[3] = 2.0 ;

  //     Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
  //     Kqcd_nj_val[1] = 1.17 ; Kqcd_nj_err[1] = 0.26 ;
  //     Kqcd_nj_val[2] = 0.30 ; Kqcd_nj_err[2] = 0.38 ;
  //     Kqcd_nj_val[3] = 0.31 ; Kqcd_nj_err[3] = 0.52 ;
  //     Kqcd_nj_val[4] = 0.31 ; Kqcd_nj_err[4] = 0.70 ;

      } else {

 ///  //-- these are for plain mdp > 0.3
 ///     Kqcd_ht_val[0] = 0.068 ;  Kqcd_ht_err[0] = 0.010 ;
 ///     Kqcd_ht_val[1] = 0.074 ;  Kqcd_ht_err[1] = 0.010 ;
 ///     Kqcd_ht_val[2] = 0.105 ;  Kqcd_ht_err[2] = 0.020 ;

 ///     Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
 ///     Kqcd_mht_val[1] = 0.28 ; Kqcd_mht_err[1] = 0.1 ;
 ///     Kqcd_mht_val[2] = 0.33 ; Kqcd_mht_err[2] = 0.3 ;
 ///     Kqcd_mht_val[3] = 0.61 ; Kqcd_mht_err[3] = 0.7 ;

 ///     Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
 ///     Kqcd_nj_val[1] = 1.38 ; Kqcd_nj_err[1] = 0.4 ;
 ///     Kqcd_nj_val[2] = 1.66 ; Kqcd_nj_err[2] = 1.0 ;
 ///     Kqcd_nj_val[3] = 1.8  ; Kqcd_nj_err[3] = 1.0 ;
 ///     Kqcd_nj_val[4] = 3.0  ; Kqcd_nj_err[4] = 2.5 ;

      //-- these are for plain mdp > 0.4
         Kqcd_ht_val[0] = 0.046 ;  Kqcd_ht_err[0] = 0.007 ;
         Kqcd_ht_val[1] = 0.041 ;  Kqcd_ht_err[1] = 0.004 ;
         Kqcd_ht_val[2] = 0.045 ;  Kqcd_ht_err[2] = 0.005 ;

         Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
         Kqcd_mht_val[1] = 0.36 ; Kqcd_mht_err[1] = 0.11 ;
         Kqcd_mht_val[2] = 0.43 ; Kqcd_mht_err[2] = 0.56 ;
         Kqcd_mht_val[3] = 0.86 ; Kqcd_mht_err[3] = 2.0 ;

         Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
         Kqcd_nj_val[1] = 1.89 ; Kqcd_nj_err[1] = 0.3 ;
         Kqcd_nj_val[2] = 2.02 ; Kqcd_nj_err[2] = 0.4 ;
         Kqcd_nj_val[3] = 2.66 ; Kqcd_nj_err[3] = 0.5 ;
         ////////Kqcd_nj_val[4] = 3.5  ; Kqcd_nj_err[4] = 2.5 ;
         Kqcd_nj_val[4] = 2.5  ; Kqcd_nj_err[4] = 1.0 ;

      }

   } // fill_kqcd_arrays

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

   void print_table( TH1* hp_qcd_pred, TH1* hp_qcd_mc, TH1* hp_allbg_mc, TGraph* gr_qcdbg_lhs ) {

      if ( hp_qcd_pred==0x0 || hp_qcd_mc==0x0 || hp_allbg_mc==0x0 || gr_qcdbg_lhs==0x0 ) {
          printf("\n\n *** missing input.\n\n") ;
          printf("  qcdpred %p   qcdmc %p   allbg %p   gr_qcdbg_lhs %p\n\n",
          hp_qcd_pred, hp_qcd_mc, hp_allbg_mc, gr_qcdbg_lhs ) ;
          return ;
      }

      if ( gr_qcdbg_lhs == 0x0 ) return ;
      double* y ;
      double* yel ;
      double* yeh ;
      y = gr_qcdbg_lhs -> GetY() ;
      yel = gr_qcdbg_lhs -> GetEYlow() ;
      yeh = gr_qcdbg_lhs -> GetEYhigh() ;

      int nbins = hp_qcd_pred -> GetNbinsX() ;

      printf("\n\n\n") ;
      printf("  search bin               QCD pred      QCD LHS    QCD MC        allbg MC\n") ;

      for ( int bi=1; bi<=nbins; bi++ ) {

         printf("  %20s :  %7.2f +/- %5.2f  |  %7.2f -%5.2f +%5.2f  |  %7.2f +/- %5.2f  |  %7.2f +/- %5.2f  |\n",
           hp_qcd_mc -> GetXaxis() -> GetBinLabel( bi ) ,
           hp_qcd_pred -> GetBinContent( bi ) ,  hp_qcd_pred -> GetBinError( bi ) ,
           y[bi-1], yel[bi-1], yeh[bi-1],
           hp_qcd_mc   -> GetBinContent( bi ) ,  hp_qcd_mc   -> GetBinError( bi ) ,
           hp_allbg_mc -> GetBinContent( bi ) ,  hp_allbg_mc -> GetBinError( bi ) 
           ) ;

      } // bi

      printf("\n\n\n") ;

   } // print_table

  //===============================================================















