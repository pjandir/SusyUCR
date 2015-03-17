

#include "draw_qcd_dphi_all_ht_met4.c"

#include "TLegend.h"
#include "TSystem.h"



   TH1F* shift_hist( TH1F* hp, float shift ) ;

   //-----

   void draw_qcd_ratio_summary_nb_njet4( bool use_dphin = true, float mdp_cut = 4.0, float ratio_max = 0.9, bool create_hists = false, bool fill_hists = false ) {

      gSystem -> Exec( " mkdir -p outputfiles " ) ;

      int njetbins = 5; //3 bins->4, 4bins->5

      char mindphi_var[100] ;
      if ( use_dphin ) {
         sprintf( mindphi_var, "minDeltaPhiN") ;
      } else {
         sprintf( mindphi_var, "minDeltaPhi") ;
      }

      if ( create_hists ) {

         draw_qcd_dphi_all_ht_met4( use_dphin, 0, mdp_cut, fill_hists, 1 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 1, mdp_cut, fill_hists, 1 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 2, mdp_cut, fill_hists, 1 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 3, mdp_cut, fill_hists, 1 ) ;

         draw_qcd_dphi_all_ht_met4( use_dphin, 0, mdp_cut, fill_hists, 2 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 1, mdp_cut, fill_hists, 2 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 2, mdp_cut, fill_hists, 2 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 3, mdp_cut, fill_hists, 2 ) ;

         draw_qcd_dphi_all_ht_met4( use_dphin, 0, mdp_cut, fill_hists, 3 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 1, mdp_cut, fill_hists, 3 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 2, mdp_cut, fill_hists, 3 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 3, mdp_cut, fill_hists, 3 ) ;

         draw_qcd_dphi_all_ht_met4( use_dphin, 0, mdp_cut, fill_hists, 4 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 1, mdp_cut, fill_hists, 4 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 2, mdp_cut, fill_hists, 4 ) ;
         //draw_qcd_dphi_all_ht_met4( use_dphin, 3, mdp_cut, fill_hists, 4 ) ;

         if ( njetbins == 5 )
           draw_qcd_dphi_all_ht_met4( use_dphin, 0, mdp_cut, fill_hists, 5 ) ;

         gDirectory -> ls( "h_ratio*" ) ;

      }

      gStyle -> SetPadBottomMargin( 0.13 ) ;
      gStyle -> SetPadLeftMargin( 0.13 ) ;
      gStyle -> SetPadRightMargin( 0.05 ) ;
      gStyle -> SetPadTopMargin( 0.05 ) ;
      char cname[1000] ;
      sprintf( cname, "can1ratio_%s", mindphi_var ) ;
      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( cname ) ;
      if ( can1 == 0x0 ) {
         char ctitle[1000] ;
         sprintf( ctitle, "H/L Ratio, %s", mindphi_var ) ;
         can1 = new TCanvas( cname, ctitle, 1800, 800 ) ;
      }
      can1 -> Clear() ;
      can1 -> Divide(4,1) ;



      for ( int njbi=1; njbi<=njetbins; njbi++ ) {

         can1 -> cd( njbi ) ;


         char hname[1000] ;
         TH1F *hp, *hps ;

         sprintf( hname, "h_ratio_%s_nb0_njet%d", mindphi_var, njbi ) ;
         hp = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hp == 0x0 ) { printf( "\n\n *** Missing hist: %s\n\n", hname ) ; return ; }
         hp -> SetMaximum( ratio_max ) ;
/*
         hps = shift_hist( hp, 0.0 ) ;
         TH1F* hp_nb0 = hps ;

         sprintf( hname, "h_ratio_%s_nb1_njet%d", mindphi_var, njbi ) ;
         hp = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hp == 0x0 ) { printf( "\n\n *** Missing hist: %s\n\n", hname ) ; return ; }
         hp -> SetMaximum( ratio_max ) ;
         hps = shift_hist( hp, 0.1 ) ;
         TH1F* hp_nb1 = hps ;

         sprintf( hname, "h_ratio_%s_nb2_njet%d", mindphi_var, njbi ) ;
         hp = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hp == 0x0 ) { printf( "\n\n *** Missing hist: %s\n\n", hname ) ; return ; }
         hp -> SetMaximum( ratio_max ) ;
         hps = shift_hist( hp, 0.2 ) ;
         TH1F* hp_nb2 = hps ;

         sprintf( hname, "h_ratio_%s_nb3_njet%d", mindphi_var, njbi ) ;
         hp = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hp == 0x0 ) { printf( "\n\n *** Missing hist: %s\n\n", hname ) ; return ; }
         hp -> SetMaximum( ratio_max ) ;
         hps = shift_hist( hp, 0.3 ) ;
         TH1F* hp_nb3 = hps ;




         hp_nb0 -> SetLineColor(1) ;
         hp_nb1 -> SetLineColor(2) ;
         hp_nb2 -> SetLineColor(4) ;
         hp_nb3 -> SetLineColor(6) ;

         hp_nb0 -> SetMarkerColor(1) ;
         hp_nb1 -> SetMarkerColor(2) ;
         hp_nb2 -> SetMarkerColor(4) ;
         hp_nb3 -> SetMarkerColor(6) ;

         TLegend* legend = new TLegend( 0.22, 0.71,  0.45, 0.87 ) ;
         legend -> SetFillColor( kWhite ) ;
         legend -> SetTextSize(0.055) ;
         legend -> AddEntry( hp_nb0, "= 0 b" ) ;
         legend -> AddEntry( hp_nb1, "= 1 b" ) ;
         legend -> AddEntry( hp_nb2, "= 2 b" ) ;
         legend -> AddEntry( hp_nb3, ">= 3 b" ) ;

         char title[100] ;
         sprintf( title, "%s   H/L Ratio", mindphi_var ) ;
         hp_nb3 -> SetYTitle( title ) ;
         hp_nb3 -> SetTitleOffset( 1.6, "y" ) ;
         hp_nb3 -> SetMaximum( ratio_max ) ;

         hp_nb3 -> Draw() ;
         hp_nb2 -> Draw("same") ;
         hp_nb1 -> Draw("same") ;
         hp_nb0 -> Draw("same") ;
*/
         hp -> SetLineColor(1);
         hp -> SetMarkerColor(1);
         TLegend* legend = new TLegend( 0.22, 0.71,  0.45, 0.87 ) ;
         legend -> SetFillColor( kWhite ) ;
         legend -> SetTextSize(0.055) ;
         legend -> AddEntry( hp, ">= 0 b" ) ;
         char title[100] ;
         sprintf( title, "%s   H/L Ratio", mindphi_var ) ;
         hp -> SetYTitle( title ) ;
         hp -> SetTitleOffset( 1.6, "y" ) ;
         hp -> SetMaximum( ratio_max ) ;
         hp -> Draw();

         legend -> Draw() ;

         gPad -> SetGridy(1) ;

      } // njbi


      can1 -> Update() ; can1 -> Draw() ;
      char fname[10000] ;
      sprintf( fname, "outputfiles/%s-ratio-summary-cut%4.2f-nb-njet.pdf", mindphi_var, mdp_cut ) ;
      can1 -> SaveAs( fname ) ;


      sprintf( fname, "outputfiles/%s-ratio-histograms-cut%4.2f-nb-njet.root", mindphi_var, mdp_cut ) ;
      printf("\n\n Saving histograms in %s\n\n", fname ) ;
      saveHist( fname, "h*" ) ;

   } // draw_qcd_ratio_summary_nb_njet4

 //=================================================================================================

   TH1F* shift_hist( TH1F* hp, float shift ) {

      if ( hp == 0x0 ) return 0x0 ;

      char hname[100] ;
      sprintf( hname, "%s_shift", hp -> GetName() ) ;

      int nbins = hp -> GetNbinsX() ;
      float xl = hp -> GetXaxis() -> GetXmin() + shift ;
      float xh = hp -> GetXaxis() -> GetXmax() + shift ;

      TH1F* rp = new TH1F( hname, "", nbins, xl, xh ) ;
      for ( int bi=1; bi<nbins; bi++ ) {
         rp -> SetBinContent( bi, hp -> GetBinContent( bi ) ) ;
         rp -> SetBinError( bi, hp -> GetBinError( bi ) ) ;
         rp -> GetXaxis() -> SetBinLabel( bi, hp -> GetXaxis() -> GetBinLabel( bi ) ) ;
      } // bi

      rp -> GetXaxis() -> LabelsOption( "v" ) ;

      rp -> SetLineWidth( hp -> GetLineWidth() ) ;
      rp -> SetLineColor( hp -> GetLineColor() ) ;
      rp -> SetMarkerStyle( hp -> GetMarkerStyle() ) ;
      rp -> SetMarkerColor( hp -> GetMarkerColor() ) ;
      rp -> SetMaximum( hp -> GetMaximum() ) ;

      return rp ;

   } // shift_hist

 //=================================================================================================




