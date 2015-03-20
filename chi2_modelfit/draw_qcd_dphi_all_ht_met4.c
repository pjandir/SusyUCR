
#include "histio.c"
#include "TSystem.h"
#include "draw_qcd_dphi_plot4.c"


   void draw_qcd_dphi_all_ht_met4( bool use_dphin = true, int nb_bin = 0, float mdp_cut = 4.0, bool fill_plots = false, int njet_bin = -1 ) {

      printf("\n\n In draw_qcd_dphi_all_ht_met4: nb=%d, nj=%d\n\n", nb_bin, njet_bin ) ;  std::cout.flush() ;

      gStyle -> SetPadBottomMargin(0.20) ;
      gStyle -> SetOptStat(0) ;

      TText text ;
      text.SetTextSize(0.060) ;

      TLine line ;
      char mindphi_var[100] ;

      if ( use_dphin ) {
         sprintf( mindphi_var, "minDeltaPhiN") ;
      } else {
         sprintf( mindphi_var, "minDeltaPhi") ;
      }

      char cname[100] ;
      sprintf( cname, "can_nb%d_%s", nb_bin, mindphi_var ) ;
      TCanvas* can1htmet = (TCanvas*) gDirectory -> FindObject( cname ) ;
      if ( can1htmet == 0x0 ) {
         char ctitle[100] ;
         sprintf( ctitle, "%s, Nb%d", mindphi_var, nb_bin ) ;
         printf("\n Creating canvas %s, %s\n", cname, ctitle ) ; std::cout.flush() ;
         can1htmet = new TCanvas( cname, ctitle, 1000, 900 ) ;
         can1htmet -> SetWindowPosition(50,50) ;
      }
      can1htmet -> Clear() ;
      can1htmet -> Divide(3,3) ;
      int ci(1) ;

  //  printf("\n\n Deleting histograms matching h_only_good_stats_qcdbins_nb*\n\n") ; std::cout.flush() ;
  //  gDirectory -> Delete( "h_only_good_stats_qcdbins_nb*" ) ;
  //  printf("\n\n Just after delete.\n\n") ; std::cout.flush() ;


      char njetbin_str[100] ;
      char njetbin_title_str[100] ;
      char njetbin_binlabel_str[100] ;
      char njetbin_fname_str[100] ;
      if ( njet_bin >= 0 ) {
         sprintf( njetbin_str, "_njet%d", njet_bin ) ;
         sprintf( njetbin_title_str, ", Njet%d", njet_bin ) ;
         sprintf( njetbin_binlabel_str, " Njet%d", njet_bin ) ;
         sprintf( njetbin_fname_str, "-njet%d", njet_bin ) ;
      } else {
         sprintf( njetbin_str, "" ) ;
         sprintf( njetbin_title_str, "" ) ;
         sprintf( njetbin_binlabel_str, "" ) ;
         sprintf( njetbin_fname_str, "" ) ;
      }

      char hfname[10000] ;
      sprintf( hfname, "outputfiles/%s-nb%d-histograms%s.root", mindphi_var, nb_bin, njetbin_fname_str ) ;

      if ( !fill_plots ) {
         printf("\n\n Loading histograms from %s\n\n", hfname ) ; std::cout.flush() ;
         ///////////////// loadHist( hfname ) ;
         char hpat[1000] ;
         sprintf( hpat, "h_qcdbins_ge300_nb%d*%s", nb_bin, njetbin_str ) ;
         printf( "\n Loading histograms matching %s\n", hpat ) ;
         loadHist( hfname, 0, hpat ) ;
         sprintf( hpat, "h_only_good_stats_qcdbins_nb%d*%s", nb_bin, njetbin_str ) ;
         printf( "\n Loading histograms matching %s\n", hpat ) ;
         loadHist( hfname, 0, hpat ) ;
      }


      int n_analysis_bins_ht(3) ;
      int n_analysis_bins_met(4) ;

      char hname[100] ;
      char htitle[1000] ;

      sprintf( hname, "h_rms_%s_nb%d%s", mindphi_var, nb_bin, njetbin_str ) ;
      sprintf( htitle, "RMS(0 mean), %s, Nb%d%s", mindphi_var, nb_bin, njetbin_title_str ) ;
      TH1F* h_rms = new TH1F( hname, htitle, n_analysis_bins_met*(n_analysis_bins_ht+1)+1, 0.5, n_analysis_bins_met*(n_analysis_bins_ht+1)+1.5 ) ;

      sprintf( hname, "h_ratio_%s_nb%d%s", mindphi_var, nb_bin, njetbin_str ) ;
      sprintf( htitle, "H/L Ratio, %s, Nb%d, >%4.2f%s", mindphi_var, nb_bin, mdp_cut, njetbin_title_str ) ;
      TH1F* h_ratio = new TH1F( hname, htitle, n_analysis_bins_met*(n_analysis_bins_ht+1)+1, 0.5, n_analysis_bins_met*(n_analysis_bins_ht+1)+1.5 ) ;

      sprintf( hname, "h_qcd_yield_all_nb%d_%s%s", nb_bin, mindphi_var, njetbin_str ) ;
      sprintf( htitle, "QCD yield, Nb%d, %s>%4.2f%s", nb_bin, mindphi_var, mdp_cut, njetbin_title_str ) ;
      TH1F* h_yield_all = new TH1F( hname, htitle, n_analysis_bins_met*(n_analysis_bins_ht+1)+1, 0.5, n_analysis_bins_met*(n_analysis_bins_ht+1)+1.5 ) ;

      sprintf( hname, "h_qcd_yield_pass_nb%d_%s%s", nb_bin, mindphi_var, njetbin_str ) ;
      sprintf( htitle, "QCD yield, Nb%d, %s>%4.2f%s", nb_bin, mindphi_var, mdp_cut, njetbin_title_str ) ;
      TH1F* h_yield_pass = new TH1F( hname, htitle, n_analysis_bins_met*(n_analysis_bins_ht+1)+1, 0.5, n_analysis_bins_met*(n_analysis_bins_ht+1)+1.5 ) ;


      int summary_bi(2) ;

         for ( int metbi=1; metbi<=n_analysis_bins_met; metbi++ ) {
      for ( int htbi=1; htbi<=n_analysis_bins_ht; htbi++ ) {


            printf("  MET%d, HT%d\n", metbi, htbi ) ; std::cout.flush() ;

            if ( fill_plots ) draw_qcd_dphi_plot4( use_dphin, nb_bin, htbi, metbi, mdp_cut, njet_bin ) ;

            sprintf( hname, "h_only_good_stats_qcdbins_nb%d_ht%d_met%d%s", nb_bin, htbi, metbi, njetbin_str ) ;
            TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
            if ( hp == 0x0 ) { printf("\n\n *** Missing histogram: %s\n\n", hname ) ; return ; }

            sprintf( hname, "h_qcdbins_ge300_nb%d_ht%d_met%d%s", nb_bin, htbi, metbi, njetbin_str ) ;
            TH1F* hp_qcdbinsge300 = (TH1F*) gDirectory -> FindObject( hname ) ;
            if ( hp_qcdbinsge300 == 0x0 ) { printf("\n\n *** Missing histogram: %s\n\n", hname ) ; return ; }

            can1htmet -> cd(ci) ;
            hp -> SetXTitle( mindphi_var ) ;
            hp -> SetTitleSize( 0.050, "x" ) ;
            hp -> SetFillColor(21) ;
            hp -> DrawCopy( "hist" ) ;
            hp -> DrawCopy( "same" ) ;

            double zm_rms_val, zm_rms_err ;
            calc_0m_rms( hp, zm_rms_val, zm_rms_err ) ;

            float rval, rerr ;
            get_ratio( hp, mdp_cut, rval, rerr ) ;

            int cut_bin = hp -> GetXaxis() -> FindBin( mdp_cut ) ;
            int max_bin = hp -> GetNbinsX() ;
            double events_above_cut_val, events_above_cut_err ;
            events_above_cut_val = hp -> IntegralAndError( cut_bin, max_bin, events_above_cut_err ) ;
            double events_all_val, events_all_err ;
            events_all_val = hp -> IntegralAndError( 1, max_bin, events_all_err ) ;

            double events_above_cut_qcdbinsge300_val, events_above_cut_qcdbinsge300_err ;
            events_above_cut_qcdbinsge300_val = hp_qcdbinsge300 -> IntegralAndError( cut_bin, max_bin, events_above_cut_qcdbinsge300_err ) ;
            double events_all_qcdbinsge300_val, events_all_qcdbinsge300_err ;
            events_all_qcdbinsge300_val = hp_qcdbinsge300 -> IntegralAndError( 1, max_bin, events_all_qcdbinsge300_err ) ;

            printf( "\n =========== Nb%d, HT%d, MET%d%s : QCD high/low ratio : %5.3f +/- %5.3f\n\n", nb_bin, htbi, metbi, njetbin_title_str, rval, rerr ) ;

            char ratio_text[100] ;
            sprintf( ratio_text, "Ratio: %5.3f +/- %5.3f\n", rval, rerr ) ;
            text.DrawTextNDC( 0.30, 0.80, ratio_text ) ;

            char zmrms_text[100] ;
            sprintf( zmrms_text, "RMS0: %.3f +/- %.3f\n", zm_rms_val, zm_rms_err ) ;
            text.DrawTextNDC( 0.30, 0.70, zmrms_text ) ;

            char events_text[100] ;
            sprintf( events_text, "evts. pass: %3.1f +/- %3.1f", events_above_cut_val, events_above_cut_err ) ;
            text.DrawTextNDC( 0.30, 0.60, events_text ) ;


            char bin_label[100] ;
            sprintf( bin_label, "HT%d MET%d%s", htbi, metbi, njetbin_binlabel_str ) ;

            h_rms -> SetBinContent( summary_bi, zm_rms_val ) ;
            h_rms -> SetBinError( summary_bi, zm_rms_err ) ;
            h_rms -> GetXaxis() -> SetBinLabel( summary_bi, bin_label ) ;

            h_ratio -> SetBinContent( summary_bi, rval ) ;
            h_ratio -> SetBinError( summary_bi, rerr ) ;
            h_ratio -> GetXaxis() -> SetBinLabel( summary_bi, bin_label ) ;

            //------
            // h_yield_all  -> SetBinContent( summary_bi, events_all_val ) ;
            // h_yield_pass -> SetBinContent( summary_bi, events_above_cut_val ) ;
            // h_yield_all  -> SetBinError( summary_bi, events_all_err ) ;
            // h_yield_pass -> SetBinError( summary_bi, events_above_cut_err ) ;
            //------
            h_yield_all  -> SetBinContent( summary_bi, events_all_qcdbinsge300_val ) ;
            h_yield_pass -> SetBinContent( summary_bi, events_above_cut_qcdbinsge300_val ) ;
            h_yield_all  -> SetBinError( summary_bi, events_all_qcdbinsge300_err ) ;
            h_yield_pass -> SetBinError( summary_bi, events_above_cut_qcdbinsge300_err ) ;
            //------

            h_yield_all  -> GetXaxis() -> SetBinLabel( summary_bi, bin_label ) ;
            h_yield_pass -> GetXaxis() -> SetBinLabel( summary_bi, bin_label ) ;

            line.SetLineColor(2) ;
            line.SetLineStyle(2) ;
            line.DrawLine( mdp_cut, 0., mdp_cut, hp->GetMaximum() ) ;

            line.SetLineColor(4) ;
            line.SetLineStyle(3) ;
            line.DrawLine( 0., 0., hp -> GetXaxis()-> GetXmax(), 0. ) ;

            can1htmet -> Update() ; can1htmet -> Draw() ;


            ci++ ;

            summary_bi ++ ;

         } // metbi
         summary_bi ++ ;
      } // htbi

      char fname[10000] ;
      sprintf( fname, "outputfiles/%s-nb%d-all-ht-met-bins-cut%4.2f%s.pdf", mindphi_var, nb_bin, mdp_cut, njetbin_fname_str ) ;
      can1htmet -> SaveAs( fname ) ;


      line.SetLineWidth(2) ;

      sprintf( cname, "can_rms_nb%d_%s", nb_bin, mindphi_var ) ;
      TCanvas* can2htmet = (TCanvas*) gDirectory -> FindObject( cname ) ;
      if ( can2htmet == 0x0 ) {
         can2htmet = new TCanvas( cname, "RMS(0 mean)", 700, 600 ) ;
         can2htmet -> SetWindowPosition(1050,150) ;
      }
      can2htmet -> Clear() ;
      can2htmet -> cd() ;
      h_rms -> GetXaxis() -> LabelsOption( "v" ) ;
      h_rms -> SetLineWidth(2) ;
      h_rms -> SetMarkerStyle(20) ;
      if ( use_dphin ) {
         h_rms -> SetMaximum( 9.0 ) ;
      } else {
         h_rms -> SetMaximum( 0.3 ) ;
      }
      h_rms -> DrawCopy() ;
      line.SetLineColor(4) ;
      line.SetLineStyle(3) ;
      line.DrawLine( h_rms -> GetXaxis() -> GetXmin(), 0., h_rms -> GetXaxis() -> GetXmax(), 0. ) ;
      gPad -> SetGridy(1) ;
      sprintf( fname, "outputfiles/%s-nb%d-rms0%s.pdf", mindphi_var, nb_bin, njetbin_fname_str ) ;
      can2htmet -> SaveAs( fname ) ;




      sprintf( cname, "can_ratio_nb%d_%s", nb_bin, mindphi_var ) ;
      TCanvas* can3htmet = (TCanvas*) gDirectory -> FindObject( cname ) ;
      if ( can3htmet == 0x0 ) {
         can3htmet = new TCanvas( cname, "H/L Ratio", 700, 600 ) ;
         can3htmet -> SetWindowPosition(1750,250) ;
      }
      can3htmet -> Clear() ;
      can3htmet -> cd() ;
      h_ratio -> GetXaxis() -> LabelsOption( "v" ) ;
      h_ratio -> SetLineWidth(2) ;
      h_ratio -> SetMarkerStyle(20) ;
      if ( use_dphin ) {
         h_ratio -> SetMaximum( 0.8 ) ;
      } else {
         h_ratio -> SetMinimum( -0.02 ) ;
         h_ratio -> SetMaximum( 0.15 ) ;
      }
      h_ratio -> DrawCopy() ;
      line.SetLineColor(4) ;
      line.SetLineStyle(3) ;
      line.DrawLine( h_ratio -> GetXaxis() -> GetXmin(), 0., h_ratio -> GetXaxis() -> GetXmax(), 0. ) ;
      gPad -> SetGridy(1) ;
      sprintf( fname, "outputfiles/%s-nb%d-ratio-cut%4.2f%s.pdf", mindphi_var, nb_bin, mdp_cut, njetbin_fname_str ) ;
      can3htmet -> SaveAs( fname ) ;





      h_yield_all  -> GetXaxis() -> LabelsOption( "v" ) ;
      h_yield_pass -> GetXaxis() -> LabelsOption( "v" ) ;

      sprintf( cname, "can_nb%d_%s_susy_yield", nb_bin, mindphi_var ) ;
      TCanvas* can4 = (TCanvas*) gDirectory -> FindObject( cname ) ;
      if ( can4 == 0x0 ) {
         char ctitle[1000] ;
         sprintf( ctitle, "QCD yield, %s, Nb%d", mindphi_var, nb_bin ) ;
         can4 = new TCanvas( cname, ctitle, 700, 500 ) ;
      }
      ////can4 -> SetWindowPosition(1800,50) ;
      can4 -> Clear() ;
      can4 -> cd() ;

      h_yield_all -> SetLineColor(3) ;
      h_yield_pass -> SetLineColor(4) ;
      h_yield_all -> SetLineWidth(2) ;
      h_yield_pass -> SetLineWidth(2) ;
      h_yield_all -> SetMinimum(0.5) ;
      h_yield_all -> Draw() ;
      h_yield_all -> Draw("hist same") ;
      h_yield_pass -> Draw("same") ;
      h_yield_pass -> Draw("hist same") ;
      gPad -> SetGridy(1) ;
      gPad -> SetLogy(1) ;
      sprintf( fname, "outputfiles/%s-nb%d-yield-cut%4.2f-qcd%s.pdf", mindphi_var, nb_bin, mdp_cut, njetbin_fname_str ) ;
      can4 -> SaveAs( fname ) ;



      printf("\n\n draw_qcd_dphi_all_ht_met4 : saving histograms in %s\n\n", hfname ) ; std::cout.flush() ;

      if ( fill_plots ) {
         char hpat[1000] ;
         sprintf( hpat, "h_*_nb%d_*%s", nb_bin, njetbin_str ) ;
         printf( "\n Saving histograms matching %s\n", hpat ) ;
         gSystem -> Exec( "mkdir -p outputfiles" ) ;
         saveHist( hfname, hpat ) ;
      }

      printf("\n\n End of draw_qcd_dphi_all_ht_met4.\n\n") ; std::cout.flush() ;


   } // draw_qcd_dphi_all_ht_met4




