
#include "TChain.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "calc_0m_rms.c"

   //-------
 ///==============================================
 //Nominal
/*
   int nhtbins = 4 ;   // first bin is not used in analysis.
   double htbins[5] = { 0., 500., 800., 1200., 20000. } ;  // first bin is not used in analysis.
   int nmetbins = 4 ;  // first bin is not used in analysis.
   double metbins[5] = { 0., 200., 500., 750., 20000. } ;  // first bin is not used in analysis.
   int nnjetbins = 4 ; // first bin is not used in analysis.
   double njetbins[5] = { 0., 3.5, 6.5, 8.5, 20. } ;  // first bin is not used in analysis.
*/
  //New nominal
   int nhtbins = 4 ;   // first bin is not used in analysis.
   double htbins[5] = { 0., 500., 800., 1200., 20000. } ;  // first bin is not used in analysis.
   int nmetbins = 5 ;  // first bin is not used in analysis.
   double metbins[6] = { 0., 200., 300., 500., 750., 20000. } ;  // first bin is not used in analysis.
   int nnjetbins = 6 ; // first bin is not used in analysis.
   double njetbins[7] = { 0., 3.5, 4.5, 5.5, 6.5, 8.5, 20. } ;  // first bin is not used in analysis.


   void add_overflow( TH1F* hp ) ;
   void get_ratio( TH1F* hp, float cut, float& rval, float& rerr ) ;

 //-------------

   void draw_qcd_dphi_plot4( bool use_dphin = true,
                             int nb_bin = 0,
                             int ht_bin = 4,
                             int met_bin = 1,
                             float mdp_cut = 4.0,
                             int njet_bin = -1 ) {

      gStyle -> SetTitleW(0.70) ;
      gStyle -> SetTitleH(0.08) ;

      TText text ;
      text.SetTextSize(0.065) ;

      TLine line ;


      char input_dir[100000] ;
      sprintf( input_dir, "/Users/pjandir/Desktop/qcd13/reducedTrees" ) ;
      //sprintf( input_dir, "/Users/pjandir/cms/ra2b/qcd13/reducedTrees" ) ;
      //sprintf( input_dir, "/Users/owen/work/cms/ra2b-2015/reducedTree-skim-nov12-2014-v73-3" ) ;

      printf("\n\n reducedTree directory : %s\n\n", input_dir ) ;

      char conditions[100] ;
      TString which = "v77";
      sprintf( conditions, "PU20bx25*%s" , which.Data()) ;

      ///////////// gDirectory->Delete("h*") ;

      int nsamples(0) ;
      char sname[50][100] ;

      //sprintf( sname[nsamples], "30to50" ) ; nsamples ++ ;
      //sprintf( sname[nsamples], "50to80" ) ; nsamples ++ ;
      //sprintf( sname[nsamples], "80to120" ) ; nsamples ++ ;
      //sprintf( sname[nsamples], "120to170" ) ; nsamples ++ ;
      //sprintf( sname[nsamples], "170to300" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "300to470" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "470to600" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "600to800" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "800to1000" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "1000to1400" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "1400to1800" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "1800to2400" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "2400to3200" ) ; nsamples ++ ;
      sprintf( sname[nsamples], "3200" ) ; nsamples ++ ;

      TChain* sch[nsamples] ;
      for ( int si=0; si<nsamples; si++ ) {
         char file_pattern[100000] ;
         int n_added ;
         int n_entries ;
         sch[si] = new TChain( "reducedTree" ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-%s*%s*.root", input_dir, sname[si], conditions ) ;
         cout << file_pattern << endl;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf("  Sample %20s : %3d files, %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
      }



      TH1F* h_mdp_all[nsamples] ;


      char mindphi_var[100] ;
      int nbins ;
      float xl, xh ;
      if ( use_dphin ) {
         sprintf( mindphi_var, "minDeltaPhiN") ;
         nbins = 200 ;
         xl = 0. ;
         xh = 100. ;
      } else {
         sprintf( mindphi_var, "minDeltaPhi") ;
         nbins = 200 ;
         xl = 0. ;
         xh = 4.0 ;
      }

      char nb_cut[1000] ;
      char ht_cut[1000] ;
      char met_cut[1000] ;
      char njet_cut[1000] ;
      //if ( nb_bin<3 ) {
      //  sprintf( nb_cut, "nbjets30==%d", nb_bin ) ; 
      //} else {
      //  sprintf( nb_cut, "nbjets30>=%d", nb_bin ) ; 
      //} 
      sprintf( nb_cut, "nbjets30>=%d", nb_bin ) ;
      sprintf( ht_cut , "HT30>%.0f && HT30<=%.0f", htbins[ht_bin], htbins[ht_bin+1] ) ;
      //sprintf( ht_cut , "HT>%.0f && HT<=%.0f", htbins[ht_bin], htbins[ht_bin+1] ) ;
      sprintf( met_cut, "MET>%.0f && MET<=%.0f", metbins[met_bin], metbins[met_bin+1] ) ;
      if ( njet_bin >= 0 ) {
         sprintf( njet_cut, "njets30>%.1f && njets30<=%.1f", njetbins[njet_bin], njetbins[njet_bin+1] ) ;
         //sprintf( njet_cut, "njets>%.1f && njets<=%.1f", njetbins[njet_bin], njetbins[njet_bin+1] ) ;
      } else {
         sprintf( njet_cut, "njets30>=0" ) ;
         //sprintf( njet_cut, "njets>=0" ) ;
      }
      printf(" Nb cut: %s\n", nb_cut ) ;
      printf(" HT bin %d,  cut: %s\n", ht_bin, ht_cut ) ;
      printf(" MET bin %d,  cut: %s\n", met_bin, met_cut ) ;
      printf(" Njet bin %d,  cut: %s\n", njet_bin, njet_cut ) ;

      char njetbin_str[100] ;
      if ( njet_bin >= 0 ) {
         sprintf( njetbin_str, "_njet%d", njet_bin ) ;
      } else {
         sprintf( njetbin_str, "" ) ;
      }

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "mdpcan1" ) ;
      if ( can1 == 0x0 ) can1 = new TCanvas( "mdpcan1", "Filling canvas", 700, 500 ) ;

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "mdpcan2" ) ;
      if ( can2 == 0x0 ) can2 = new TCanvas( "mdpcan2", "mdp, QCD bin samples", 900, 800 ) ;
      can2 -> Clear() ; can2 -> Divide(4,4) ;

      can2 -> SetWindowPosition(50,50) ;
      can1 -> SetWindowPosition(950,50) ;

      for ( int si=0; si<nsamples; si++ ) {

         char hname[100] ;
         char htitle[1000] ;
         char arg1[1000] ;
         char allcuts[10000] ;


         sprintf( hname, "h_mdp_all_%s_nb%d_ht%d_met%d%s", sname[si], nb_bin, ht_bin, met_bin, njetbin_str ) ;
         sprintf( htitle, "%s, all, %s", mindphi_var, sname[si] ) ;
         h_mdp_all[si] = new TH1F( hname, htitle,   nbins, xl, xh  ) ;
         h_mdp_all[si] -> Sumw2() ;
         TString mdpvar = mindphi_var;
         mdpvar += "_MHT";
         //if ( mdpvar == "minDeltaPhiN" ) mdpvar += "_pt30";
         sprintf( arg1, "%s>>%s", mdpvar.Data(), hname ) ;
         sprintf( allcuts, "(nElectrons==0&&nMuons==0&&(%s)&&(%s)&&(%s)&&(%s))*weight3*5000", nb_cut, ht_cut, met_cut, njet_cut ) ;
         can1 -> cd() ;
         sch[si] -> Draw( arg1, allcuts ) ; can1 -> Update() ; can1 -> Draw() ;
         add_overflow( h_mdp_all[si] ) ;
         double hmax = h_mdp_all[si] -> GetMaximum() ;
         h_mdp_all[si] -> SetMaximum( 1.2*hmax ) ;
         h_mdp_all[si] -> SetMinimum( -0.1*hmax ) ;

         can2 -> cd(si+1) ;
         h_mdp_all[si] -> Draw() ;


         float rval, rerr ;
         get_ratio( h_mdp_all[si], mdp_cut, rval, rerr ) ;

         printf( "%30s : %5.3f +/- %5.3f\n", hname, rval, rerr ) ;
         char ratio_text[100] ;
         sprintf( ratio_text, "%5.3f +/- %5.3f\n", rval, rerr ) ;

         text.DrawTextNDC( 0.30, 0.80, ratio_text ) ;

         line.SetLineColor(2) ;
         line.SetLineStyle(2) ;
         line.DrawLine( mdp_cut, 0., mdp_cut, h_mdp_all[si]->GetMaximum() ) ;

         line.SetLineColor(4) ;
         line.SetLineStyle(3) ;
         line.DrawLine( 0., 0., h_mdp_all[si]->GetXaxis() -> GetXmax(), 0. ) ;

         can2 -> Update() ; can2 -> Draw() ;

      } // si




    //----------

      TH1F* h_all_qcdbins ;
      {
         for ( int si=0; si<nsamples; si++ ) {
            if ( si==0 ) {
               char hname[100] ;
               sprintf( hname, "h_all_qcdbins_nb%d_ht%d_met%d%s", nb_bin, ht_bin, met_bin, njetbin_str ) ;
               h_all_qcdbins = (TH1F*) h_mdp_all[si]->Clone( hname ) ;
               char htitle[1000] ;
               if ( njet_bin >= 0 ) {
                  sprintf( htitle, "%s, Njet%d, Nb%d, HT%d, MET%d, all QCD bins", mindphi_var, njet_bin, nb_bin, ht_bin, met_bin ) ;
               } else {
                  sprintf( htitle, "%s, Nb%d, HT%d, MET%d, all QCD bins", mindphi_var, nb_bin, ht_bin, met_bin ) ;
               }
               h_all_qcdbins -> SetTitle( htitle ) ;
            } else {
               h_all_qcdbins -> Add( h_mdp_all[si] ) ;
            }
         } // si

         can2 -> cd(nsamples+1) ;
         double hmax = h_all_qcdbins -> GetMaximum() ;
         h_all_qcdbins -> SetMaximum( 1.2*hmax ) ;
         h_all_qcdbins -> SetMinimum( -0.1*hmax ) ;
         h_all_qcdbins -> Draw() ;


         float rval, rerr ;
         get_ratio( h_all_qcdbins, mdp_cut, rval, rerr ) ;

         printf( "All QCD bins : %5.3f +/- %5.3f\n", rval, rerr ) ;
         char ratio_text[100] ;
         sprintf( ratio_text, "%5.3f +/- %5.3f\n", rval, rerr ) ;

         text.DrawTextNDC( 0.30, 0.80, ratio_text ) ;

         line.SetLineColor(2) ;
         line.SetLineStyle(2) ;
         line.DrawLine( mdp_cut, 0., mdp_cut, h_all_qcdbins->GetMaximum() ) ;

         line.SetLineColor(4) ;
         line.SetLineStyle(3) ;
         line.DrawLine( 0., 0., h_all_qcdbins->GetXaxis() -> GetXmax(), 0. ) ;

         can2 -> Update() ; can2 -> Draw() ;
      }

    //----------

      TH1F* h_qcdbins_ge300 ;
      {
         //for ( int si=5; si<nsamples; si++ ) {
         //   if ( si==5 ) {
         for ( int si=0; si<nsamples; si++ ) {
            if ( si==0 ) {
               char hname[100] ;
               sprintf( hname, "h_qcdbins_ge300_nb%d_ht%d_met%d%s", nb_bin, ht_bin, met_bin, njetbin_str ) ;
               h_qcdbins_ge300 = (TH1F*) h_mdp_all[si]->Clone( hname ) ;
               char htitle[1000] ;
               if ( njet_bin >= 0 ) {
                  sprintf( htitle, "%s, Njet%d, Nb%d, HT%d, MET%d, QCD bins >=300", mindphi_var, njet_bin, nb_bin, ht_bin, met_bin ) ;
               } else {
                  sprintf( htitle, "%s, Nb%d, HT%d, MET%d, QCD bins >=300", mindphi_var, nb_bin, ht_bin, met_bin ) ;
               }
               h_qcdbins_ge300 -> SetTitle( htitle ) ;
            } else {
               h_qcdbins_ge300 -> Add( h_mdp_all[si] ) ;
            }
         } // si

         can2 -> cd(nsamples+1) ;
         double hmax = h_qcdbins_ge300 -> GetMaximum() ;
         h_qcdbins_ge300 -> SetMaximum( 1.2*hmax ) ;
         h_qcdbins_ge300 -> SetMinimum( -0.1*hmax ) ;
      }



    //----------

      TH1F* h_only_good_stats_qcdbins ;
      {
         printf("\n\n") ;
         int min_entries_for_usage(50) ; // Orig 100
         for ( int si=0; si<nsamples; si++ ) {
            if ( si==0 ) {
               char hname[100] ;
               sprintf( hname, "h_only_good_stats_qcdbins_nb%d_ht%d_met%d%s", nb_bin, ht_bin, met_bin, njetbin_str ) ;
               h_only_good_stats_qcdbins = (TH1F*) h_mdp_all[si]->Clone( hname ) ;
               if ( h_mdp_all[si] -> GetEntries() < min_entries_for_usage ) { h_only_good_stats_qcdbins -> Reset() ; }
               //if ( h_mdp_all[si] -> GetEntries() < min_entries_for_usage ) { cout << "FAILED min_entries! Clearing" << endl ; h_only_good_stats_qcdbins -> Clear() ; cout << "  " << h_only_good_stats_qcdbins->GetEntries() << "  "  << endl; h_only_good_stats_qcdbins->Print("all"); }
               //else { cout << "It passed.. " << h_mdp_all[si]->GetEntries() << "  " << min_entries_for_usage << endl;}
               char htitle[1000] ;
               if ( njet_bin >= 0 ) {
                  sprintf( htitle, "%s, Njet%d, Nb%d, HT%d, MET%d, good stats QCD bins", mindphi_var, njet_bin, nb_bin, ht_bin, met_bin ) ;
               } else {
                  sprintf( htitle, "%s, Nb%d, HT%d, MET%d, good stats QCD bins", mindphi_var, nb_bin, ht_bin, met_bin ) ;
               }
               h_only_good_stats_qcdbins -> SetName( hname ) ;
               h_only_good_stats_qcdbins -> SetTitle( htitle ) ;
            } else {
               if ( h_mdp_all[si] -> GetEntries() > min_entries_for_usage ) {
                  printf(" %30s has %.0f entries.  Including it.\n", h_mdp_all[si] -> GetName(), h_mdp_all[si] -> GetEntries() ) ;
                  h_only_good_stats_qcdbins -> Add( h_mdp_all[si] ) ;
               }
            }
         } // si

         can2 -> cd(nsamples+2) ;
         double hmax = h_only_good_stats_qcdbins -> GetMaximum() ;
         h_only_good_stats_qcdbins -> SetMaximum( 1.2*hmax ) ;
         h_only_good_stats_qcdbins -> SetMinimum( -0.1*hmax ) ;
         h_only_good_stats_qcdbins -> Draw() ;

         TH1F* h_overlay = (TH1F*) h_all_qcdbins -> Clone( "h_overlay" ) ;
         h_overlay -> SetLineColor(33) ;
         h_overlay -> Draw( "same" ) ;
         h_only_good_stats_qcdbins -> Draw("same") ;

         float rval, rerr ;
         get_ratio( h_only_good_stats_qcdbins, mdp_cut, rval, rerr ) ;

         printf( "Only QCD bins with >%d entries : %5.3f +/- %5.3f\n", min_entries_for_usage, rval, rerr ) ;
         char ratio_text[100] ;
         sprintf( ratio_text, "Ratio: %5.3f +/- %5.3f\n", rval, rerr ) ;
         text.DrawTextNDC( 0.30, 0.60, ratio_text ) ;

         double zm_rms_val, zm_rms_err ;
         calc_0m_rms( h_only_good_stats_qcdbins, zm_rms_val, zm_rms_err ) ;
         char zmrms_text[100] ;
         sprintf( zmrms_text, "RMS0: %.3f +/- %0.3f\n", zm_rms_val, zm_rms_err ) ;
         text.DrawTextNDC( 0.30, 0.50, zmrms_text ) ;

         line.SetLineColor(2) ;
         line.SetLineStyle(2) ;
         line.DrawLine( mdp_cut, 0., mdp_cut, h_only_good_stats_qcdbins->GetMaximum() ) ;

         line.SetLineColor(4) ;
         line.SetLineStyle(3) ;
         line.DrawLine( 0., 0., h_only_good_stats_qcdbins->GetXaxis() -> GetXmax(), 0. ) ;


         can2 -> Update() ; can2 -> Draw() ;
      }

      char fname[10000] ;
      if ( njet_bin >= 0 ) {
         sprintf( fname, "outputfiles/%s-njet%d-nb%d-ht%d-met%d-cut%4.2f.pdf", mindphi_var, njet_bin, nb_bin, ht_bin, met_bin, mdp_cut ) ;
      } else {
         sprintf( fname, "outputfiles/%s-nb%d-ht%d-met%d-cut%4.2f.pdf", mindphi_var, nb_bin, ht_bin, met_bin, mdp_cut ) ;
      }
      can2 -> SaveAs( fname ) ;





      for ( int si=0; si<nsamples; si++ ) { delete sch[si] ; }

   } // draw_qcd_dphi_plot4

   //================================================================================================

   void add_overflow( TH1F* hp ) {
      if ( hp == 0x0 ) return ;
      int nb = hp -> GetNbinsX() ;
      float new_last_content = hp -> GetBinContent( nb ) + hp -> GetBinContent( nb+1 ) ;
      float new_last_err = sqrt( pow( hp -> GetBinError( nb), 2 ) + pow( hp -> GetBinContent( nb+1), 2 ) ) ;
      hp -> SetBinContent( nb, new_last_content ) ;
      hp -> SetBinError( nb, new_last_err ) ;
   } // add_overflow

  //==================================================================================================

   void get_ratio( TH1F* hp, float cut, float& rval, float& rerr ) {

      rval = 0. ;
      rerr = 0. ;
      if ( hp == 0x0 ) return;
      int cbi(0) ;
      TAxis* xaxis = hp -> GetXaxis() ;
      int nbins = hp->GetNbinsX() ;
      for ( int bi=1; bi<=nbins; bi++ ) {
         if ( xaxis -> GetBinLowEdge(bi) >= cut ) {
            cbi = bi ;
            break ;
         }
      } // bi
      //////////// printf("  Cut of %.2f is bin %d (low edge = %.2f)\n", cut, cbi, xaxis -> GetBinLowEdge(cbi) ) ;
      double ivl, iel, ivh, ieh ;
      ivl = hp -> IntegralAndError( 1, cbi-1, iel ) ;
      ivh = hp -> IntegralAndError( cbi, nbins, ieh ) ;
      if ( ivl > 0 ) rval = ivh / ivl ;
      if ( ivh > 0 && ivl > 0 ) rerr = rval * sqrt( pow(iel/ivl,2) + pow(ieh/ivh,2) ) ;
      return ;

   }

  //==================================================================================================


