
#include "TChain.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "histio.c"


  //--------------------------------------------------------------------------------------------------------

   int nhtbins = 4 ;   // first bin is not used in analysis.
   double htbins[5] = { 0., 500., 800., 1200., 20000. } ;  // first bin is not used in analysis.

   //////int nmhtbins = 6 ;  // first two bins are not used in analysis.
   //////double mhtbins[7] = { 0., 150., 200., 300., 500., 750., 20000. } ;    // first two bins are not used in analysis.
   int nmhtbins ;
   double mhtbins[10] ;

   int nnbbins = 4 ;
   double nbbins[5] = { -0.5, 0.5, 1.5, 2.5, 20. } ;

   int nnjetbins = 6 ; // first bin is not used in analysis.
   double njetbins[7] = { 0., 3.5, 4.5, 5.5, 6.5, 8.5, 20. } ;  // first bin is not used in analysis.

  //--------------------------------------------------------------------------------------------------------






   void fill_bg_hists_tm4( float integrated_lumi_ipb = 10000., bool sum_over_nb=false, bool include_mht0=false ) {

      float amin_lumi = 52.181 ;

      if ( include_mht0 ) {
         nmhtbins = 6 ;
         mhtbins[0] = 0. ;
         mhtbins[1] = 150. ;
         mhtbins[2] = 200. ;
         mhtbins[3] = 300. ;
         mhtbins[4] = 500. ;
         mhtbins[5] = 750. ;
         mhtbins[6] = 20000. ;
      } else {
         nmhtbins = 5 ;
         mhtbins[0] = 0. ;
         mhtbins[1] = 200. ;
         mhtbins[2] = 300. ;
         mhtbins[3] = 500. ;
         mhtbins[4] = 750. ;
         mhtbins[5] = 20000. ;
      }

      if ( sum_over_nb ) {
         nnbbins = 1 ;
         nbbins[1] = 20. ;
      }

      gDirectory -> Delete( "h*" ) ;

      new TH1F( "h_binning_ht_bg", "HT binning, fill_bg_hists_tm4", nhtbins, htbins ) ;
      new TH1F( "h_binning_mht_bg", "MHT binning, fill_bg_hists_tm4", nmhtbins, mhtbins ) ;
      new TH1F( "h_binning_nb_bg", "Nb binning, fill_bg_hists_tm4", nnbbins, nbbins ) ;
      new TH1F( "h_binning_njets_bg", "Njets binning, fill_bg_hists_tm4", nnjetbins, njetbins ) ;

      int n_selections = 4 ;
      char selname[4][100] ;
      sprintf( selname[0], "zl" ) ;
      sprintf( selname[1], "sl" ) ;
      sprintf( selname[2], "ldp" ) ;
      sprintf( selname[3], "slldp" ) ;

      char input_dir[100000] ;
      /////////sprintf( input_dir, "/Users/owen/work/cms/ra2b-2015/tree-maker/amin-aug6-2015/skim/" ) ;
      sprintf( input_dir, "current-tree-dir/" ) ;

      TChain* sch[100] ;
      int n_samples(0) ;
      char sname[100][100] ;
      float k_factor[100] ;


      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_fill_bg_hists" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_fill_bg_hists", "BG MHT vs HT", 700, 500 ) ;
      can -> Clear() ;


      {  int si ; char file_pattern[100000] ; int n_added ; int n_entries ;

         si = n_samples ;
         sprintf( sname[si], "ttbar" ) ;
         ///////////sprintf( file_pattern, "%s/ttbar.root", input_dir ) ;
         sprintf( file_pattern, "%s/TTJets*.root", input_dir ) ;
         sch[si] = new TChain( "PreSelection" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         k_factor[si] = 1.0 ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "wjets" ) ;
         ///////////sprintf( file_pattern, "%s/wjets.root", input_dir ) ;
         sprintf( file_pattern, "%s/WJets*.root", input_dir ) ;
         sch[si] = new TChain( "PreSelection" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 4), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 4 ) return ;
         if ( n_entries <= 0 ) return ;
         /////////////k_factor[si] = 1.23 ;
         k_factor[si] = 1.00 ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "znunu" ) ;
         ////////////sprintf( file_pattern, "%s/znunu.root", input_dir ) ;
         sprintf( file_pattern, "%s/ZJetsToNuNu*.root", input_dir ) ;
         sch[si] = new TChain( "PreSelection" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 4), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 4 ) return ;
         if ( n_entries <= 0 ) return ;
         /////////////k_factor[si] = 1.27 ;
         k_factor[si] = 1.00 ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "sngltop" ) ;
         sch[si] = new TChain( "PreSelection" ) ;
         ////////////sprintf( file_pattern, "%s/singletop.root", input_dir ) ;
         sprintf( file_pattern, "%s/ST_*.root", input_dir ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 5), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 5 ) return ;
         if ( n_entries <= 0 ) return ;
         k_factor[si] = 1.0 ;
         n_samples ++ ;

         n_added = 0 ;
         si = n_samples ;
         sprintf( sname[si], "qcd" ) ;
         sch[si] = new TChain( "PreSelection" ) ;
         ///////////sprintf( file_pattern, "%s/qcd-pt-ge300.root", input_dir ) ;
         sprintf( file_pattern, "%s/QCD*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 10), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 10 ) return ;
         if ( n_entries <= 0 ) return ;
         k_factor[si] = 1.0 ;
         n_samples ++ ;

      } // scoping bracket


      for ( int selind=0; selind<n_selections; selind++ ) {

         char leptons_cut[1000] ;
         char ldp_cut[1000] ;

         if ( strcmp( selname[selind], "zl" ) == 0 ) {
            sprintf( leptons_cut, "(MuonsNum+ElectronsNum) == 0 && (isoElectronTracks+isoMuonTracks+isoPionTracks)==0" ) ;
            sprintf( ldp_cut, "(DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3)" ) ;
         } else if ( strcmp( selname[selind], "sl" ) == 0 ) {
            sprintf( leptons_cut, "(MuonsNum+ElectronsNum) == 1" ) ;
            sprintf( ldp_cut, "(DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3)" ) ;
         } else if ( strcmp( selname[selind], "ldp" ) == 0 ) {
            sprintf( leptons_cut, "(MuonsNum+ElectronsNum) == 0 && (isoElectronTracks+isoMuonTracks+isoPionTracks)==0" ) ;
            sprintf( ldp_cut, "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3)" ) ;
         } else if ( strcmp( selname[selind], "slldp" ) == 0 ) {
            sprintf( leptons_cut, "(MuonsNum+ElectronsNum) == 1" ) ;
            sprintf( ldp_cut, "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3)" ) ;
         } else {
            printf("\n\n Unknown selection name %d: %s\n\n", selind, selname[selind] ) ; return ;
         }

         char sel_cut[1000] ;
         sprintf( sel_cut, "(%s) && (%s)", leptons_cut, ldp_cut ) ;

         printf("  Selection %d : %5s : %s\n", selind, selname[selind], sel_cut ) ;

         for ( int si=0; si<n_samples; si++ ) {

            for ( int nbi=0; nbi<nnbbins; nbi++ ) {

               char nb_cut[1000] ;
               sprintf( nb_cut, "BTags > %.1f && BTags <= %.1f", nbbins[nbi], nbbins[nbi+1] ) ;

               printf("  Nb %d : %s\n", nbi, nb_cut ) ;

               char hname[1000] ;
               sprintf( hname, "h_njvsmhtvsht_nb%d_%s_%s", nbi, selname[selind], sname[si] ) ;

               char all_cuts[10000] ;
               //////////////sprintf( all_cuts, "((JetID>0 && (isoElectronTracks+isoMuonTracks+isoPionTracks)==0) && (%s) && (%s))*(Weight/4000.)*(%.3f)*%.0f", sel_cut, nb_cut, k_factor[si], integrated_lumi_ipb ) ;
               sprintf( all_cuts, "((JetID>0) && (%s) && (%s))*Weight*(%.3f/%.3f)*(%.3f)", sel_cut, nb_cut, integrated_lumi_ipb, amin_lumi, k_factor[si] ) ;

               printf("        %s : sel%d, sample%d, Nb%d, %s\n", hname, selind, si, nbi, all_cuts ) ;

               TH3F* hp = new TH3F( hname, hname, nhtbins, htbins, nmhtbins, mhtbins, nnjetbins, njetbins ) ;
               hp -> Sumw2() ;

               char arg1[1000] ;
               sprintf( arg1, "NJets:MHT:HT>>%s", hname ) ;

               sch[si] -> Draw( arg1, all_cuts, "box" ) ;
               can -> Update() ;


            } // nbi.

         } // si.

      } // selind.



      gSystem -> Exec( "mkdir -p outputfiles" ) ;

      char outfile[1000] ;
      if ( sum_over_nb ) {
         sprintf( outfile, "outputfiles/fill-bg-hists-tm4-nbsum.root" ) ;
      } else {
         sprintf( outfile, "outputfiles/fill-bg-hists-tm4.root" ) ;
      }

      printf("\n\n Saving histograms in %s\n\n", outfile ) ;
      saveHist( outfile, "h*" ) ;





   } // fill_bg_hists4.

   //================================================================================================




