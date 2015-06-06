
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






   void fill_bg_hists4( float integrated_lumi_ipb = 4000., bool sum_over_nb=false, bool include_mht0=false ) {

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

      new TH1F( "h_binning_ht_bg", "HT binning, fill_bg_hists2", nhtbins, htbins ) ;
      new TH1F( "h_binning_mht_bg", "MHT binning, fill_bg_hists2", nmhtbins, mhtbins ) ;
      new TH1F( "h_binning_nb_bg", "Nb binning, fill_bg_hists2", nnbbins, nbbins ) ;
      new TH1F( "h_binning_njets_bg", "Njets binning, fill_bg_hists2", nnjetbins, njetbins ) ;

      int n_selections = 4 ;
      char selname[4][100] ;
      sprintf( selname[0], "zl" ) ;
      sprintf( selname[1], "sl" ) ;
      sprintf( selname[2], "ldp" ) ;
      sprintf( selname[3], "slldp" ) ;

      char input_dir[100000] ;
      sprintf( input_dir, "current-reducedTree-dir" ) ;
      //////sprintf( input_dir, "/Users/owen/work/cms/ra2b-2015/reducedTree-skim-may05-2015" ) ;
      //////sprintf( input_dir, "/Users/owen/work/cms/ra2b-2015/reducedTree-skim-june01-2015" ) ;

      TChain* sch[100] ;
      int n_samples(0) ;
      char sname[100][100] ;


      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_fill_bg_hists" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_fill_bg_hists", "BG MHT vs HT", 700, 500 ) ;
      can -> Clear() ;


      {  int si ; char file_pattern[100000] ; int n_added ; int n_entries ;

         si = n_samples ;
         sprintf( sname[si], "ttbar" ) ;
         sprintf( file_pattern, "%s/*TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "wjets" ) ;
         sprintf( file_pattern, "%s/*WJetsToLNu_HT-*_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 4), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 4 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "znunu" ) ;
         sprintf( file_pattern, "%s/*ZJetsToNuNu_HT-*_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 4), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 4 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "sngltop" ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         sprintf( file_pattern, "%s/*TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola*.root", input_dir ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 6), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 6 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

         n_added = 0 ;
         si = n_samples ;
         sprintf( sname[si], "qcd" ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-300to470_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-470to600_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-600to800_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-800to1000_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-1000to1400_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-1400to1800_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-1800to2400_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-2400to3200_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         sprintf( file_pattern, "%s/*QCD_Pt-3200_Tune4C_13TeV_pythia8_*.root", input_dir ) ;
         n_added += sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 9), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 9 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

      } // scoping bracket


      for ( int selind=0; selind<n_selections; selind++ ) {

         char leptons_cut[1000] ;
         char ldp_cut[1000] ;

         if ( strcmp( selname[selind], "zl" ) == 0 ) {
            sprintf( leptons_cut, "nElectrons == 0 && nMuons == 0" ) ;
            sprintf( ldp_cut, "minDeltaPhiN_pt30 > 6" ) ;
         } else if ( strcmp( selname[selind], "sl" ) == 0 ) {
            sprintf( leptons_cut, "(nElectrons + nMuons) == 1 && MT_Wlep<100" ) ;
            sprintf( ldp_cut, "minDeltaPhiN_pt30 > 6" ) ;
         } else if ( strcmp( selname[selind], "ldp" ) == 0 ) {
            sprintf( leptons_cut, "nElectrons == 0 && nMuons == 0" ) ;
            sprintf( ldp_cut, "minDeltaPhiN_pt30 < 6" ) ;
         } else if ( strcmp( selname[selind], "slldp" ) == 0 ) {
            sprintf( leptons_cut, "(nElectrons + nMuons) == 1 && MT_Wlep<100" ) ;
            sprintf( ldp_cut, "minDeltaPhiN_pt30 < 6" ) ;
         } else {
            printf("\n\n Unknown selection name %d: %s\n\n", selind, selname[selind] ) ; return ;
         }

         char sel_cut[1000] ;
         sprintf( sel_cut, "(%s) && (%s)", leptons_cut, ldp_cut ) ;

         printf("  Selection %d : %5s : %s\n", selind, selname[selind], sel_cut ) ;

         for ( int si=0; si<n_samples; si++ ) {

            for ( int nbi=0; nbi<nnbbins; nbi++ ) {

               char nb_cut[1000] ;
               sprintf( nb_cut, "nbjets30 > %.1f && nbjets30 <= %.1f", nbbins[nbi], nbbins[nbi+1] ) ;

               printf("  Nb %d : %s\n", nbi, nb_cut ) ;

               char hname[1000] ;
               sprintf( hname, "h_njvsmhtvsht_nb%d_%s_%s", nbi, selname[selind], sname[si] ) ;

               char all_cuts[10000] ;
               sprintf( all_cuts, "((badjetFilter && PBNRcode>0 && nisotk_15_mT==0) && (%s) && (%s))*weight3*%.0f", sel_cut, nb_cut, integrated_lumi_ipb ) ;
               ////////sprintf( all_cuts, "((badjetFilter && PBNRcode>0) && (%s) && (%s))*weight3*%.0f", sel_cut, nb_cut, integrated_lumi_ipb ) ;

               printf("        %s : sel%d, sample%d, Nb%d, %s\n", hname, selind, si, nbi, all_cuts ) ;

               TH3F* hp = new TH3F( hname, hname, nhtbins, htbins, nmhtbins, mhtbins, nnjetbins, njetbins ) ;
               hp -> Sumw2() ;

               char arg1[1000] ;
               sprintf( arg1, "njets30:MHT:HT30>>%s", hname ) ;

               sch[si] -> Draw( arg1, all_cuts, "box" ) ;
               can -> Update() ;


            } // nbi.

         } // si.

      } // selind.



      gSystem -> Exec( "mkdir -p outputfiles" ) ;

      char outfile[1000] ;
      if ( sum_over_nb ) {
         sprintf( outfile, "outputfiles/fill-bg-hists4-nbsum.root" ) ;
      } else {
         sprintf( outfile, "outputfiles/fill-bg-hists4.root" ) ;
      }

      printf("\n\n Saving histograms in %s\n\n", outfile ) ;
      saveHist( outfile, "h*" ) ;





   } // fill_bg_hists4.

   //================================================================================================




