
#include "TChain.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
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






   void fill_sig_hists4( float integrated_lumi_ipb = 4000., bool sum_over_nb = false, bool include_mht0=false ) {

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

      new TH1F( "h_binning_ht_sig", "HT binning, fill_sig_hists2", nhtbins, htbins ) ;
      new TH1F( "h_binning_mht_sig", "MHT binning, fill_sig_hists2", nmhtbins, mhtbins ) ;
      new TH1F( "h_binning_nb_sig", "Nb binning, fill_sig_hists2", nnbbins, nbbins ) ;
      new TH1F( "h_binning_njets_sig", "Njets binning, fill_sig_hists2", nnjetbins, njetbins ) ;

      int n_selections = 4 ;
      char selname[4][100] ;
      sprintf( selname[0], "zl" ) ;
      sprintf( selname[1], "sl" ) ;
      sprintf( selname[2], "ldp" ) ;
      sprintf( selname[3], "slldp" ) ;


      char input_dir[100000] ;
      sprintf( input_dir, "current-reducedTree-dir" ) ;
      //////sprintf( input_dir, "/Users/owen/work/cms/ra2b-2015/reducedTree-skim-may05-2015" ) ;
      ///////sprintf( input_dir, "/Users/owen/work/cms/ra2b-2015/reducedTree-skim-june01-2015" ) ;

      TChain* sch[100] ;
      int n_samples(0) ;
      char sname[100][100] ;

      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_fill_sig_hists" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_fill_sig_hists", "Signal MHT vs HT", 700, 500 ) ;
      can -> Clear() ;


      {  int si ; char file_pattern[100000] ; int n_added ; int n_entries ;

         si = n_samples ;
         sprintf( sname[si], "t1bbbbH" ) ;
         sprintf( file_pattern, "%s/*.SMS-T1bbbb_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;


         si = n_samples ;
         sprintf( sname[si], "t1bbbbC" ) ;
         sprintf( file_pattern, "%s/*.SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;


         si = n_samples ;
         sprintf( sname[si], "t1ttttH" ) ;
         sprintf( file_pattern, "%s/*.SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "t1ttttC" ) ;
         sprintf( file_pattern, "%s/*.SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;



         si = n_samples ;
         sprintf( sname[si], "t1qqqqH" ) ;
         sprintf( file_pattern, "%s/*.SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;

         si = n_samples ;
         sprintf( sname[si], "t1qqqqC" ) ;
         sprintf( file_pattern, "%s/*.SMS-T1qqqq_2J_mGl-1000_mLSP-800_Tune4C_13TeV-madgraph-tauola*.root", input_dir ) ;
         sch[si] = new TChain( "reducedTree" ) ;
         n_added = sch[si] -> Add( file_pattern ) ;
         n_entries = sch[si] -> GetEntries() ;
         printf( "  Sample %20s : %3d files (should be 1), %9d unweighted events.\n", sname[si], n_added, n_entries ) ;
         if ( n_added != 1 ) return ;
         if ( n_entries <= 0 ) return ;
         n_samples ++ ;



      } // scoping bracket


      for ( int selind=0; selind<n_selections; selind++ ) {

         char leptons_cut[1000] ;
         char ldp_cut[1000] ;

         if ( strcmp( selname[selind], "zl" ) == 0 ) {
            sprintf( leptons_cut, "nElectrons == 0 && nMuons == 0" ) ;
            //////sprintf( ldp_cut, "minDeltaPhiN_pt30 > 6" ) ;
            ////////sprintf( ldp_cut, "minDeltaPhi30 > 0.4" ) ;
            ////////sprintf( ldp_cut, "(deltaPhi1_20>0.5 && deltaPhi2_20>0.5 && deltaPhi3_20>0.3)" ) ;
            sprintf( ldp_cut, "(deltaPhi1_MHT>0.5 && deltaPhi2_MHT>0.5 && deltaPhi3_MHT>0.3)" ) ;
         } else if ( strcmp( selname[selind], "sl" ) == 0 ) {
            sprintf( leptons_cut, "(nElectrons + nMuons) == 1 && MT_Wlep<100" ) ;
            //////sprintf( ldp_cut, "minDeltaPhiN_pt30 > 6" ) ;
            ////////sprintf( ldp_cut, "minDeltaPhi30 > 0.4" ) ;
            ////////sprintf( ldp_cut, "(deltaPhi1_20>0.5 && deltaPhi2_20>0.5 && deltaPhi3_20>0.3)" ) ;
            sprintf( ldp_cut, "(deltaPhi1_MHT>0.5 && deltaPhi2_MHT>0.5 && deltaPhi3_MHT>0.3)" ) ;
         } else if ( strcmp( selname[selind], "ldp" ) == 0 ) {
            sprintf( leptons_cut, "nElectrons == 0 && nMuons == 0" ) ;
            //////sprintf( ldp_cut, "minDeltaPhiN_pt30 < 6" ) ;
            ////////sprintf( ldp_cut, "minDeltaPhi30 < 0.4" ) ;
            ////////sprintf( ldp_cut, "(!(deltaPhi1_20>0.5 && deltaPhi2_20>0.5 && deltaPhi3_20>0.3))" ) ;
            sprintf( ldp_cut, "(deltaPhi1_MHT<0.5 || deltaPhi2_MHT<0.5 || deltaPhi3_MHT<0.3)" ) ;
         } else if ( strcmp( selname[selind], "slldp" ) == 0 ) {
            sprintf( leptons_cut, "(nElectrons + nMuons) == 1 && MT_Wlep<100" ) ;
            //////sprintf( ldp_cut, "minDeltaPhiN_pt30 < 6" ) ;
            ////////sprintf( ldp_cut, "minDeltaPhi30 < 0.4" ) ;
            ///////sprintf( ldp_cut, "(!(deltaPhi1_20>0.5 && deltaPhi2_20>0.5 && deltaPhi3_20>0.3))" ) ;
            sprintf( ldp_cut, "(deltaPhi1_MHT<0.5 || deltaPhi2_MHT<0.5 || deltaPhi3_MHT<0.3)" ) ;
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



      char outfile[1000] ;
      if ( sum_over_nb ) {
         sprintf( outfile, "outputfiles/fill-sig-hists4-nbsum.root" ) ;
      } else {
         sprintf( outfile, "outputfiles/fill-sig-hists4.root" ) ;
      }

      printf("\n\n Saving histograms in %s\n\n", outfile ) ;
      saveHist( outfile, "h*" ) ;



   } // fill_sig_hists4.

   //================================================================================================




