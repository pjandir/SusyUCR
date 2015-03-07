
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "TAxis.h"
#include "TLine.h"
#include "TText.h"

#include <iostream>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;

   void significance_by_sb( const char* wsfile = "outputfiles/ws-t1bbbbH.root",
                            const char* outfile = "outputfiles/significance-per-bin-t1bbbbH.pdf",
                            float ymax = 2.5,
                            bool check_signif_with_all_bins = false ) {

      gStyle->SetOptStat(0) ;

      TFile* wstf = new TFile( wsfile ) ;
      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );


      const RooArgSet* sbIndexList = ws -> set( "sbIndexList" ) ;
      if ( sbIndexList == 0x0 ) { printf("\n\n *** Workspace missing sbIndexList.\n\n") ; return ; }
      RooLinkedListIter sb_index_iter = sbIndexList -> iterator() ;
      printf("\n\n List of all Search Bin indices:\n") ;
      char sb_name[1000][100] ;
      int total_sb(0) ;
      while ( RooConstVar* sb_index = (RooConstVar*) sb_index_iter.Next() ) {
         printf("  %3.0f : %s\n", sb_index->getVal(), sb_index->GetName() ) ;
         TString name( sb_index->GetName() ) ;
         name.ReplaceAll( "sb_index_","") ;
         sprintf( sb_name[total_sb], "%s %3.0f", name.Data(), sb_index->getVal()+1 ) ;
         total_sb ++ ;
      }
      printf("\n\n") ;


      char pname[100] ;

      RooAbsPdf* likelihood = ws->pdf("likelihood") ;
      if ( likelihood == 0x0 ) { printf("\n\n *** can't find likelihood in workspace.\n\n" ) ; return ; }

      const RooArgList lh_pdf_list = ((RooProdPdf*)likelihood) -> pdfList() ;

      int n_sb(0) ;
      int sbi_list[200] ;
      printf("  SB indices of PDFs in likelihood\n" ) ;
      RooLinkedListIter pdf_iter = lh_pdf_list.iterator() ;
      while ( RooAbsPdf* pdf = (RooAbsPdf*) pdf_iter.Next() ) {
         sprintf( pname, "%s_sb_index", pdf->GetName() ) ;
         const RooConstVar* rv_pdf_sb_index = (const RooConstVar*) ws->obj( pname ) ;
         if ( rv_pdf_sb_index == 0x0 ) { printf("\n\n *** can't find %s in workspace.\n\n", pname ) ; return ; }
         int sbi = TMath::Nint( rv_pdf_sb_index -> getVal() ) ;
         if ( sbi >= 0 ) {
            bool found(false) ;
            for ( int bi=0; bi<n_sb; bi++ ) {
               if ( sbi == sbi_list[bi] ) { found = true ; break ; }
            }
            if ( !found ) {
               sbi_list[n_sb] = sbi ;
               n_sb ++ ;
            }
         }
         printf("  %35s : SB index = %3.0f", pdf->GetName(), rv_pdf_sb_index->getVal() ) ;
         printf("\n") ;
      }





      RooDataSet* rds = (RooDataSet*) ws->obj( "observed_rds" ) ;
      cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }


      double all_bins_signif(0.) ;
      if ( check_signif_with_all_bins ) {
         rv_sig_strength -> setConstant( kFALSE ) ;

         printf("\n  All bins fit with sig_strength floating.\n") ;
         RooFitResult* fitResult = likelihood -> fitTo( *rds, Save(true), Optimize(0), PrintLevel(0), Hesse(true), Strategy(1) ) ;
         double minNllSusyFloat = fitResult->minNll() ;
         double susy_ss_atMinNll = rv_sig_strength -> getVal() ;

         RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
         RooMsgService::instance().getStream(1).removeTopic(Fitting) ;


         rv_sig_strength -> setVal( 0. ) ;
         rv_sig_strength -> setConstant( kTRUE ) ;

         printf("\n  Fit with sig_strength fixed to zero.\n") ;
         RooFitResult* fitResult_sp = likelihood -> fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(-1) ) ;

         double minNll_sp = fitResult_sp->minNll() ;

         double test_stat_val = 2.*( minNll_sp - minNllSusyFloat ) ;
         delete fitResult_sp ;

         if ( test_stat_val > 0 ) all_bins_signif = sqrt( test_stat_val ) ;
         printf("\n\n === Significance with all bins:  %6.3f\n\n", all_bins_signif ) ;
      }










      TH1F* h_signif_only_active = new TH1F( "h_signif_only_active", "Sensitivity per search bin", n_sb, 0.5, n_sb+0.5 ) ;
      TH1F* h_signif = new TH1F( "h_signif", "Sensitivity per search bin", 72, 0.5, 72.5 ) ;

      double signif_vals[1000] ;

      for ( int bi=0; bi<n_sb; bi++ ) {

         int sbi = sbi_list[bi] ;

         RooLinkedListIter pdf_iter = lh_pdf_list.iterator() ;
         RooArgSet pdf_subset ;
         printf("\n\n\n  ========= Finding PDFs in likelihood with SB = %3df\n", sbi ) ;
         while ( RooAbsPdf* pdf = (RooAbsPdf*) pdf_iter.Next() ) {
            sprintf( pname, "%s_sb_index", pdf->GetName() ) ;
            const RooConstVar* rv_pdf_sb_index = (const RooConstVar*) ws->obj( pname ) ;
            if ( rv_pdf_sb_index == 0x0 ) { printf("\n\n *** can't find %s in workspace.\n\n", pname ) ; return ; }
            int this_sbi = TMath::Nint( rv_pdf_sb_index -> getVal() ) ;
            if ( this_sbi == sbi || this_sbi < 0 ) {
               pdf_subset.add( *pdf ) ;
               printf("  %35s : SB index = %3.0f\n", pdf->GetName(), rv_pdf_sb_index->getVal() ) ;
            }
         }

         sprintf( pname, "likelihood_sb_%d", sbi ) ;
         RooProdPdf new_likelihood( pname, pname, pdf_subset ) ;
         new_likelihood.Print() ;



         rv_sig_strength -> setConstant( kFALSE ) ;

         printf("\n  Fit with sig_strength floating.\n") ;
         RooFitResult* fitResult = new_likelihood.fitTo( *rds, Save(true), Optimize(0), PrintLevel(0), Hesse(true), Strategy(1) ) ;
         double minNllSusyFloat = fitResult->minNll() ;
         double susy_ss_atMinNll = rv_sig_strength -> getVal() ;

         RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
         RooMsgService::instance().getStream(1).removeTopic(Fitting) ;



         rv_sig_strength -> setVal( 0. ) ;
         rv_sig_strength -> setConstant( kTRUE ) ;

         printf("\n  Fit with sig_strength fixed to zero.\n") ;
         RooFitResult* fitResult_sp = new_likelihood.fitTo( *rds, Save(true), Optimize(0), Hesse(true), Strategy(1), PrintLevel(0) ) ;

         double minNll_sp = fitResult_sp->minNll() ;

         double test_stat_val = 2.*( minNll_sp - minNllSusyFloat ) ;
         delete fitResult_sp ;

         double signif(0.) ;
         if ( test_stat_val > 0 ) signif = sqrt( test_stat_val ) ;

         signif_vals[bi] = signif ;

         char bin_label[100] ;

         h_signif_only_active -> SetBinContent( bi+1, signif ) ;
         sprintf( bin_label, "%d", sbi+1 ) ;
         h_signif_only_active -> GetXaxis() -> SetBinLabel( bi+1, bin_label ) ;

         h_signif -> SetBinContent( sbi+1, signif ) ;
         sprintf( bin_label, "%s", sb_name[sbi] ) ;
         h_signif -> GetXaxis() -> SetBinLabel( sbi+1, bin_label ) ;

         printf( " Single-bin-result,  %30s : test_stat = %6.3f,  signif = %5.2f\n\n\n", sb_name[sbi], test_stat_val, signif ) ;


      } // bi

      printf("\n\n") ;

      double combined_signif2(0.) ;
      for ( int bi=0; bi<n_sb; bi++ ) {
         combined_signif2 += pow( signif_vals[bi], 2. ) ;
      } // bi

      double combined_signif = sqrt( combined_signif2 ) ;
      printf("\n\n Estimated combined significance: %6.3f\n\n\n", combined_signif ) ;
      if ( check_signif_with_all_bins )  printf("\n\n === Significance with all bins:  %6.3f\n\n", all_bins_signif ) ;




      gStyle -> SetPadBottomMargin( 0.45 ) ;
      gStyle -> SetPadRightMargin( 0.03 ) ;
      gStyle -> SetOptTitle(0) ;

      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_signficance_by_sb" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_signficance_by_sb", "Significance by bin", 1300, 650 ) ;

      if ( ymax > 0 ) h_signif -> SetMaximum( ymax ) ;
      h_signif -> SetYTitle( "Significance" ) ;
      h_signif -> SetLabelOffset( 0.04, "x" ) ;
      h_signif -> SetFillColor(11) ;
      h_signif -> SetLineWidth(2) ;
      h_signif -> Draw() ;
      gPad -> SetGridy(1) ;




      TText* ttext = new TText() ;
      ttext -> SetTextSize( 0.035 ) ;

      char text[1000] ;

      sprintf( text, "Input file: %s", wsfile ) ;
      ttext -> SetTextColor( kBlue ) ;
      ttext -> DrawTextNDC( 0.03, 0.95, text ) ;

      sprintf( text, "Combined significance: %6.3f", combined_signif ) ;
      ttext -> SetTextColor( kRed ) ;
      ttext -> DrawTextNDC( 0.75, 0.95, text ) ;

      TLine* line = new TLine() ;

      line -> SetLineStyle(1) ;
      line -> SetLineWidth(2) ;
      line -> SetLineColor( kBlue ) ;
      ttext -> SetTextColor( kBlue ) ;
      ttext -> SetTextAlign(22) ;
      ttext -> SetTextSize( 0.030 ) ;
      for ( int i=1; i<=12; i++ ) {
         line -> DrawLine( 6*i+0.5, -0.85*ymax, 6*i+0.5, 0.87*ymax ) ;
         sprintf( text, "Nb%d", (i-1)%4 ) ;
         ttext -> DrawText( 6*i-2.5, 0.84*ymax, text ) ;
      }
      line -> SetLineStyle(1) ;
      line -> SetLineColor( kRed ) ;
      line -> DrawLine( 24.5, -ymax, 24.5, ymax ) ;
      line -> DrawLine( 48.5, -ymax, 48.5, ymax ) ;
      ttext -> SetTextColor( kRed ) ;
      ttext -> SetTextSize( 0.038 ) ;
      ttext -> DrawText( 12.5, 0.92*ymax, "Njet1" ) ;
      ttext -> DrawText( 36.5, 0.92*ymax, "Njet2" ) ;
      ttext -> DrawText( 60.5, 0.92*ymax, "Njet3" ) ;




      can -> SaveAs( outfile ) ;


   } // significance_by_sb

