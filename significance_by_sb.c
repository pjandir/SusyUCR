
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

   void significance_by_sb( const char* wsfile = "outputfiles/ws-t1bbbbH.root" ) {

      gStyle->SetOptStat(0) ;

      TFile* wstf = new TFile( wsfile ) ;
      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );



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
         RooFitResult* fitResult_sp = new_likelihood.fitTo( *rds, Save(true), Optimize(0), Hesse(false), Strategy(1), PrintLevel(-1) ) ;

         double minNll_sp = fitResult_sp->minNll() ;

         double test_stat_val = 2.*( minNll_sp - minNllSusyFloat ) ;
         delete fitResult_sp ;

         double signif(0.) ;
         if ( test_stat_val > 0 ) signif = sqrt( test_stat_val ) ;

         signif_vals[bi] = signif ;

         printf( " Single-bin-result,  SB %3d : test_stat = %6.3f,  signif = %5.2f\n\n\n", sbi, test_stat_val, signif ) ;


      } // bi

      printf("\n\n") ;

      double combined_signif2(0.) ;
      for ( int bi=0; bi<n_sb; bi++ ) {
         combined_signif2 += pow( signif_vals[bi], 2. ) ;
      } // bi

      double combined_signif = sqrt( combined_signif2 ) ;
      printf("\n\n Estimated combined significance: %6.3f\n\n\n", combined_signif ) ;



   } // significance_by_sb

