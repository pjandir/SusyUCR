
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

   void fix_pars_to_current_val( const RooAbsCollection& plist ) ;
   void fix_pars( const RooAbsCollection& plist, float val ) ;
   void free_pars( const RooAbsCollection& plist ) ;

   void qcdfit( const char* wsfile = "outputfiles/ws-sb-sameas-fb-qcdfloat-t1bbbbH.root",
                            bool fix_nuisance_pars = true,
                            bool fix_bg_mu_pars = true
                            ) {

      gStyle->SetOptStat(0) ;

      TFile* wstf = new TFile( wsfile ) ;
      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      if ( ws == 0x0 ) { printf("\n\n *** no workspace in %s!\n\n", wsfile ) ; return ; }

      ws -> Print() ;


      if ( fix_nuisance_pars ) {
         const RooAbsCollection* all_nuisance_pars =  ws -> set( "all_nuisance_pars" ) ;
         if ( all_nuisance_pars == 0x0 ) { printf("\n\n *** Workspace missing all_nuisance_pars list.\n\n") ; return ; }
         fix_pars( *all_nuisance_pars, 0. ) ;
      }

      if ( fix_bg_mu_pars ) {
         const RooAbsCollection* all_bg_mu_pars =  ws -> set( "all_bg_mu_pars" ) ;
         if ( all_bg_mu_pars == 0x0 ) { printf("\n\n *** Workspace missing all_bg_mu_pars list.\n\n") ; return ; }
         fix_pars_to_current_val( *all_bg_mu_pars ) ;
      }

  //  {
  //     RooRealVar* rv_qcdpar ;

  //     rv_qcdpar = (RooRealVar*) ws -> obj( "Kqcd_mht2" ) ;
  //     if ( rv_qcdpar == 0x0 ) { printf("\n\n ** missing par.\n\n") ; return ; }
  //     rv_qcdpar -> setConstant( kTRUE ) ;

  //     rv_qcdpar = (RooRealVar*) ws -> obj( "Kqcd_mht3" ) ;
  //     if ( rv_qcdpar == 0x0 ) { printf("\n\n ** missing par.\n\n") ; return ; }
  //     rv_qcdpar -> setConstant( kTRUE ) ;

  //     rv_qcdpar = (RooRealVar*) ws -> obj( "Kqcd_mht4" ) ;
  //     if ( rv_qcdpar == 0x0 ) { printf("\n\n ** missing par.\n\n") ; return ; }
  //     rv_qcdpar -> setConstant( kTRUE ) ;

  //     rv_qcdpar = (RooRealVar*) ws -> obj( "Kqcd_njet2" ) ;
  //     if ( rv_qcdpar == 0x0 ) { printf("\n\n ** missing par.\n\n") ; return ; }
  //     rv_qcdpar -> setConstant( kTRUE ) ;

  //     rv_qcdpar = (RooRealVar*) ws -> obj( "Kqcd_njet3" ) ;
  //     if ( rv_qcdpar == 0x0 ) { printf("\n\n ** missing par.\n\n") ; return ; }
  //     rv_qcdpar -> setConstant( kTRUE ) ;

  //  }



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
      int sbi_list[900] ;
      double pdfprod(1.) ;
      double sumlogpdf(0.) ;
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
         printf("  %35s : SB index = %3.0f,  value = %8.6f, -ln(val) = %15.4f\n", pdf->GetName(), rv_pdf_sb_index->getVal(), pdf->getVal(), -1.*log(pdf->getVal()) ) ;
         pdfprod = pdfprod * (pdf->getVal()) ;
         sumlogpdf += -1.*log(pdf->getVal()) ;
         printf("      pdf prod = %g ,  sum -ln(pdf) = %g \n", pdfprod, sumlogpdf ) ;
      }

      printf("\n\n\n ======== PDF prod = %g ,    sum -ln(pdf) = %g\n\n", pdfprod, sumlogpdf ) ;





      RooDataSet* rds = (RooDataSet*) ws->obj( "observed_rds" ) ;
      cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }



      printf("\n\n Evaluating negative log likelihood.\n") ;
      RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;
      printf("  Nll value = %g\n\n", nll -> getVal() ) ;


      rv_sig_strength -> setVal( 1.0 ) ;
      rv_sig_strength -> setConstant( kTRUE ) ;

      RooFitResult* fitResult = likelihood -> fitTo( *rds, Save(true), Optimize(0), PrintLevel(3), Hesse(true), Strategy(1) ) ;
      double minNllSusyFloat = fitResult->minNll() ;
      double susy_ss_atMinNll = rv_sig_strength -> getVal() ;







   } // qcdfit


  //---------

   void fix_pars_to_current_val( const RooAbsCollection& plist ) {

      RooLinkedListIter iter = plist.iterator() ;
      while ( RooRealVar* rv = (RooRealVar*) iter.Next() ) {
         //printf(" Fixing %s\n", rv->GetName() ) ;
         rv->setConstant( kTRUE ) ;
      }

   } // fix_pars_to_current_val

  //---------

   void fix_pars( const RooAbsCollection& plist, float val ) {

      RooLinkedListIter iter = plist.iterator() ;
      while ( RooRealVar* rv = (RooRealVar*) iter.Next() ) {
         //printf(" Fixing %s\n", rv->GetName() ) ;
         rv->setVal( val ) ;
         rv->setConstant( kTRUE ) ;
      }

   } // fix_pars

  //---------

   void free_pars( const RooAbsCollection& plist ) {

      RooLinkedListIter iter = plist.iterator() ;
      while ( RooRealVar* rv = (RooRealVar*) iter.Next() ) {
         //printf(" Floating %s\n", rv->GetName() ) ;
         rv->setConstant( kFALSE ) ;
      }

   } // free_pars

  //---------







