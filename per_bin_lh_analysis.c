
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
#include "TLine.h"
#include "TText.h"
#include "TGraphAsymmErrors.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooMinuit.h"

#include <iostream>
#include <sstream>

  RooWorkspace* lws ;
  TCanvas* cscan ;

  using namespace RooFit;
  using namespace RooStats;

  TString* scan_dir ;

  //using std::cout ;
  //using std::endl ;

   void scan_bg( RooAbsPdf& lh, RooAbsPdf* ignore_pdf, RooDataSet* rds, const char* sb_name,
                 double& val, double& err_low_stat_only, double& err_high_stat_only, double& err_low, double& err_high ) ;
   void fix_pars( RooArgList& plist, float val ) ;
   void free_pars( RooArgList& plist ) ;


  //-----------

   void per_bin_lh_analysis( const char* wsfile = "outputfiles/ws-t1bbbbH.root",
                            float ymax = 2.5,
                            bool check_signif_with_all_bins = false ) {


      cscan = 0x0 ;

      scan_dir = new TString( wsfile ) ;
      scan_dir -> ReplaceAll( ".root", "-scans" ) ;
      char command[1000] ;
      sprintf( command, "mkdir -p %s", scan_dir->Data() ) ;
      gSystem -> Exec( command ) ;


      gStyle->SetOptStat(0) ;

      TFile* wstf = new TFile( wsfile ) ;
      lws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );

      TString root_outfile( wsfile ) ;
      root_outfile.ReplaceAll("ws-","per-bin-lh-analysis-") ;
      TFile* tf_output = new TFile( root_outfile, "RECREATE" ) ;



      const RooArgSet* sbIndexList = lws -> set( "sbIndexList" ) ;
      if ( sbIndexList == 0x0 ) { printf("\n\n *** Workspace missing sbIndexList.\n\n") ; return ; }
      RooLinkedListIter sb_index_iter = sbIndexList -> iterator() ;
      printf("\n\n List of all Search Bin indices:\n") ;
      char sb_name[1000][100] ;
      char sb_label[1000][100] ;
      int total_sb(0) ;
      while ( RooConstVar* sb_index = (RooConstVar*) sb_index_iter.Next() ) {
         printf("  %3.0f : %s\n", sb_index->getVal(), sb_index->GetName() ) ;
         TString name( sb_index->GetName() ) ;
         name.ReplaceAll( "sb_index_","") ;
         sprintf( sb_name[total_sb], "%s", name.Data() ) ;
         sprintf( sb_label[total_sb], "%s %3.0f", name.Data(), sb_index->getVal()+1 ) ;
         total_sb ++ ;
      }
      printf("\n\n") ;


      char pname[100] ;

      RooAbsPdf* likelihood = lws->pdf("likelihood") ;
      if ( likelihood == 0x0 ) { printf("\n\n *** can't find likelihood in workspace.\n\n" ) ; return ; }

      const RooArgList lh_pdf_list = ((RooProdPdf*)likelihood) -> pdfList() ;

      int n_sb(0) ;
      int sbi_list[200] ;
      printf("  SB indices of PDFs in likelihood\n" ) ;
      RooLinkedListIter pdf_iter = lh_pdf_list.iterator() ;
      while ( RooAbsPdf* pdf = (RooAbsPdf*) pdf_iter.Next() ) {
         sprintf( pname, "%s_sb_index", pdf->GetName() ) ;
         const RooConstVar* rv_pdf_sb_index = (const RooConstVar*) lws->obj( pname ) ;
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





      RooDataSet* rds = (RooDataSet*) lws->obj( "observed_rds" ) ;
      cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooRealVar* rv_sig_strength = lws->var("sig_strength") ;
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
      TH1F* h_sig = new TH1F( "h_sig", "Signal at nominal Xsec", 72, 0.5, 72.5 ) ;
      double bg_x[100] ;
      double bg_exl[100] ;
      double bg_exh[100] ;
      double bg_val[100] ;
      double bg_err_low_stat_only[100] ;
      double bg_err_high_stat_only[100] ;
      double bg_err_low[100] ;
      double bg_err_high[100] ;
      for ( int i=0; i<72; i++ ) {
         bg_exl[i] = 0. ;
         bg_exh[i] = 0. ;
         bg_val[i] = 0. ;
         bg_err_low_stat_only[i] = 0. ;
         bg_err_high_stat_only[i] = 0. ;
         bg_err_low[i] = 0. ;
         bg_err_high[i] = 0. ;
      }

      double signif_vals[1000] ;

      for ( int bi=0; bi<n_sb; bi++ ) {

         int sbi = sbi_list[bi] ;

         RooLinkedListIter pdf_iter = lh_pdf_list.iterator() ;
         RooArgSet pdf_subset ;
         RooAbsPdf* ignore_pdf(0x0) ;
         printf("\n\n\n  ========= Finding PDFs in likelihood with SB = %3df\n", sbi ) ;
         while ( RooAbsPdf* pdf = (RooAbsPdf*) pdf_iter.Next() ) {
            sprintf( pname, "%s_sb_index", pdf->GetName() ) ;
            const RooConstVar* rv_pdf_sb_index = (const RooConstVar*) lws->obj( pname ) ;
            if ( rv_pdf_sb_index == 0x0 ) { printf("\n\n *** can't find %s in workspace.\n\n", pname ) ; return ; }
            int this_sbi = TMath::Nint( rv_pdf_sb_index -> getVal() ) ;
            if ( this_sbi == sbi || this_sbi < 0 ) {
               pdf_subset.add( *pdf ) ;
               printf("  %35s : SB index = %3.0f\n", pdf->GetName(), rv_pdf_sb_index->getVal() ) ;
               char zl_poisson_pdf_name[100] ;
               sprintf( zl_poisson_pdf_name, "pdf_zl_%s", sb_name[sbi] ) ;
               if ( strcmp( zl_poisson_pdf_name, pdf->GetName() ) == 0 ) {
                  printf("    This one is the ZL Poisson pdf   %s\n", pdf->GetName() ) ;
                  ignore_pdf = pdf ;
               }
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


         sprintf( pname, "mu_sig0_zl_%s", sb_name[sbi] ) ;
         RooAbsReal* rv_sig0 = (RooAbsReal*) lws -> obj( pname ) ;
         if ( rv_sig0 == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; return ; }




         char bin_label[100] ;

         h_signif_only_active -> SetBinContent( bi+1, signif ) ;
         sprintf( bin_label, "%d", sbi+1 ) ;
         h_signif_only_active -> GetXaxis() -> SetBinLabel( bi+1, bin_label ) ;

         h_signif -> SetBinContent( sbi+1, signif ) ;
         sprintf( bin_label, "%s", sb_label[sbi] ) ;
         h_signif -> GetXaxis() -> SetBinLabel( sbi+1, bin_label ) ;

         printf( " Single-bin-result,  %30s : test_stat = %6.3f,  signif = %5.2f\n\n\n", sb_label[sbi], test_stat_val, signif ) ;

         h_sig -> SetBinContent( sbi+1, rv_sig0->getVal() ) ;
         h_sig -> GetXaxis() -> SetBinLabel( sbi+1, bin_label ) ;

         rv_sig_strength -> setVal( 1. ) ;
         rv_sig_strength -> setConstant( kTRUE ) ;
         double val, err_low_stat_only, err_high_stat_only, err_low, err_high ;
         scan_bg( new_likelihood, ignore_pdf, rds, sb_name[sbi], val, err_low_stat_only, err_high_stat_only, err_low, err_high ) ;

         bg_x[bi] = sbi+1. ;
         bg_val[bi] = val ;
         bg_err_low_stat_only[bi] = err_low_stat_only ;
         bg_err_high_stat_only[bi] = err_high_stat_only ;
         bg_err_low[bi] = err_low ;
         bg_err_high[bi] = err_high ;


      } // bi

      printf("\n\n") ;

      double combined_signif2(0.) ;
      for ( int bi=0; bi<n_sb; bi++ ) {
         combined_signif2 += pow( signif_vals[bi], 2. ) ;
      } // bi

      double combined_signif = sqrt( combined_signif2 ) ;
      printf("\n\n Estimated combined significance: %6.3f\n\n\n", combined_signif ) ;
      if ( check_signif_with_all_bins )  printf("\n\n === Significance with all bins:  %6.3f\n\n", all_bins_signif ) ;


      TGraphAsymmErrors* g_bg_stat_only = new TGraphAsymmErrors( n_sb, bg_x, bg_val, bg_exl, bg_exh, bg_err_low_stat_only, bg_err_high_stat_only ) ;
      g_bg_stat_only -> SetName( "g_bg_stat_only" ) ;
      g_bg_stat_only -> SetMarkerStyle(20) ;
      TGraphAsymmErrors* g_bg           = new TGraphAsymmErrors( n_sb, bg_x, bg_val, bg_exl, bg_exh, bg_err_low, bg_err_high ) ;
      g_bg -> SetName( "g_bg" ) ;
      g_bg -> SetMarkerStyle(20) ;





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
      h_signif -> DrawCopy() ;
      gPad -> SetGridy(1) ;

      h_signif -> Write() ;



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

      TString outfile_signif( wsfile ) ;
      outfile_signif.ReplaceAll( "ws-", "significance-per-bin-" ) ;
      outfile_signif.ReplaceAll( ".root", ".pdf" ) ;
      can -> SaveAs( outfile_signif ) ;



     //---------

      ymax = 2.5* h_sig -> GetMaximum() ;
      h_sig -> SetMaximum( ymax ) ;

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2_signficance_by_sb" ) ;
      if ( can2 == 0x0 ) can2 = new TCanvas( "can2_signficance_by_sb", "Signal by bin", 1300, 650 ) ;
      h_sig -> SetYTitle( "Events" ) ;
      h_sig -> SetLabelOffset( 0.04, "x" ) ;
      h_sig -> SetFillColor( kMagenta ) ;

      h_sig -> DrawCopy() ;
      g_bg_stat_only -> Draw("P") ;
      g_bg -> Draw("P") ;

      h_sig -> Write() ;
      g_bg -> Write() ;
      g_bg_stat_only -> Write() ;

      sprintf( text, "Input file: %s", wsfile ) ;
      ttext -> SetTextAlign(11) ;
      ttext -> SetTextColor( kBlue ) ;
      ttext -> DrawTextNDC( 0.03, 0.95, text ) ;

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
      gPad -> SetGridy(1) ;

      TString outfile_events( wsfile ) ;
      outfile_events.ReplaceAll( "ws-", "events-per-bin-" ) ;
      outfile_events.ReplaceAll( ".root", ".pdf" ) ;
      can2 -> SaveAs( outfile_events ) ;


      tf_output -> Close() ;

   } // per_bin_lh_analysis

  //---------

   void scan_bg( RooAbsPdf& lh, RooAbsPdf* ignore_pdf, RooDataSet* rds, const char* sb_name,
                 double& val, double& err_low_stat_only, double& err_high_stat_only, double& err_low, double& err_high ) {

      val = 0. ;
      err_low_stat_only = 0. ;
      err_high_stat_only = 0. ;
      err_low = 0. ;
      err_high = 0. ;

      printf("\n\n === Scanning BG for %s\n\n", sb_name ) ;


      RooArgList nuisance_pars ;

      const RooArgList lh_pdf_list = ((RooProdPdf&)lh).pdfList() ;

      RooLinkedListIter pdf_iter = lh_pdf_list.iterator() ;
      while ( RooAbsPdf* pdf = (RooAbsPdf*) pdf_iter.Next() ) {
         printf(" pdf: %s", pdf->GetName() ) ;
         TString pdfname( pdf->GetName() ) ;
         if ( pdfname.BeginsWith( "pdf_prim_" ) ) {
            printf(" * nuisance par pdf") ;
            TString npname( pdfname ) ;
            npname.ReplaceAll( "pdf_", "" ) ;
            RooAbsReal* rv_np = (RooAbsReal*) lws -> obj( npname ) ;
            if ( rv_np == 0x0 ) {
               printf("\n\n *** can't find NP with name %s\n\n", npname.Data() ) ;
               return ;
            }
            nuisance_pars.add( *rv_np ) ;
         }
         printf("\n") ;
      }
      nuisance_pars.Print() ;



      char pname[100] ;
      sprintf( pname, "mu_allbg_zl_%s", sb_name ) ;
      RooAbsReal* rv_bg = (RooAbsReal*) lws -> obj( pname ) ;
      if ( rv_bg == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }

      RooFitResult* fitResult = lh.fitTo( *rds, Save(true), Optimize(0), PrintLevel(0), Hesse(true), Strategy(1) ) ;

      float best_val = rv_bg -> getVal() ;
      printf("\n   scan_bg: Best value for %40s : %7.2f\n", rv_bg -> GetName(), rv_bg -> getVal() ) ;
      if ( best_val < 0 ) { printf("\n\n *** wtf?\n\n") ; return ; }

      val = best_val ;

      int n_scan_points(25) ;
      float scan_low  = best_val - 7*sqrt(best_val) ;
      if ( scan_low < 0 ) scan_low = 0. ;
      float scan_low_step = (best_val - scan_low) / n_scan_points ;
      float scan_high = best_val + 7*sqrt(best_val) ;
      float scan_high_step = 7*sqrt(best_val) / n_scan_points ;

      RooRealVar rv_scan_bg_val( "rv_scan_bg_val", "rv_scan_bg_val", rv_bg -> getVal(), scan_low, scan_high ) ;

      RooAbsReal* nll = lh.createNLL( *rds, Verbose(true) ) ;

      float penalty(1) ;
      if ( rv_bg->getVal() > 10 ) {
         //penalty = 10./sqrt(rv_bg->getVal()) ;
         penalty = 5./sqrt(rv_bg->getVal()) ;
      } else if ( rv_bg->getVal() > 5 ) {
         penalty = 50 ;
      } else {
         penalty = 200 ;
      }
      RooRealVar rv_penalty_weight( "rv_penalty_weight", "rv_penalty_weight", penalty ) ;

      RooFormulaVar* new_minuit_var = new RooFormulaVar( "new_minuit_var", "@0+@4*(@1-@2)*(@1-@2)+log(@3)",
                  RooArgList( *nll, rv_scan_bg_val, *rv_bg, *ignore_pdf, rv_penalty_weight ) ) ;

      RooMinuit* rminuit = new RooMinuit( *new_minuit_var ) ;
      rminuit->setPrintLevel(-1) ;
      rminuit->setNoWarn() ;

      RooFormulaVar nll_scan_var( "nll_scan_var", "@0+log(@1)", RooArgList( *nll, *ignore_pdf ) ) ;

      double scan_min_nll(1e9) ;
      double scan_min_nll_stat_only(1e9) ;

      int n_scan_low(0) ;
      double scan_low_x[100] ;
      double scan_low_y[100] ;
      double scan_low_x_stat_only[100] ;
      double scan_low_y_stat_only[100] ;

      for ( int si=0; si<n_scan_points; si++ ) {


         float val = best_val - si*scan_low_step ;
         printf("       Setting %s to %.2f\n", rv_bg -> GetName(), val ) ;
         rv_scan_bg_val.setVal( val ) ;
         rv_scan_bg_val.setConstant( kTRUE ) ;

         free_pars( nuisance_pars ) ;

         printf("    Fit with NPs free.\n" ) ;
         rminuit->migrad() ;
         rminuit->hesse() ;
         RooFitResult* rfr = rminuit->save() ;
         //rfr -> Print("v") ;
         float nll_val = nll_scan_var.getVal() ;
         if ( nll_val < scan_min_nll ) scan_min_nll = nll_val ;
         scan_low_x[n_scan_low] = rv_bg->getVal() ;
         scan_low_y[n_scan_low] = nll_val ;

         fix_pars( nuisance_pars, 0. ) ;

         printf("    Fit with NPs fixed.\n" ) ;
         rminuit->migrad() ;
         rminuit->hesse() ;
         RooFitResult* rfr_stat_only = rminuit->save() ;
         //rfr_stat_only -> Print("v") ;
         float nll_val_stat_only = nll_scan_var.getVal() ;
         if ( nll_val_stat_only < scan_min_nll_stat_only ) scan_min_nll_stat_only = nll_val_stat_only ;
         scan_low_x_stat_only[n_scan_low] = rv_bg->getVal() ;
         scan_low_y_stat_only[n_scan_low] = nll_val_stat_only ;

         free_pars( nuisance_pars ) ;

         printf("    %40s Scan result:  BG = %7.2f,  nll = %.5f,  %.5f\n", sb_name, val, nll_val, nll_val_stat_only ) ;

         n_scan_low++ ;

         if ( (nll_val-scan_min_nll) > 3. ) break ;

         delete rfr ;
         delete rfr_stat_only ;

      } // si

      int n_scan_high(0) ;
      double scan_high_x[100] ;
      double scan_high_y[100] ;
      double scan_high_x_stat_only[100] ;
      double scan_high_y_stat_only[100] ;

      for ( int si=0; si<n_scan_points; si++ ) {

         float val = best_val + si*scan_high_step ;
         printf("       Setting %s to %.2f\n", rv_bg -> GetName(), val ) ;
         rv_scan_bg_val.setVal( val ) ;
         rv_scan_bg_val.setConstant( kTRUE ) ;

         free_pars( nuisance_pars ) ;

         printf("    Fit with NPs free.\n" ) ;
         rminuit->migrad() ;
         rminuit->hesse() ;
         RooFitResult* rfr = rminuit->save() ;
         //rfr -> Print("v") ;
         float nll_val = nll_scan_var.getVal() ;
         if ( nll_val < scan_min_nll ) scan_min_nll = nll_val ;
         scan_high_x[n_scan_high] = rv_bg->getVal() ;
         scan_high_y[n_scan_high] = nll_val ;

         fix_pars( nuisance_pars, 0. ) ;

         printf("    Fit with NPs fixed.\n" ) ;
         rminuit->migrad() ;
         rminuit->hesse() ;
         RooFitResult* rfr_stat_only = rminuit->save() ;
         //rfr_stat_only -> Print("v") ;
         float nll_val_stat_only = nll_scan_var.getVal() ;
         if ( nll_val_stat_only < scan_min_nll_stat_only ) scan_min_nll_stat_only = nll_val_stat_only ;
         scan_high_x_stat_only[n_scan_high] = rv_bg->getVal() ;
         scan_high_y_stat_only[n_scan_high] = nll_val_stat_only ;

         free_pars( nuisance_pars ) ;

         printf("    %40s Scan point result:  BG = %7.2f,  nll = %.5f,  %.5f\n", sb_name, val, nll_val, nll_scan_var.getVal() ) ;

         n_scan_high++ ;

         if ( (nll_val-scan_min_nll) > 3. ) break ;

         delete rfr ;
         delete rfr_stat_only ;


      } // si

      int gr_n(0) ;
      double gr_x[100] ;
      double gr_y[100] ;
      double gr_x_stat_only[100] ;
      double gr_y_stat_only[100] ;
      for ( int i=(n_scan_low-1); i>=0; i-- ) {
         gr_x[gr_n] = scan_low_x[i] ;
         gr_y[gr_n] = scan_low_y[i] - scan_min_nll ;
         gr_x_stat_only[gr_n] = scan_low_x_stat_only[i] ;
         gr_y_stat_only[gr_n] = scan_low_y_stat_only[i] - scan_min_nll_stat_only ;
         printf("     %s  Scan data  %3d : x=%8.2f,  y=%6.3f\n", sb_name, gr_n, gr_x[gr_n], gr_y[gr_n] ) ;
         gr_n ++ ;
      }
      for ( int i=1; i<n_scan_high; i++ ) {
         gr_x[gr_n] = scan_high_x[i] ;
         gr_y[gr_n] = scan_high_y[i] - scan_min_nll ;
         gr_x_stat_only[gr_n] = scan_high_x_stat_only[i] ;
         gr_y_stat_only[gr_n] = scan_high_y_stat_only[i] - scan_min_nll_stat_only ;
         printf("     %s  Scan data  %3d : x=%8.2f,  y=%6.3f\n", sb_name, gr_n, gr_x[gr_n], gr_y[gr_n] ) ;
         gr_n ++ ;
      }


      TGraph* graph = new TGraph( gr_n, gr_x, gr_y ) ;
      TGraph* graph_stat_only = new TGraph( gr_n, gr_x_stat_only, gr_y_stat_only ) ;

      for ( int i=0; i<100; i++ ) {
         double x_up = val + i*(gr_x[gr_n-1]-val)/100. ;
         double y_up = graph -> Eval( x_up, 0, "S" ) ;
         double y_up_stat_only = graph_stat_only -> Eval( x_up, 0, "S" ) ;
         if ( err_high <= 0 && y_up >= 1. ) err_high = x_up - val ;
         if ( err_high_stat_only <= 0 && y_up_stat_only >= 1. ) err_high_stat_only = x_up - val ;
      } // i
      for ( int i=0; i<100; i++ ) {
         double x_down = val - i*(val-gr_x[0])/100. ;
         double y_down = graph -> Eval( x_down, 0, "S" ) ;
         double y_down_stat_only = graph_stat_only -> Eval( x_down, 0, "S" ) ;
         if ( err_low <= 0 && y_down >= 1. ) err_low = val - x_down ;
         if ( err_low_stat_only <= 0 && y_down_stat_only >= 1. ) err_low_stat_only = val - x_down ;
      } // i
      if ( err_low <= 0 ) err_low = val - 0. ;
      if ( err_low_stat_only <= 0 ) err_low_stat_only = val - 0. ;

      printf("\n\n   %40s : Scan result   val = %7.2f + (%7.2f, %7.2f) - (%7.2f, %7.2f)\n\n\n\n",
         sb_name, val, err_high_stat_only, err_high, err_low_stat_only, err_low ) ;

      if ( cscan == 0x0 ) cscan = new TCanvas( "cscan", "Scan", 500, 600 ) ;

      char gname[100] ;
      sprintf( gname, "scan_bg_%s", sb_name ) ;
      graph -> SetName( gname ) ;
      sprintf( gname, "scan_bg_stat_only_%s", sb_name ) ;
      graph_stat_only -> SetName( gname ) ;

      char htitle[1000] ;
      sprintf( htitle, "%s profile likelihood scan: -2ln(L/Lm)", rv_bg->GetName() ) ;
      char hname[1000] ;
      sprintf( hname, "h_scan_%s", rv_bg->GetName() ) ;
      TH1F* hscan = new TH1F( hname, htitle, 10, scan_low, scan_high ) ;
      hscan->SetMinimum(0.) ;
      hscan->SetMaximum(3.0) ;


      hscan->DrawCopy() ;
      cscan->Update() ;
      graph->SetLineColor(4) ;
      graph->SetLineWidth(3) ;
      graph_stat_only->SetLineColor(2) ;
      graph_stat_only->SetLineWidth(3) ;
      graph->Draw("CP") ;
      graph_stat_only->Draw("CP") ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;
      cscan->Update() ;

      char savename[10000] ;
      sprintf( savename, "%s/scan-%s.pdf", scan_dir->Data(), rv_bg->GetName() ) ;
      cscan -> SaveAs( savename ) ;

      graph -> Write() ;
      graph_stat_only -> Write() ;

   } // scan_bg

  //---------

   void fix_pars( RooArgList& plist, float val ) {

      RooLinkedListIter iter = plist.iterator() ;
      while ( RooRealVar* rv = (RooRealVar*) iter.Next() ) {
         //printf(" Fixing %s\n", rv->GetName() ) ;
         rv->setVal( val ) ;
         rv->setConstant( kTRUE ) ;
      }

   } // fix_pars

  //---------

   void free_pars( RooArgList& plist ) {

      RooLinkedListIter iter = plist.iterator() ;
      while ( RooRealVar* rv = (RooRealVar*) iter.Next() ) {
         //printf(" Floating %s\n", rv->GetName() ) ;
         rv->setConstant( kFALSE ) ;
      }

   } // free_pars

  //---------







