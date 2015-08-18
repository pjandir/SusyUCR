
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
#include "TH2.h"
#include "TCanvas.h"
#include "TMatrixT.h"
#include "TSystem.h"
#include "THStack.h"
#include "TLegend.h"

#include "histio.c"

#include <iostream>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;

   void fix_pars_to_current_val( const RooAbsCollection& plist ) ;
   void fix_pars( const RooAbsCollection& plist, float val ) ;
   void free_pars( const RooAbsCollection& plist ) ;

   void get_Rqcd_val_and_error( int hbi, int mbi, int nji, float& Rqcd_val, float& Rqcd_err, RooFitResult* rfr, RooArgList& ral_kqcd_pars ) ;

   TH1F* make_subset_hist( TH1F* hp_in, const char* bin_name_substring ) ;

   void draw_stack( const char* bin_name_substring = "", int yaxis_option = 1 ) ;
   void draw_nonqcdsub( const char* bin_name_substring = "" ) ;
   void draw_pull( const char* bin_name_substring = "" ) ;

   TH1F* get_hist( const char* hname ) ;

   char output_dir[10000] ;

   TCanvas* can ;

   RooWorkspace* ws_pointer ;

  //---------------

   void run_qcdlhfit1( const char* wsfile = "outputfiles/ws-kqcd-lhfit-perfect-qcd-closure-random-nobs-with-constraints.root",
                       float fixed_sig_strength = 0.,
                       bool make_all_plots = true,
                       bool fix_nuisance_pars = false,
                       bool fix_bg_mu_pars = false
                      ) {

      char output_file[10000] ;

      char command[10000] ;
      sprintf( command, "basename %s", wsfile ) ;
      TString wsfilenopath = gSystem->GetFromPipe( command ) ;
      wsfilenopath.ReplaceAll(".root","") ;

      sprintf( command, "dirname %s", wsfile ) ;
      printf("\n\n Executing %s\n\n", command ) ;
      TString wsdir = gSystem->GetFromPipe( command ) ;
      sprintf( output_dir, "%s/qcdlhfit-results-%s", wsdir.Data(), wsfilenopath.Data() ) ;
      sprintf( command, "mkdir -p %s", output_dir ) ;
      gSystem -> Exec( command ) ;


      gStyle->SetOptStat(0) ;
      gStyle -> SetEndErrorSize(5) ;

      TFile* wstf = new TFile( wsfile ) ;
      ws_pointer = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      if ( ws_pointer == 0x0 ) { printf("\n\n *** no workspace in %s!\n\n", wsfile ) ; return ; }

      ws_pointer -> Print() ;


      if ( fix_nuisance_pars ) {
         const RooAbsCollection* all_nuisance_pars =  ws_pointer -> set( "all_nuisance_pars" ) ;
         if ( all_nuisance_pars == 0x0 ) { printf("\n\n *** Workspace missing all_nuisance_pars list.\n\n") ; return ; }
         fix_pars( *all_nuisance_pars, 0. ) ;
      }

      if ( fix_bg_mu_pars ) {
         const RooAbsCollection* all_bg_mu_pars =  ws_pointer -> set( "all_bg_mu_pars" ) ;
         if ( all_bg_mu_pars == 0x0 ) { printf("\n\n *** Workspace missing all_bg_mu_pars list.\n\n") ; return ; }
         fix_pars_to_current_val( *all_bg_mu_pars ) ;
      }






      char pname[100] ;

      RooAbsPdf* likelihood = ws_pointer->pdf("likelihood") ;
      if ( likelihood == 0x0 ) { printf("\n\n *** can't find likelihood in workspace.\n\n" ) ; return ; }

      const RooArgList lh_pdf_list = ((RooProdPdf*)likelihood) -> pdfList() ;

      double pdfprod(1.) ;
      double sumlogpdf(0.) ;
      printf("  PDFs in likelihood\n" ) ;
      RooLinkedListIter pdf_iter = lh_pdf_list.iterator() ;
      while ( RooAbsPdf* pdf = (RooAbsPdf*) pdf_iter.Next() ) {
         printf("  %35s :  value = %8.6f, -ln(val) = %15.4f\n", pdf->GetName(), pdf->getVal(), -1.*log(pdf->getVal()) ) ;
         pdfprod = pdfprod * (pdf->getVal()) ;
         sumlogpdf += -1.*log(pdf->getVal()) ;
         printf("      pdf prod = %g ,  sum -ln(pdf) = %g \n", pdfprod, sumlogpdf ) ;
      }

      printf("\n\n\n ======== PDF prod = %g ,    sum -ln(pdf) = %g\n\n", pdfprod, sumlogpdf ) ;





      RooDataSet* rds = (RooDataSet*) ws_pointer->obj( "observed_rds" ) ;
      printf( "\n\n\n  ===== RooDataSet ====================\n\n") ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooRealVar* rv_sig_strength = ws_pointer->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }



      printf("\n\n Evaluating negative log likelihood.\n") ;
      RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;
      printf("  Nll value = %g\n\n", nll -> getVal() ) ;


      rv_sig_strength -> setVal( fixed_sig_strength ) ;
      rv_sig_strength -> setConstant( kTRUE ) ;

      RooFitResult* fitResult = likelihood -> fitTo( *rds, Save(true), Optimize(0), PrintLevel(3), Hesse(true), Strategy(1) ) ;
      double minNllSusyFloat = fitResult->minNll() ;
      double susy_ss_atMinNll = rv_sig_strength -> getVal() ;





      TMatrixDSym cov_mat = fitResult -> covarianceMatrix() ;
      for ( int i=0; i<10; i++ ) {
         printf( "  row %2d : ", i ) ;
         for ( int j=0; j<10; j++ ) {
            printf("  %12.8f  ", cov_mat[i][j] ) ;
         } // j
         printf("\n") ;
      } // i


      gStyle -> SetPadRightMargin(0.20) ;
      gStyle -> SetPadBottomMargin(0.20) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetPaintTextFormat("5.3f") ;
      gDirectory -> cd("Rint:/") ;
      TCanvas* can_cor_mat = (TCanvas*) gDirectory -> FindObject( "can_cor_mat" ) ;
      if ( can_cor_mat == 0x0 ) {
         can_cor_mat = new TCanvas( "can_cor_mat", "Correlation matrix", 900, 900 ) ;
      }
      can_cor_mat -> cd() ;
      TH2* h_cor_mat = fitResult -> correlationHist() ;
   // h_cor_mat -> GetXaxis() -> SetRange(1,10) ;
   // h_cor_mat -> GetYaxis() -> SetRange((h_cor_mat -> GetNbinsY()-9), h_cor_mat -> GetNbinsY()) ;
      h_cor_mat -> Draw("colz") ;
   // h_cor_mat -> Draw("text same") ;
      can_cor_mat -> Update() ; can_cor_mat -> Draw() ;





      FILE* ofp_kqcd_table(0x0) ;
      sprintf( output_file, "%s/kqcd-parameter-fit-results.tex", output_dir ) ;
      if ( (ofp_kqcd_table=fopen( output_file, "w" ))==NULL ) {
         printf("\n\n *** Problem opening %s\n\n", output_file ) ;
         return ;
      }

      fprintf( ofp_kqcd_table, "\\begin{tabular}{|l||c|c|}\n" ) ;
      fprintf( ofp_kqcd_table, "\\hline\n" ) ;
      fprintf( ofp_kqcd_table, " Kqcd par &  Value  &  Constraint  \\\\\n" ) ;
      fprintf( ofp_kqcd_table, "\\hline\n" ) ;
      fprintf( ofp_kqcd_table, "\\hline\n" ) ;

      FILE* ofp_kqcd_fit_results(0x0) ;
      sprintf( output_file, "%s/kqcd-parameter-fit-results.txt", output_dir ) ;
      if ( (ofp_kqcd_fit_results=fopen( output_file, "w" ))==NULL ) {
         printf("\n\n *** Problem opening %s\n\n", output_file ) ;
         return ;
      }


      RooArgList float_pars_final = fitResult -> floatParsFinal() ;
      RooArgList kqcd_pars ;
      for ( int pi=0; pi<float_pars_final.getSize(); pi++ ) {
         TString pname( float_pars_final.at(pi)->GetName() ) ;
         TString latex_pname = pname ;
         latex_pname.ReplaceAll("_"," ") ;
         if ( pname.Index( "Kqcd" ) == 0 ) {
            kqcd_pars.add( *(float_pars_final.at(pi)) ) ;
            printf("  %2d : %30s : final val = %.4f, error = %.4f\n",
               pi, float_pars_final.at(pi)->GetName(),
               ((RooAbsReal*)float_pars_final.at(pi))->getVal(),
               ((RooRealVar*)float_pars_final.at(pi))->getError()
               ) ;
            fprintf( ofp_kqcd_table, " %20s  &  $%5.3f \\pm %5.3f$  &  \\\\\n",
               latex_pname.Data(), ((RooAbsReal*)float_pars_final.at(pi))->getVal(), ((RooRealVar*)float_pars_final.at(pi))->getError()
               ) ;
            fprintf( ofp_kqcd_fit_results, " %20s  %6.4f +/- %6.4f\n",
               pname.Data(), ((RooAbsReal*)float_pars_final.at(pi))->getVal(), ((RooRealVar*)float_pars_final.at(pi))->getError() ) ;
         }
         if ( pname.Index( "prim_Kqcd" ) == 0 ) {
            kqcd_pars.add( *(float_pars_final.at(pi)) ) ;
            printf("  %2d : %30s : final val = %.4f, error = %.4f\n",
               pi, float_pars_final.at(pi)->GetName(),
               ((RooAbsReal*)float_pars_final.at(pi))->getVal(),
               ((RooRealVar*)float_pars_final.at(pi))->getError()
               ) ;
            TString pname2 = pname ;
            pname2.ReplaceAll("prim_","") ;
            latex_pname.ReplaceAll("prim","") ;
            RooAbsReal* rar_par = (RooAbsReal*) ws_pointer -> obj( pname2 ) ;
            if ( rar_par == 0x0 ) { printf("\n\n *** can't find obj %s\n\n", pname2.Data() ) ; return ; }
            char pname3[100] ;
            sprintf( pname3, "mean_%s", pname2.Data() ) ;
            RooAbsReal* rar_mean = (RooAbsReal*) ws_pointer -> obj( pname3 ) ;
            if ( rar_mean == 0x0 ) { printf("\n\n *** can't find obj %s\n\n", pname3 ) ; return ; }
            sprintf( pname3, "sigma_%s", pname2.Data() ) ;
            RooAbsReal* rar_sigma = (RooAbsReal*) ws_pointer -> obj( pname3 ) ;
            if ( rar_sigma == 0x0 ) { printf("\n\n *** can't find obj %s\n\n", pname3 ) ; return ; }
            fprintf( ofp_kqcd_table, " %20s  &  $%5.3f$   &  $%5.3f \\pm %5.3f$  \\\\\n",
             latex_pname.Data(), rar_par->getVal(), rar_mean -> getVal(), rar_sigma -> getVal() ) ;
            fprintf( ofp_kqcd_fit_results, " %20s  %6.4f  constrained by %6.4f +/- %6.4f\n",
                pname2.Data(), rar_par ->getVal(), rar_mean -> getVal(), rar_sigma -> getVal() ) ;
         }
      } // pi

      fprintf( ofp_kqcd_table, "\\hline\n" ) ;
      fprintf( ofp_kqcd_table, "\\end{tabular}\n" ) ;
      fclose( ofp_kqcd_table ) ;

      fclose( ofp_kqcd_fit_results ) ;


      //------

      vector<TH1F*> hist_list ;

      gDirectory -> cd("Rint:/") ;
      gDirectory -> Delete( "h_*" ) ;

      TH1F* h_rqcd = new TH1F( "h_rqcd", "Rqcd", 55, 0.5, 55.5 ) ; hist_list.push_back( h_rqcd ) ;

      TH1F* h_nobs_ldp = new TH1F( "h_nobs_ldp", "Nobs, LDP", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nobs_ldp ) ;
      TH1F* h_nobs_zl = new TH1F( "h_nobs_zl", "Nobs, zl", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nobs_zl ) ;

      TH1F* h_model_ldp = new TH1F( "h_model_ldp", "model, LDP", 55, 0.5, 55.5 ) ; hist_list.push_back( h_model_ldp ) ;
      TH1F* h_model_zl = new TH1F( "h_model_zl", "model, zl", 55, 0.5, 55.5 ) ; hist_list.push_back( h_model_zl ) ;

      TH1F* h_nonqcd_ldp = new TH1F( "h_nonqcd_ldp", "nonQCD, LDP", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nonqcd_ldp ) ;
      TH1F* h_nonqcd_zl = new TH1F( "h_nonqcd_zl", "nonQCD, ZL", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nonqcd_zl ) ;

      TH1F* h_nonqcdsub_ldp_statonly = new TH1F( "h_nonqcdsub_ldp_statonly", "Nobs-nonQCD, LDP", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nonqcdsub_ldp_statonly ) ;
      TH1F* h_nonqcdsub_ldp = new TH1F( "h_nonqcdsub_ldp", "Nobs-nonQCD, LDP", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nonqcdsub_ldp ) ;

      TH1F* h_nonqcdsub_zl_statonly = new TH1F( "h_nonqcdsub_zl_statonly", "Nobs-nonQCD, ZL", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nonqcdsub_zl_statonly ) ;
      TH1F* h_nonqcdsub_zl = new TH1F( "h_nonqcdsub_zl", "Nobs-nonQCD, ZL", 55, 0.5, 55.5 ) ; hist_list.push_back( h_nonqcdsub_zl ) ;

      TH1F* h_lostlep_model_zl = new TH1F( "h_lostlep_model_zl", "lostlep ZL, model fit result", 55, 0.5, 55.5 ) ; hist_list.push_back( h_lostlep_model_zl ) ;
      TH1F* h_hadtau_model_zl = new TH1F( "h_hadtau_model_zl", "hadtau ZL, model fit result", 55, 0.5, 55.5 ) ; hist_list.push_back( h_hadtau_model_zl ) ;
      TH1F* h_znunu_model_zl = new TH1F( "h_znunu_model_zl", "znunu ZL, model fit result", 55, 0.5, 55.5 ) ; hist_list.push_back( h_znunu_model_zl ) ;
      TH1F* h_qcd_model_zl = new TH1F( "h_qcd_model_zl", "QCD ZL, model fit result", 55, 0.5, 55.5 ) ; hist_list.push_back( h_qcd_model_zl ) ;
      TH1F* h_qcd_model_ldp = new TH1F( "h_qcd_model_ldp", "QCD LDP, model fit result", 55, 0.5, 55.5 ) ; hist_list.push_back( h_qcd_model_ldp ) ;

      TH1F* h_sig0_zl = new TH1F( "h_sig0_zl", "SUSY at nominal Xsec, ZL", 55, 0.5, 55.5 ) ; hist_list.push_back( h_sig0_zl ) ;

      TH1F* h_pull_ldp = new TH1F( "h_pull_ldp", "(data-fit)/sqrt(data), LDP", 55, 0.5, 55.5 ) ; hist_list.push_back( h_pull_ldp ) ;
      TH1F* h_pull_zl = new TH1F( "h_pull_zl", "(data-fit)/sqrt(data), ZL", 55, 0.5, 55.5 ) ; hist_list.push_back( h_pull_zl ) ;


      //------

      {
         int bi(1) ;
         int ht1bi(1) ;
         for ( int nji=1; nji<=5; nji++ ) {
            for ( int mbi=1; mbi<=4; mbi++ ) {
               for ( int hbi=1; hbi<=3; hbi++ ) {

                  if ( hbi==1 && mbi==4 ) continue ;

                  char bin_name[1000] ;
                  sprintf( bin_name, "FB-Njet%d-Nbsum-MHT%d-HT%d", nji, mbi, hbi ) ;

                  float Rqcd_val, Rqcd_err ;
                  get_Rqcd_val_and_error( hbi, mbi, nji, Rqcd_val, Rqcd_err, fitResult,  kqcd_pars ) ;


                  sprintf( pname, "Nldp_%s", bin_name ) ;
                  RooAbsReal* rar_nldp = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rar_nldp == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
                  float nldp_val = rar_nldp -> getVal() ;

                  sprintf( pname, "Nzl_%s", bin_name ) ;
                  RooAbsReal* rar_nzl = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rar_nzl == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
                  float nzl_val = rar_nzl -> getVal() ;

                  printf(" %2d %30s  :  Nldp = %8.1f   Nzl = %8.1f\n", bi, bin_name, nldp_val, nzl_val ) ;


                  sprintf( pname, "n_ldp_%s", bin_name ) ;
                  RooAbsReal* rar_model_nldp = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rar_model_nldp == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
                  float model_nldp_val = rar_model_nldp -> getVal() ;

                  sprintf( pname, "n_zl_%s", bin_name ) ;
                  RooAbsReal* rar_model_nzl = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rar_model_nzl == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
                  float model_nzl_val = rar_model_nzl -> getVal() ;






                  double lostlep_ldp_val0(-1.) ;
                  sprintf( pname, "mean_mu_lostlep_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_lostlep_ldp_val0 = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_lostlep_ldp_val0 != 0x0 ) lostlep_ldp_val0 = rv_lostlep_ldp_val0 -> getVal() ;

                  double lostlep_ldp_val(-1.) ;
                  sprintf( pname, "mu_lostlep_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_lostlep_ldp_val = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_lostlep_ldp_val != 0x0 ) lostlep_ldp_val = rv_lostlep_ldp_val -> getVal() ;

                  double lostlep_ldp_err(-1.) ;
                  sprintf( pname, "sigma_mu_lostlep_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_lostlep_ldp_err = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_lostlep_ldp_err != 0x0 ) lostlep_ldp_err = rv_lostlep_ldp_err -> getVal() ;

                  printf("  %2d  %30s : lost lep, mean0 = %6.1f +/- %5.1f,  val = %6.1f\n",
                     bi, bin_name, lostlep_ldp_val0, lostlep_ldp_err, lostlep_ldp_val ) ;


                  double hadtau_ldp_val0(-1.) ;
                  sprintf( pname, "mean_mu_hadtau_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_hadtau_ldp_val0 = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_hadtau_ldp_val0 != 0x0 ) hadtau_ldp_val0 = rv_hadtau_ldp_val0 -> getVal() ;

                  double hadtau_ldp_val(-1.) ;
                  sprintf( pname, "mu_hadtau_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_hadtau_ldp_val = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_hadtau_ldp_val != 0x0 ) hadtau_ldp_val = rv_hadtau_ldp_val -> getVal() ;

                  double hadtau_ldp_err(-1.) ;
                  sprintf( pname, "sigma_mu_hadtau_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_hadtau_ldp_err = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_hadtau_ldp_err != 0x0 ) hadtau_ldp_err = rv_hadtau_ldp_err -> getVal() ;

                  printf("  %2d  %30s : had tau,  mean0 = %6.1f +/- %5.1f,  val = %6.1f\n",
                     bi, bin_name, hadtau_ldp_val0, hadtau_ldp_err, hadtau_ldp_val ) ;


                  double znunu_ldp_val0(-1.) ;
                  sprintf( pname, "mean_mu_znunu_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_znunu_ldp_val0 = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_znunu_ldp_val0 != 0x0 ) znunu_ldp_val0 = rv_znunu_ldp_val0 -> getVal() ;

                  double znunu_ldp_val(-1.) ;
                  sprintf( pname, "mu_znunu_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_znunu_ldp_val = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_znunu_ldp_val != 0x0 ) znunu_ldp_val = rv_znunu_ldp_val -> getVal() ;

                  double znunu_ldp_err(-1.) ;
                  sprintf( pname, "sigma_mu_znunu_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_znunu_ldp_err = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_znunu_ldp_err != 0x0 ) znunu_ldp_err = rv_znunu_ldp_err -> getVal() ;

                  printf("  %2d  %30s : znunu,    mean0 = %6.1f +/- %5.1f,  val = %6.1f\n",
                     bi, bin_name, znunu_ldp_val0, znunu_ldp_err, znunu_ldp_val ) ;


                  double sig0_ldp_val(-1.) ;
                  sprintf( pname, "mu_sig0_ldp_%s", bin_name ) ;
                  RooAbsReal* rv_sig0_ldp = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_sig0_ldp != 0x0 ) sig0_ldp_val = rv_sig0_ldp -> getVal() ;






                  double lostlep_zl_val0(-1.) ;
                  sprintf( pname, "mean_mu_lostlep_zl_%s", bin_name ) ;
                  RooAbsReal* rv_lostlep_zl_val0 = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_lostlep_zl_val0 != 0x0 ) lostlep_zl_val0 = rv_lostlep_zl_val0 -> getVal() ;

                  double lostlep_zl_val(-1.) ;
                  sprintf( pname, "mu_lostlep_zl_%s", bin_name ) ;
                  RooAbsReal* rv_lostlep_zl_val = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_lostlep_zl_val != 0x0 ) lostlep_zl_val = rv_lostlep_zl_val -> getVal() ;

                  double lostlep_zl_err(-1.) ;
                  sprintf( pname, "sigma_mu_lostlep_zl_%s", bin_name ) ;
                  RooAbsReal* rv_lostlep_zl_err = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_lostlep_zl_err != 0x0 ) lostlep_zl_err = rv_lostlep_zl_err -> getVal() ;

                  printf("  %2d  %30s : lost lep, mean0 = %6.1f +/- %5.1f,  val = %6.1f\n",
                     bi, bin_name, lostlep_zl_val0, lostlep_zl_err, lostlep_zl_val ) ;


                  double hadtau_zl_val0(-1.) ;
                  sprintf( pname, "mean_mu_hadtau_zl_%s", bin_name ) ;
                  RooAbsReal* rv_hadtau_zl_val0 = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_hadtau_zl_val0 != 0x0 ) hadtau_zl_val0 = rv_hadtau_zl_val0 -> getVal() ;

                  double hadtau_zl_val(-1.) ;
                  sprintf( pname, "mu_hadtau_zl_%s", bin_name ) ;
                  RooAbsReal* rv_hadtau_zl_val = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_hadtau_zl_val != 0x0 ) hadtau_zl_val = rv_hadtau_zl_val -> getVal() ;

                  double hadtau_zl_err(-1.) ;
                  sprintf( pname, "sigma_mu_hadtau_zl_%s", bin_name ) ;
                  RooAbsReal* rv_hadtau_zl_err = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_hadtau_zl_err != 0x0 ) hadtau_zl_err = rv_hadtau_zl_err -> getVal() ;

                  printf("  %2d  %30s : had tau,  mean0 = %6.1f +/- %5.1f,  val = %6.1f\n",
                     bi, bin_name, hadtau_zl_val0, hadtau_zl_err, hadtau_zl_val ) ;


                  double znunu_zl_val0(-1.) ;
                  sprintf( pname, "mean_mu_znunu_zl_%s", bin_name ) ;
                  RooAbsReal* rv_znunu_zl_val0 = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_znunu_zl_val0 != 0x0 ) znunu_zl_val0 = rv_znunu_zl_val0 -> getVal() ;

                  double znunu_zl_val(-1.) ;
                  sprintf( pname, "mu_znunu_zl_%s", bin_name ) ;
                  RooAbsReal* rv_znunu_zl_val = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_znunu_zl_val != 0x0 ) znunu_zl_val = rv_znunu_zl_val -> getVal() ;

                  double znunu_zl_err(-1.) ;
                  sprintf( pname, "sigma_mu_znunu_zl_%s", bin_name ) ;
                  RooAbsReal* rv_znunu_zl_err = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_znunu_zl_err != 0x0 ) znunu_zl_err = rv_znunu_zl_err -> getVal() ;

                  printf("  %2d  %30s : znunu,    mean0 = %6.1f +/- %5.1f,  val = %6.1f\n",
                     bi, bin_name, znunu_zl_val0, znunu_zl_err, znunu_zl_val ) ;


                  double sig0_zl_val(-1.) ;
                  sprintf( pname, "mu_sig0_zl_%s", bin_name ) ;
                  RooAbsReal* rv_sig0_zl = (RooAbsReal*) ws_pointer -> obj( pname ) ;
                  if ( rv_sig0_zl != 0x0 ) sig0_zl_val = rv_sig0_zl -> getVal() ;




                  float nonqcd_ldp_val = lostlep_ldp_val + hadtau_ldp_val + znunu_ldp_val ;
                  float nonqcd_ldp_err = sqrt( lostlep_ldp_err*lostlep_ldp_err + hadtau_ldp_err*hadtau_ldp_err + znunu_ldp_err*znunu_ldp_err ) ;

                  float nonqcd_zl_val = lostlep_zl_val + hadtau_zl_val + znunu_zl_val ;
                  float nonqcd_zl_err = sqrt( lostlep_zl_err*lostlep_zl_err + hadtau_zl_err*hadtau_zl_err + znunu_zl_err*znunu_zl_err ) ;

                  float nonqcdsub_ldp_val = nldp_val - nonqcd_ldp_val ;
                  float nonqcdsub_ldp_err = sqrt( nldp_val + nonqcd_ldp_err*nonqcd_ldp_err ) ;

                  float nonqcdsub_zl_val = nzl_val - nonqcd_zl_val ;
                  float nonqcdsub_zl_err = sqrt( nzl_val + nonqcd_zl_err*nonqcd_zl_err ) ;



                  double qcd_ldp_val(-1.) ;
                  double qcd_ldp_err(-1.) ;
                  sprintf( pname, "mu_qcd_ldp_%s", bin_name ) ;
                  RooRealVar* rv_qcd_ldp = (RooRealVar*) ws_pointer -> obj( pname ) ;
                  if ( rv_qcd_ldp == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
                  qcd_ldp_val = rv_qcd_ldp -> getVal() ;
                  qcd_ldp_err = rv_qcd_ldp -> getError() ;

                  double qcd_zl_val = qcd_ldp_val * Rqcd_val ;
                  double qcd_zl_err(0.) ;
                  if ( qcd_ldp_val > 0. && Rqcd_val > 0. ) {
                     qcd_zl_err = qcd_zl_val * sqrt( pow( qcd_ldp_err/qcd_ldp_val, 2 ) + pow( Rqcd_err/Rqcd_val, 2 ) ) ;
                  }
                  printf(" %2d  %30s : QCD LDP %7.1f +/- %5.1f,   QCD ZL %7.1f +/- %5.1f\n",
                     bi, bin_name, qcd_ldp_val, qcd_ldp_err, qcd_zl_val, qcd_zl_err ) ;


                  double pull_zl_val(0.) ;
                  if ( nzl_val > 0 ) {
                     pull_zl_val = (nzl_val - model_nzl_val) / sqrt(nzl_val) ;
                  }
                  double pull_ldp_val(0.) ;
                  if ( nldp_val > 0 ) {
                     pull_ldp_val = (nldp_val - model_nldp_val) / sqrt(nldp_val) ;
                  }



                  h_rqcd -> SetBinContent( bi, Rqcd_val ) ;
                  h_rqcd -> SetBinError( bi, Rqcd_err ) ;

                  h_nobs_ldp -> SetBinContent( bi, nldp_val ) ;
                  h_nobs_zl  -> SetBinContent( bi, nzl_val ) ;

                  h_model_ldp -> SetBinContent( bi, model_nldp_val ) ;
                  h_model_zl  -> SetBinContent( bi, model_nzl_val ) ;

                  h_nonqcd_ldp -> SetBinContent( bi, nonqcd_ldp_val ) ;
                  h_nonqcd_ldp -> SetBinError( bi, nonqcd_ldp_err ) ;

                  h_nonqcd_zl -> SetBinContent( bi, nonqcd_zl_val ) ;
                  h_nonqcd_zl -> SetBinError( bi, nonqcd_zl_err ) ;

                  h_nonqcdsub_zl -> SetBinContent( bi, nonqcdsub_zl_val ) ;
                  h_nonqcdsub_zl -> SetBinError( bi, nonqcdsub_zl_err ) ;

                  h_nonqcdsub_zl_statonly -> SetBinContent( bi, nonqcdsub_zl_val ) ;
                  h_nonqcdsub_zl_statonly -> SetBinError( bi, sqrt(nzl_val) ) ;

                  h_nonqcdsub_ldp -> SetBinContent( bi, nonqcdsub_ldp_val ) ;
                  h_nonqcdsub_ldp -> SetBinError( bi, nonqcdsub_ldp_err ) ;

                  h_nonqcdsub_ldp_statonly -> SetBinContent( bi, nonqcdsub_ldp_val ) ;
                  h_nonqcdsub_ldp_statonly -> SetBinError( bi, sqrt(nldp_val) ) ;

                  h_lostlep_model_zl -> SetBinContent( bi, lostlep_zl_val ) ;
                  h_lostlep_model_zl -> SetBinError( bi, lostlep_zl_err ) ;

                  h_hadtau_model_zl -> SetBinContent( bi, hadtau_zl_val ) ;
                  h_hadtau_model_zl -> SetBinError( bi, hadtau_zl_err ) ;

                  h_znunu_model_zl -> SetBinContent( bi, znunu_zl_val ) ;
                  h_znunu_model_zl -> SetBinError( bi, znunu_zl_err ) ;

                  h_qcd_model_zl -> SetBinContent( bi, qcd_zl_val ) ;
                  h_qcd_model_zl -> SetBinError( bi, qcd_zl_err ) ;

                  h_qcd_model_ldp -> SetBinContent( bi, qcd_ldp_val ) ;
                  h_qcd_model_ldp -> SetBinError( bi, qcd_ldp_err ) ;

                  h_pull_zl -> SetBinContent( bi, pull_zl_val ) ;
                  h_pull_ldp -> SetBinContent( bi, pull_ldp_val ) ;

                  h_sig0_zl -> SetBinContent( bi, sig0_zl_val ) ;


                  char bin_label[100] ;
                  sprintf( bin_label, "%s %2d", bin_name, bi ) ;
                  for ( unsigned long hi=0; hi<hist_list.size(); hi++ ) { hist_list.at(hi) -> GetXaxis() -> SetBinLabel( bi, bin_label ) ; }

                  bi ++ ;
               }
            } // mbi
         } // nji
      }

      for ( unsigned long hi=0; hi<hist_list.size(); hi++ ) { hist_list.at(hi) -> GetXaxis() -> LabelsOption("v") ; }



      vector<TString> bin_name_substrings ;
      TString bin_name_substring ;
      bin_name_substring = "" ; bin_name_substrings.push_back( bin_name_substring ) ;
      if ( make_all_plots ) {
         bin_name_substring = "-HT1" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "-HT2" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "-HT3" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "MHT1" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "MHT2" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "MHT3" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "MHT4" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "Njet1" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "Njet2" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "Njet3" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "Njet4" ; bin_name_substrings.push_back( bin_name_substring ) ;
         bin_name_substring = "Njet5" ; bin_name_substrings.push_back( bin_name_substring ) ;
      }

      for ( unsigned long hi=0; hi<hist_list.size(); hi++ ) {
         for ( unsigned long si=0; si<bin_name_substrings.size(); si++ ) {
            make_subset_hist( hist_list.at(hi), bin_name_substrings.at(si).Data() ) ;
         } // si
      } // hi






      FILE* ofp_qcd_table1(0x0) ;
      sprintf( output_file, "%s/qcd-event-yield-fit-results1.tex", output_dir ) ;
      if ( (ofp_qcd_table1=fopen( output_file, "w" ))==NULL ) {
         printf("\n\n *** Problem opening %s\n\n", output_file ) ;
         return ;
      }

      FILE* ofp_qcd_table2(0x0) ;
      sprintf( output_file, "%s/qcd-event-yield-fit-results2.tex", output_dir ) ;
      if ( (ofp_qcd_table2=fopen( output_file, "w" ))==NULL ) {
         printf("\n\n *** Problem opening %s\n\n", output_file ) ;
         return ;
      }

      FILE* ofp_qcd_table3(0x0) ;
      sprintf( output_file, "%s/qcd-event-yield-fit-results3.tex", output_dir ) ;
      if ( (ofp_qcd_table3=fopen( output_file, "w" ))==NULL ) {
         printf("\n\n *** Problem opening %s\n\n", output_file ) ;
         return ;
      }

      fprintf( ofp_qcd_table1, "\\begin{tabular}{|l||c|c|c||c|c|c|}\n" ) ;
      fprintf( ofp_qcd_table1, "\\hline\n" ) ;
      fprintf( ofp_qcd_table1, "  &  \\multicolumn{3}{|c||}{  LDP }   &  \\multicolumn{3}{|c|}{ ZL }  \\\\\n" ) ;
      fprintf( ofp_qcd_table1, " Bin &  Nobs  &  Fit  &  Fit QCD   &  Nobs  &  Fit  &  Fit QCD  \\\\\n" ) ;
      fprintf( ofp_qcd_table1, "\\hline\n" ) ;
      fprintf( ofp_qcd_table1, "\\hline\n" ) ;

      fprintf( ofp_qcd_table2, "\\begin{tabular}{|l||c|c|c||c|c|c|}\n" ) ;
      fprintf( ofp_qcd_table2, "\\hline\n" ) ;
      fprintf( ofp_qcd_table2, "  &  \\multicolumn{3}{|c||}{  LDP }   &  \\multicolumn{3}{|c|}{ ZL }  \\\\\n" ) ;
      fprintf( ofp_qcd_table2, " Bin &  Nobs  &  Non-QCD &  Fit QCD   &  Nobs  &  Non-QCD  &  Fit QCD  \\\\\n" ) ;
      fprintf( ofp_qcd_table2, "\\hline\n" ) ;
      fprintf( ofp_qcd_table2, "\\hline\n" ) ;

      fprintf( ofp_qcd_table3, "\\begin{tabular}{|l||c|c|c|c|c|}\n" ) ;
      fprintf( ofp_qcd_table3, "\\hline\n" ) ;
      fprintf( ofp_qcd_table3, "  &  \\multicolumn{2}{|c||}{  LDP } &    &  \\multicolumn{2}{|c|}{ ZL }  \\\\\n" ) ;
      fprintf( ofp_qcd_table3, " Bin &  Nobs-NonQCD &  Fit QCD   &  Rqcd  &  Nobs-NonQCD  &  Fit QCD  \\\\\n" ) ;
      fprintf( ofp_qcd_table3, "\\hline\n" ) ;
      fprintf( ofp_qcd_table3, "\\hline\n" ) ;

      for ( int bi=1; bi<=55; bi++ ) {

         char bin_name[100] ;
         sprintf( bin_name, "%s", h_nobs_ldp -> GetXaxis() -> GetBinLabel( bi ) ) ;

         fprintf( ofp_qcd_table1,  " %20s &  $%7.0f$  &  $%7.1f$  &  $%7.1f \\pm %5.1f$ &  $%7.0f$  &  $%7.1f$  &  $%7.1f \\pm %5.1f$  \\\\\n",
              bin_name,
              h_nobs_ldp -> GetBinContent(bi),
              h_model_ldp -> GetBinContent(bi),
              h_qcd_model_ldp -> GetBinContent(bi),  h_qcd_model_ldp -> GetBinError(bi),
              h_nobs_zl -> GetBinContent(bi),
              h_model_zl -> GetBinContent(bi),
              h_qcd_model_zl -> GetBinContent(bi),  h_qcd_model_zl -> GetBinError(bi) ) ;

         fprintf( ofp_qcd_table2,  " %20s &  $%7.0f$  &  $%7.1f \\pm %5.1f$  &  $%7.1f \\pm %5.1f$ &  $%7.0f$  &  $%7.1f \\pm %5.1f$  &  $%7.1f \\pm %5.1f$  \\\\\n",
              bin_name,
              h_nobs_ldp -> GetBinContent(bi),
              h_nonqcd_ldp -> GetBinContent(bi),  h_nonqcd_ldp -> GetBinError(bi),
              h_qcd_model_ldp -> GetBinContent(bi),  h_qcd_model_ldp -> GetBinError(bi),
              h_nobs_zl -> GetBinContent(bi),
              h_nonqcd_zl -> GetBinContent(bi),  h_nonqcd_zl -> GetBinError(bi),
              h_qcd_model_zl -> GetBinContent(bi),  h_qcd_model_zl -> GetBinError(bi) ) ;

         fprintf( ofp_qcd_table3,  " %20s &  $%7.1f \\pm %5.1f$  &  $%7.1f \\pm %5.1f$ &    $%5.3f \\pm %5.3f$  &  $%7.1f \\pm %5.1f$  &  $%7.1f \\pm %5.1f$  \\\\\n",
              bin_name,
              h_nonqcdsub_ldp -> GetBinContent(bi),  h_nonqcdsub_ldp -> GetBinError(bi),
              h_qcd_model_ldp -> GetBinContent(bi),  h_qcd_model_ldp -> GetBinError(bi),
              h_rqcd -> GetBinContent(bi),          h_rqcd -> GetBinError(bi),
              h_nonqcdsub_zl -> GetBinContent(bi),  h_nonqcdsub_zl -> GetBinError(bi),
              h_qcd_model_zl -> GetBinContent(bi),  h_qcd_model_zl -> GetBinError(bi) ) ;

         if ( bi%11 == 0 && bi<55 ) fprintf( ofp_qcd_table1,  " \\hline\n" ) ;
         if ( bi%11 == 0 && bi<55 ) fprintf( ofp_qcd_table2,  " \\hline\n" ) ;
         if ( bi%11 == 0 && bi<55 ) fprintf( ofp_qcd_table3,  " \\hline\n" ) ;

      } // bi

      fprintf( ofp_qcd_table1, "\\hline\n" ) ;
      fprintf( ofp_qcd_table1, "\\end{tabular}\n" ) ;
      fclose( ofp_qcd_table1 ) ;

      fprintf( ofp_qcd_table2, "\\hline\n" ) ;
      fprintf( ofp_qcd_table2, "\\end{tabular}\n" ) ;
      fclose( ofp_qcd_table2 ) ;

      fprintf( ofp_qcd_table3, "\\hline\n" ) ;
      fprintf( ofp_qcd_table3, "\\end{tabular}\n" ) ;
      fclose( ofp_qcd_table3 ) ;






      gStyle -> SetPadBottomMargin(0.37 ) ;

      h_rqcd -> SetMarkerStyle(20) ;
      h_rqcd -> SetMaximum(0.50) ;
      h_rqcd -> SetMinimum(-0.10) ;

      gDirectory -> cd("Rint:/") ;
      TCanvas* can_rqcd = (TCanvas*) gDirectory -> FindObject( "can_rqcd" ) ;
      if ( can_rqcd == 0x0 ) {
         can_rqcd = new TCanvas( "can_rqcd", "Rqcd", 1000, 700 ) ;
      }
      can_rqcd -> cd() ;

      h_rqcd -> Draw() ;
      gPad -> SetGridy(1) ;

      sprintf( output_file, "%s/fit-rqcd.pdf", output_dir ) ;
      can_rqcd -> SaveAs( output_file ) ;




      gDirectory -> cd("Rint:/") ;
      TCanvas* can_qcd = (TCanvas*) gDirectory -> FindObject( "can_qcd" ) ;
      if ( can_qcd == 0x0 ) {
         can_qcd = new TCanvas( "can_qcd", "QCD, ZL fit result", 1000, 700 ) ;
      }
      can_qcd -> cd() ;

      h_qcd_model_zl -> Draw() ;

      sprintf( output_file, "%s/fit-qcd.pdf", output_dir ) ;
      can_qcd -> SaveAs( output_file ) ;







      gDirectory -> cd("Rint:/") ;
      TCanvas* can_stack =(TCanvas*) gDirectory -> FindObject( "can_stack" ) ;
      if ( can_stack == 0x0 ) {
         can_stack = new TCanvas( "can_stack", "", 1200, 800 ) ;
      }
      can_stack -> cd() ;
      can = can_stack ;
      for ( unsigned long si=0; si<bin_name_substrings.size(); si++ ) {
         for ( int yi=1; yi<=4; yi++ ) {
            draw_stack( bin_name_substrings.at(si).Data(), yi ) ;
         }
      } // si


      gDirectory -> cd("Rint:/") ;
      TCanvas* can_nonqcdsub =(TCanvas*) gDirectory -> FindObject( "can_nonqcdsub" ) ;
      if ( can_nonqcdsub == 0x0 ) {
         can_nonqcdsub = new TCanvas( "can_nonqcdsub", "", 1200, 800 ) ;
      }
      can_nonqcdsub -> cd() ;
      can = can_nonqcdsub ;
      for ( unsigned long si=0; si<bin_name_substrings.size(); si++ ) {
         draw_nonqcdsub( bin_name_substrings.at(si).Data() ) ;
      } // si



      gDirectory -> cd("Rint:/") ;
      TCanvas* can_pull =(TCanvas*) gDirectory -> FindObject( "can_pull" ) ;
      if ( can_pull == 0x0 ) {
         can_pull = new TCanvas( "can_pull", "", 1200, 800 ) ;
      }
      can_pull -> cd() ;
      can = can_pull ;
      for ( unsigned long si=0; si<bin_name_substrings.size(); si++ ) {
         draw_pull( bin_name_substrings.at(si).Data() ) ;
      } // si



      if ( make_all_plots ) {
         sprintf( command, "cp qcdlhfit-plots.tex %s", output_dir ) ;
         gSystem -> Exec( command ) ;
         sprintf( command, "cd %s ; pdflatex qcdlhfit-plots.tex ; cd -", output_dir ) ;
         gSystem -> Exec( command ) ;
      } else {
         sprintf( command, "cp qcdlhfit-plots-brief.tex %s", output_dir ) ;
         gSystem -> Exec( command ) ;
         sprintf( command, "cd %s ; pdflatex qcdlhfit-plots-brief.tex ; cd -", output_dir ) ;
         gSystem -> Exec( command ) ;
      }

      char hist_file[10000] ;
      sprintf( hist_file, "%s/histograms.root", output_dir ) ;
      printf("\n\n Saving histograms in %s\n\n", hist_file ) ;
      saveHist( hist_file, "h*" ) ;


   } // run_qcdlhfit1

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //---------

   TH1F* get_hist( const char* hname ) {

      TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** Missing histogram : %s \n\n", hname ) ;
         gSystem -> Exit(-1) ;
      }
      return hp ;

   } // get_hist

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

   void get_Rqcd_val_and_error( int hbi, int mbi, int nji, float& Rqcd_val, float& Rqcd_err, RooFitResult* rfr, RooArgList& ral_kqcd_pars ) {

      Rqcd_val = -1 ;
      Rqcd_err = -1 ;

      if ( rfr == 0x0 ) return ;

      char kqcd_ht_name[100] ;
      sprintf( kqcd_ht_name, "Kqcd_ht%d", hbi ) ;

      char kqcd_mht_name[100] ;
      sprintf( kqcd_mht_name, "Kqcd_mht%d", mbi ) ;

      char kqcd_njet_name[100] ;
      sprintf( kqcd_njet_name, "Kqcd_njet%d", nji ) ;

      int npars(0) ;

      RooRealVar* rrv_kqcd_ht = (RooRealVar*) ral_kqcd_pars.find( kqcd_ht_name ) ;
      if ( rrv_kqcd_ht == 0x0 ) return ;
      double kqcd_ht_val = rrv_kqcd_ht -> getVal() ;
      double kqcd_ht_err = rrv_kqcd_ht -> getError() ;
      npars++ ;

      double kqcd_mht_val = 1. ;
      double kqcd_mht_err = 0. ;
      RooRealVar* rrv_kqcd_mht(0x0) ;
      if ( mbi > 1 ) {
         rrv_kqcd_mht = (RooRealVar*) ral_kqcd_pars.find( kqcd_mht_name ) ;
         if ( rrv_kqcd_mht == 0x0 ) {
            char pname[100] ;
            sprintf( pname, "prim_%s", kqcd_mht_name ) ;
            rrv_kqcd_mht = (RooRealVar*) ral_kqcd_pars.find( pname ) ;
            if ( rrv_kqcd_mht == 0x0 ) { printf("\n\n *** Missing %s and %s parameters. \n\n", kqcd_mht_name, pname ) ; return ; }
            RooAbsReal* rar_kqcd_mht = (RooAbsReal*) ws_pointer -> obj( kqcd_mht_name ) ;
            if ( rar_kqcd_mht == 0x0 ) { printf("\n\n *** %s constrained but can't find par in ws.\n\n", kqcd_mht_name ) ; return ; }
            sprintf( pname, "sigma_%s", kqcd_mht_name ) ;
            RooAbsReal* rar_kqcd_mht_sigma = (RooAbsReal*) ws_pointer -> obj( pname ) ;
            if ( rar_kqcd_mht_sigma == 0x0 ) { printf("\n\n *** %s constrained but can't find sigma in ws.\n\n", kqcd_mht_name ) ; return ; }
            kqcd_mht_val = rar_kqcd_mht -> getVal() ;
            kqcd_mht_err = rar_kqcd_mht_sigma -> getVal() ;
         } else {
            kqcd_mht_val = rrv_kqcd_mht -> getVal() ;
            kqcd_mht_err = rrv_kqcd_mht -> getError() ;
         }
         npars++ ;
      }

      double kqcd_njet_val = 1. ;
      double kqcd_njet_err = 0. ;
      RooRealVar* rrv_kqcd_njet(0x0) ;
      if ( nji > 1 ) {
         rrv_kqcd_njet = (RooRealVar*) ral_kqcd_pars.find( kqcd_njet_name ) ;
         if ( rrv_kqcd_njet == 0x0 ) {
            char pname[100] ;
            sprintf( pname, "prim_%s", kqcd_njet_name ) ;
            rrv_kqcd_njet = (RooRealVar*) ral_kqcd_pars.find( pname ) ;
            if ( rrv_kqcd_njet == 0x0 ) { printf("\n\n *** Missing %s and %s parameters. \n\n", kqcd_njet_name, pname ) ; return ; }
            RooAbsReal* rar_kqcd_njet = (RooAbsReal*) ws_pointer -> obj( kqcd_njet_name ) ;
            if ( rar_kqcd_njet == 0x0 ) { printf("\n\n *** %s constrained but can't find par in ws.\n\n", kqcd_njet_name ) ; return ; }
            sprintf( pname, "sigma_%s", kqcd_njet_name ) ;
            RooAbsReal* rar_kqcd_njet_sigma = (RooAbsReal*) ws_pointer -> obj( pname ) ;
            if ( rar_kqcd_njet_sigma == 0x0 ) { printf("\n\n *** %s constrained but can't find sigma in ws.\n\n", kqcd_njet_name ) ; return ; }
            kqcd_njet_val = rar_kqcd_njet -> getVal() ;
            kqcd_njet_err = rar_kqcd_njet_sigma -> getVal() ;
         } else {
            kqcd_njet_val = rrv_kqcd_njet -> getVal() ;
            kqcd_njet_err = rrv_kqcd_njet -> getError() ;
         }
         npars++ ;
      }

      Rqcd_val = kqcd_ht_val * kqcd_mht_val * kqcd_njet_val ;

      TMatrixT<double> pd_col_vec( npars, 1 ) ;
      TMatrixT<double> pd_row_vec( 1, npars ) ;
      TMatrixT<double> cov_mat( npars, npars ) ;

    //--- Terms for HT par.
      pd_col_vec(0,0) = Rqcd_val / kqcd_ht_val ;
      pd_row_vec(0,0) = Rqcd_val / kqcd_ht_val ;
      cov_mat(0,0) = kqcd_ht_err * kqcd_ht_err ;

      if ( mbi > 1 && nji == 1 ) {
    //--- Terms for MHT par.
         pd_col_vec(1,0) = Rqcd_val / kqcd_mht_val ;
         pd_row_vec(0,1) = Rqcd_val / kqcd_mht_val ;
         cov_mat(0,1) = kqcd_ht_err * kqcd_mht_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_mht ) ;
         cov_mat(1,0) = kqcd_ht_err * kqcd_mht_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_mht ) ;
         ///cov_mat(0,1) = 0. ;
         ///cov_mat(1,0) = 0. ;
         cov_mat(1,1) = kqcd_mht_err * kqcd_mht_err ;
      } else if ( nji > 1 && mbi == 1 ) {
    //--- Terms for Njet par.
         pd_col_vec(1,0) = Rqcd_val / kqcd_njet_val ;
         pd_row_vec(0,1) = Rqcd_val / kqcd_njet_val ;
         cov_mat(0,1) = kqcd_ht_err * kqcd_njet_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_njet ) ;
         cov_mat(1,0) = kqcd_ht_err * kqcd_njet_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_njet ) ;
         ///cov_mat(0,1) = 0. ;
         ///cov_mat(1,0) = 0. ;
         cov_mat(1,1) = kqcd_njet_err * kqcd_njet_err ;
      } else if ( nji > 1 && mbi > 1 ) {
    //--- Terms for MHT par.
         pd_col_vec(1,0) = Rqcd_val / kqcd_mht_val ;
         pd_row_vec(0,1) = Rqcd_val / kqcd_mht_val ;
         cov_mat(0,1) = kqcd_ht_err * kqcd_mht_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_mht ) ;
         cov_mat(1,0) = kqcd_ht_err * kqcd_mht_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_mht ) ;
         ///cov_mat(0,1) = 0. ;
         ///cov_mat(1,0) = 0. ;
         cov_mat(1,1) = kqcd_mht_err * kqcd_mht_err ;
    //--- Terms for Njet par.
         pd_col_vec(2,0) = Rqcd_val / kqcd_njet_val ;
         pd_row_vec(0,2) = Rqcd_val / kqcd_njet_val ;
         cov_mat(0,2) = kqcd_ht_err * kqcd_njet_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_njet ) ;
         cov_mat(2,0) = kqcd_ht_err * kqcd_njet_err * rfr->correlation( *rrv_kqcd_ht, *rrv_kqcd_njet ) ;
         ////cov_mat(0,2) = 0. ;
         ////cov_mat(2,0) = 0. ;
         cov_mat(2,2) = kqcd_njet_err * kqcd_njet_err ;

         cov_mat(1,2) = kqcd_mht_err * kqcd_njet_err * rfr->correlation( *rrv_kqcd_mht, *rrv_kqcd_njet ) ;
         cov_mat(2,1) = kqcd_mht_err * kqcd_njet_err * rfr->correlation( *rrv_kqcd_mht, *rrv_kqcd_njet ) ;
         ////cov_mat(1,2) = 0. ;
         ////cov_mat(2,1) = 0. ;
      }

      //cov_mat.Print() ;
      //pd_col_vec.Print() ;
      //pd_row_vec.Print() ;

      TMatrixT<double> cov_times_pd_col( npars, 1 ) ;
      cov_times_pd_col.Mult( cov_mat, pd_col_vec ) ;
      TMatrixT<double> pd_row_times_prod( 1, 1 ) ;
      pd_row_times_prod.Mult( pd_row_vec, cov_times_pd_col ) ;
      Rqcd_err = sqrt( pd_row_times_prod(0,0) ) ;

      printf("   nji=%d, mbi=%d, hbi=%d :  Rqcd_val = kqcd_njet_val * kqcd_mht_val * kqcd_ht_val = %.4f * %.4f * %.4f = %.4f +/- %.4f\n",
         nji, mbi, hbi, kqcd_njet_val, kqcd_mht_val, kqcd_ht_val, Rqcd_val, Rqcd_err ) ;


   } // get_Rqcd_val_and_error

  //---------

   TH1F* make_subset_hist( TH1F* hp_in, const char* bin_name_substring ) {

      if ( hp_in == 0x0 ) return 0x0 ;

      char newname[1000] ;
      TString clean_substring( bin_name_substring ) ;
      clean_substring.ReplaceAll("-","") ;
      clean_substring.ReplaceAll("_","") ;
      sprintf( newname, "%s_%s", hp_in -> GetName(), clean_substring.Data() ) ;
      char newtitle[1000] ;
      sprintf( newtitle, "%s, %s bins", hp_in -> GetTitle(), clean_substring.Data() ) ;


      TH1F* hp_return = (TH1F*) hp_in -> Clone( newname ) ;
      hp_return -> Reset() ;
      hp_return -> SetTitle( newtitle ) ;

      int n_keep_bins(0) ;
      int keep_bins[500] ;

      for ( int bi=1; bi<=hp_in->GetNbinsX(); bi++ ) {
         TString bname( hp_in -> GetXaxis() -> GetBinLabel( bi ) ) ;
         if ( bname.Index( bin_name_substring ) >= 0 ) {
            keep_bins[n_keep_bins] = bi ;
            ////printf("  Keeping bin %2d : %s\n", bi, bname.Data() ) ;
            n_keep_bins ++ ;
         }
      } // bi

      hp_return -> SetBins( n_keep_bins, 0.5, n_keep_bins+0.5 ) ;

      for ( int bi=1; bi<=n_keep_bins; bi++ ) {
         hp_return -> SetBinContent( bi, hp_in -> GetBinContent( keep_bins[bi-1] ) ) ;
         hp_return -> SetBinError( bi, hp_in -> GetBinError( keep_bins[bi-1] ) ) ;
         hp_return -> GetXaxis() -> SetBinLabel( bi, hp_in -> GetXaxis() -> GetBinLabel( keep_bins[bi-1] ) ) ;
      } // bi

      return hp_return ;

   } // make_subset_hist





  //---------

   void draw_stack( const char* bin_name_substring, int yaxis_option ) {

      char hname[1000] ;

      TString clean_substring( bin_name_substring ) ;
      clean_substring.ReplaceAll("-","") ;
      clean_substring.ReplaceAll("_","") ;

      char histend[100] ;
      if ( strlen( bin_name_substring ) > 0 ) {
         sprintf( histend, "_%s", clean_substring.Data() ) ;
      } else {
         sprintf( histend, "" ) ;
      }

      sprintf( hname, "h_nobs_zl%s", histend ) ;
      TH1F* h_nobs = get_hist( hname ) ;

      sprintf( hname, "h_lostlep_model_zl%s", histend ) ;
      TH1F* h_lostlep = get_hist( hname ) ;
      sprintf( hname, "h_hadtau_model_zl%s", histend ) ;
      TH1F* h_hadtau = get_hist( hname ) ;
      sprintf( hname, "h_znunu_model_zl%s", histend ) ;
      TH1F* h_znunu = get_hist( hname ) ;
      sprintf( hname, "h_qcd_model_zl%s", histend ) ;
      TH1F* h_qcd = get_hist( hname ) ;

      sprintf( hname, "h_model_zl%s", histend ) ;
      TH1F* h_all = get_hist( hname ) ;

      sprintf( hname, "h_sig0_zl%s", histend ) ;
      TH1F* h_sig0 = get_hist( hname ) ;

      h_sig0 -> SetLineWidth(3) ;
      h_sig0 -> SetLineColor(kMagenta-3) ;

      h_qcd -> SetFillColor( kRed-9 ) ;
      h_znunu -> SetFillColor( kGreen-7 ) ;
      h_lostlep -> SetFillColor( kBlue-10 ) ;
      h_hadtau -> SetFillColor( kCyan-10 ) ;

      h_nobs -> SetMarkerStyle(20) ;


      sprintf( hname, "h_stack%s", histend ) ;
      //gDirectory -> Delete( hname ) ;
      THStack* hs = new THStack( hname, hname ) ;
      hs -> Add( h_qcd ) ;
      hs -> Add( h_znunu ) ;
      hs -> Add( h_hadtau ) ;
      hs -> Add( h_lostlep ) ;

      TLegend* legend = new TLegend( 0.82, 0.35, 0.98, 0.90 ) ;
      legend -> AddEntry( h_lostlep, "Lost Lep" ) ;
      legend -> AddEntry( h_hadtau, "Had Tau" ) ;
      legend -> AddEntry( h_znunu, "Znunu" ) ;
      legend -> AddEntry( h_qcd, "QCD" ) ;
      legend -> AddEntry( h_nobs, "data (fake)" ) ;

      TLegend* legend_with_sig0 = new TLegend( 0.82, 0.35, 0.98, 0.90 ) ;
      legend_with_sig0 -> AddEntry( h_lostlep, "Lost Lep" ) ;
      legend_with_sig0 -> AddEntry( h_hadtau, "Had Tau" ) ;
      legend_with_sig0 -> AddEntry( h_znunu, "Znunu" ) ;
      legend_with_sig0 -> AddEntry( h_qcd, "QCD" ) ;
      legend_with_sig0 -> AddEntry( h_nobs, "data (fake)" ) ;
      legend_with_sig0 -> AddEntry( h_sig0, "SUSY" ) ;

      h_nobs -> SetMinimum(0.1) ;

      if ( yaxis_option == 3 ) { h_nobs -> SetMaximum(200) ; }
      if ( yaxis_option == 4 ) { h_nobs -> SetMaximum(50) ; }

      h_nobs -> Draw("e") ;
      hs -> Draw( "hist same" ) ;
      hs -> Draw( "axis same" ) ;
      hs -> Draw( "axig same" ) ;
      hs -> Draw( "same" ) ;
      h_nobs -> Draw("e same") ;
      legend -> Draw() ;

      if ( yaxis_option == 2 ) { gPad -> SetLogy(1) ; } else { gPad -> SetLogy(0) ; }

      char file_end[100] ;
      if ( strlen( bin_name_substring ) > 0 ) {
         sprintf( file_end, "%s", bin_name_substring ) ;
      } else {
         sprintf( file_end, "allbins" ) ;
      }

      char yoptstring[100] ;
      if ( yaxis_option == 1 ) sprintf( yoptstring, "liny" ) ;
      if ( yaxis_option == 2 ) sprintf( yoptstring, "logy" ) ;
      if ( yaxis_option == 3 ) sprintf( yoptstring, "zoom1" ) ;
      if ( yaxis_option == 4 ) sprintf( yoptstring, "zoom2" ) ;

      char output_file[10000] ;
      sprintf( output_file, "%s/stack-%s-%s.pdf", output_dir, file_end, yoptstring ) ;
      can -> SaveAs( output_file ) ;

      h_sig0 -> Draw("hist same") ;
      h_nobs -> Draw("e same") ;
      legend_with_sig0 -> Draw() ;

      sprintf( output_file, "%s/stack-%s-with-sig0-%s.pdf", output_dir, file_end, yoptstring ) ;
      can -> SaveAs( output_file ) ;

   } // draw_stack

  //---------

   void draw_nonqcdsub( const char* bin_name_substring ) {

      char hname[1000] ;

      TString clean_substring( bin_name_substring ) ;
      clean_substring.ReplaceAll("-","") ;
      clean_substring.ReplaceAll("_","") ;

      char histend[100] ;
      if ( strlen( bin_name_substring ) > 0 ) {
         sprintf( histend, "_%s", clean_substring.Data() ) ;
      } else {
         sprintf( histend, "" ) ;
      }

      sprintf( hname, "h_nonqcdsub_zl_statonly%s", histend ) ;
      TH1F* h_nonqcdsub_statonly = get_hist( hname ) ;

      sprintf( hname, "h_nonqcdsub_zl%s", histend ) ;
      TH1F* h_nonqcdsub = get_hist( hname ) ;

      sprintf( hname, "h_qcd_model_zl%s", histend ) ;
      TH1F* h_qcd = get_hist( hname ) ;

      sprintf( hname, "h_sig0_zl%s", histend ) ;
      TH1F* h_sig0 = get_hist( hname ) ;

      h_sig0 -> SetLineWidth(3) ;
      h_sig0 -> SetLineColor(kMagenta-3) ;

      h_nonqcdsub_statonly -> SetMarkerStyle(20) ;
      h_nonqcdsub -> SetMarkerStyle(20) ;

      char file_end[100] ;
      if ( strlen( bin_name_substring ) > 0 ) {
         sprintf( file_end, "%s", bin_name_substring ) ;
      } else {
         sprintf( file_end, "allbins" ) ;
      }

      TLegend* legend = new TLegend( 0.82, 0.35, 0.98, 0.90 ) ;
      legend -> AddEntry( h_qcd, "Fit QCD") ;
      legend -> AddEntry( h_nonqcdsub, "Data - NonQCD") ;

      TLegend* legend_with_sig0 = new TLegend( 0.82, 0.35, 0.98, 0.90 ) ;
      legend_with_sig0 -> AddEntry( h_qcd, "Fit QCD") ;
      legend_with_sig0 -> AddEntry( h_nonqcdsub, "Data - NonQCD") ;
      legend_with_sig0 -> AddEntry( h_sig0, "SUSY") ;


      h_nonqcdsub_statonly -> Draw("e1") ;
      h_nonqcdsub -> Draw("e1 same") ;
      h_qcd -> Draw( "hist same" ) ;
      h_qcd -> Draw( "axis same" ) ;
      h_qcd -> Draw( "axig same" ) ;
      h_nonqcdsub_statonly -> Draw("e1 same") ;
      h_nonqcdsub -> Draw("e1 same") ;
      legend -> Draw() ;

      char output_file[10000] ;
      sprintf( output_file, "%s/nonqcdsub-%s.pdf", output_dir, file_end ) ;
      can -> SaveAs( output_file ) ;

      h_sig0 -> Draw( "same hist" ) ;
      h_nonqcdsub -> Draw("e1 same") ;
      legend_with_sig0 -> Draw() ;

      sprintf( output_file, "%s/nonqcdsub-%s-with-sig0.pdf", output_dir, file_end ) ;
      can -> SaveAs( output_file ) ;


   } // draw_nonqcdsub

  //---------


   void draw_pull( const char* bin_name_substring ) {

      char hname[1000] ;

      TString clean_substring( bin_name_substring ) ;
      clean_substring.ReplaceAll("-","") ;
      clean_substring.ReplaceAll("_","") ;

      char histend[100] ;
      if ( strlen( bin_name_substring ) > 0 ) {
         sprintf( histend, "_%s", clean_substring.Data() ) ;
      } else {
         sprintf( histend, "" ) ;
      }

      sprintf( hname, "h_pull_zl%s", histend ) ;
      TH1F* h_pull = get_hist( hname ) ;

      TLegend* legend = new TLegend( 0.82, 0.35, 0.98, 0.90 ) ;
      legend -> AddEntry( h_pull, "Fit pull" ) ;

      h_pull -> SetFillColor(11) ;
      h_pull -> SetMinimum(-5) ;
      h_pull -> SetMaximum( 5) ;

      h_pull -> Draw("hist") ;
      h_pull -> Draw("axis hist same") ;
      h_pull -> Draw("axig hist same") ;
      gPad -> SetGridy(1) ;
      legend -> Draw() ;

      char file_end[100] ;
      if ( strlen( bin_name_substring ) > 0 ) {
         sprintf( file_end, "%s", bin_name_substring ) ;
      } else {
         sprintf( file_end, "allbins" ) ;
      }

      char output_file[10000] ;
      sprintf( output_file, "%s/pull-%s.pdf", output_dir, file_end ) ;
      can -> SaveAs( output_file ) ;


   } // draw_pull

  //---------




