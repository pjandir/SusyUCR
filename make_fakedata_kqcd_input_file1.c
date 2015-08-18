
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"

#include <fstream>

      float qcd_kht_val[10] ;
      float qcd_kht_err[10] ;

      float qcd_knb_val[10] ;
      float qcd_knb_err[10] ;

      float qcd_knjet_val[10] ;
      float qcd_knjet_err[10] ;

      float qcd_kmht_val[10] ;
      float qcd_kmht_err[10] ;

   float find_line_val( ifstream& ifs, const char* key ) ;
   bool  find_line_val_err( ifstream& ifs, const char* key, float& val, float& err ) ;

   void  read_kqcd_pars( const char* kqcd_fitconfig_file ) ;

  //------------

   void make_fakedata_kqcd_input_file1(
          bool perfect_qcd_closure = false,
          bool randomize_nobs = false,
          int  random_seed = 0,
          bool make_plots = true,
          const char* kqcd_output_fakedata_file = "outputfiles/kqcd-input-fakedata.txt",
          const char* kqcd_input_lostlep_file = "outputfiles/kqcd-input-lostlep.txt",
          const char* kqcd_input_hadtau_file = "outputfiles/kqcd-input-hadtau.txt",
          const char* kqcd_input_znunu_file = "outputfiles/kqcd-input-znunu.txt",
          const char* kqcd_input_qcd_file = "outputfiles/kqcd-input-qcdmc.txt",
          const char* kqcd_fitconfig_file = "outputfiles/kqcd-fitconfig-with-constraints.txt"
       ) {

      TRandom3 tran ;
      tran.SetSeed( random_seed ) ;

      if ( perfect_qcd_closure ) read_kqcd_pars( kqcd_fitconfig_file ) ;

      gStyle -> SetPadBottomMargin( 0.35 ) ;
      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadGridY(1) ;

      ifstream ifs_lostlep ;
      ifstream ifs_hadtau ;
      ifstream ifs_znunu ;
      ifstream ifs_qcd ;

      ifs_lostlep.open( kqcd_input_lostlep_file ) ;
      if ( !ifs_lostlep.good() ) { printf("\n\n *** Bad input lostlep file: %s\n\n", kqcd_input_lostlep_file ) ; return ; }
      ifs_hadtau.open( kqcd_input_hadtau_file ) ;
      if ( !ifs_hadtau.good() ) { printf("\n\n *** Bad input hadtau file: %s\n\n", kqcd_input_hadtau_file ) ; return ; }
      ifs_znunu.open( kqcd_input_znunu_file ) ;
      if ( !ifs_znunu.good() ) { printf("\n\n *** Bad input znunu file: %s\n\n", kqcd_input_znunu_file ) ; return ; }
      ifs_qcd.open( kqcd_input_qcd_file ) ;
      if ( !ifs_qcd.good() ) { printf("\n\n *** Bad input qcd file: %s\n\n", kqcd_input_qcd_file ) ; return ; }

      TH1F* h_kqcd_input_data_lowdphi = new TH1F( "h_kqcd_input_data_lowdphi", "KQCD fit input, low delta phi, data (fake)", 55, 0.5, 55.5 ) ;
      TH1F* h_kqcd_input_data_highdphi = new TH1F( "h_kqcd_input_data_highdphi", "KQCD fit input, high delta phi, data (fake)", 55, 0.5, 55.5 ) ;

      FILE* ofp ;
      if ( (ofp=fopen( kqcd_output_fakedata_file, "w" ) )==NULL ) { printf("\n\n *** Problem opening output file: %s\n\n", kqcd_output_fakedata_file ) ; return ; }

      while ( ifs_qcd.good() ) {

         TString ts ;
         ts.ReadLine( ifs_qcd ) ;
         if ( !ifs_qcd.good() ) { break ; }
         int qcd_bin_index ;
         char qcd_bin_name[100] ;
         float qcd_ldp_val, qcd_ldp_err, qcd_zl_val, qcd_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &qcd_bin_index, qcd_bin_name, &qcd_ldp_val, &qcd_ldp_err, &qcd_zl_val, &qcd_zl_err ) ;

         if ( perfect_qcd_closure ) {
            int fb_nji(-1), fb_hbi(-1), fb_mbi(-1) ;
            sscanf( qcd_bin_name, "FB-Njet%d-Nbsum-MHT%d-HT%d", &fb_nji, &fb_mbi, &fb_hbi ) ;
            printf( " %s, fb_nji=%d, fb_mbi=%d, fb_hbi=%d\n", qcd_bin_name, fb_nji, fb_mbi, fb_hbi ) ;
            if ( fb_nji<0 || fb_hbi<0 || fb_mbi<0 ) { printf("\n\n *** Bad bin name.\n\n\n") ; return ; }
            float qcd_ratio = qcd_kht_val[fb_hbi] * qcd_knjet_val[fb_nji] * qcd_kmht_val[fb_mbi] ;
            printf("  QCD ratio = Kqcd_ht * Kqcd_nj * Kqcd_mht = %.3f * %.3f * %.3f = %.4f\n", qcd_kht_val[fb_hbi], qcd_knjet_val[fb_nji], qcd_kmht_val[fb_mbi], qcd_ratio ) ;
            float model_qcd_zl_val = qcd_ldp_val * qcd_ratio ;
            printf("  QCD model prediction:  Nldp * R = %6.1f * %.4f = %6.1f\n", qcd_ldp_val, qcd_ratio, model_qcd_zl_val ) ;
            printf(" qcd     bin index %2d,  name %s : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f  (model ZL %6.1f)\n", qcd_bin_index, qcd_bin_name, qcd_ldp_val, qcd_ldp_err, qcd_zl_val, qcd_zl_err, model_qcd_zl_val ) ;
            qcd_zl_val = model_qcd_zl_val ;
         } else {
            printf(" qcd     bin index %2d,  name %s : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", qcd_bin_index, qcd_bin_name, qcd_ldp_val, qcd_ldp_err, qcd_zl_val, qcd_zl_err ) ;
         }




         ts.ReadLine( ifs_lostlep ) ;
         int lostlep_bin_index ;
         char lostlep_bin_name[100] ;
         float lostlep_ldp_val, lostlep_ldp_err, lostlep_zl_val, lostlep_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &lostlep_bin_index, lostlep_bin_name, &lostlep_ldp_val, &lostlep_ldp_err, &lostlep_zl_val, &lostlep_zl_err ) ;
         if ( lostlep_bin_index != qcd_bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", qcd_bin_index, lostlep_bin_index ) ; return ; }
         if ( strcmp( qcd_bin_name, lostlep_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", qcd_bin_name, lostlep_bin_name ) ; return ; }
         printf(" lostlep bin index %2d,  name %s : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", lostlep_bin_index, lostlep_bin_name, lostlep_ldp_val, lostlep_ldp_err, lostlep_zl_val, lostlep_zl_err ) ;

         ts.ReadLine( ifs_hadtau ) ;
         int hadtau_bin_index ;
         char hadtau_bin_name[100] ;
         float hadtau_ldp_val, hadtau_ldp_err, hadtau_zl_val, hadtau_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &hadtau_bin_index, hadtau_bin_name, &hadtau_ldp_val, &hadtau_ldp_err, &hadtau_zl_val, &hadtau_zl_err ) ;
         if ( hadtau_bin_index != qcd_bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", qcd_bin_index, hadtau_bin_index ) ; return ; }
         if ( strcmp( qcd_bin_name, hadtau_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", qcd_bin_name, hadtau_bin_name ) ; return ; }
         printf(" hadtau  bin index %2d,  name %s : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", hadtau_bin_index, hadtau_bin_name, hadtau_ldp_val, hadtau_ldp_err, hadtau_zl_val, hadtau_zl_err ) ;

         ts.ReadLine( ifs_znunu ) ;
         int znunu_bin_index ;
         char znunu_bin_name[100] ;
         float znunu_ldp_val, znunu_ldp_err, znunu_zl_val, znunu_zl_err ;
         sscanf( ts.Data(), "%d %s %f +/- %f %f +/- %f", &znunu_bin_index, znunu_bin_name, &znunu_ldp_val, &znunu_ldp_err, &znunu_zl_val, &znunu_zl_err ) ;
         if ( znunu_bin_index != qcd_bin_index ) { printf("\n\n *** Inconsistent bin indices: %d %d\n\n", qcd_bin_index, znunu_bin_index ) ; return ; }
         if ( strcmp( qcd_bin_name, znunu_bin_name ) != 0 ) { printf("\n\n *** Inconsistent bin names: %s %s\n\n", qcd_bin_name, znunu_bin_name ) ; return ; }
         printf(" znunu   bin index %2d,  name %s : LDP %7.1f +/- %5.1f,   ZL %6.1f +/- %5.1f\n", znunu_bin_index, znunu_bin_name, znunu_ldp_val, znunu_ldp_err, znunu_zl_val, znunu_zl_err ) ;

         float ldp_val = lostlep_ldp_val + hadtau_ldp_val + znunu_ldp_val + qcd_ldp_val ;
         float zl_val  = lostlep_zl_val  + hadtau_zl_val  + znunu_zl_val  + qcd_zl_val  ;

         if ( randomize_nobs ) {
            double mean ;
            mean = ldp_val ;
            ldp_val = tran.Poisson( mean ) ;
            mean = zl_val ;
            zl_val = tran.Poisson( mean ) ;
         }

         h_kqcd_input_data_lowdphi  -> SetBinContent( qcd_bin_index, ldp_val ) ;
         h_kqcd_input_data_highdphi -> SetBinContent( qcd_bin_index, zl_val ) ;

         char binlabel[100] ;
         sprintf( binlabel, "%2d %s", qcd_bin_index, qcd_bin_name ) ;
         h_kqcd_input_data_lowdphi  -> GetXaxis() -> SetBinLabel( qcd_bin_index, binlabel ) ;
         h_kqcd_input_data_highdphi -> GetXaxis() -> SetBinLabel( qcd_bin_index, binlabel ) ;

         fprintf( ofp, " %2d %30s  %6d  %6d\n", qcd_bin_index, qcd_bin_name, TMath::Nint( ldp_val ), TMath::Nint( zl_val ) ) ;

         printf("\n") ;

      } // reading qcd file

      fclose( ofp ) ;

      if ( !make_plots ) return ;

      h_kqcd_input_data_lowdphi  -> GetXaxis() -> LabelsOption( "v" ) ;
      h_kqcd_input_data_highdphi -> GetXaxis() -> LabelsOption( "v" ) ;

      h_kqcd_input_data_lowdphi -> SetMarkerStyle(20) ;
      h_kqcd_input_data_highdphi -> SetMarkerStyle(20) ;

      h_kqcd_input_data_lowdphi -> SetMinimum(0.1) ;
      h_kqcd_input_data_highdphi -> SetMinimum(0.1) ;


      char pdffile[10000] ;

      TCanvas* can = new TCanvas( "can_make_fakedata_input_files1", "Data inputs (fake)", 1100, 900 ) ;

      can -> Divide( 1, 2 ) ;

      can -> cd(1) ;
      h_kqcd_input_data_lowdphi -> Draw( "e" ) ;
      can -> cd(2) ;
      h_kqcd_input_data_highdphi -> Draw( "e" ) ;

      sprintf( pdffile, "outputfiles/kqcd-input-fakedata-liny.pdf" ) ;
      can -> Update() ; can -> Draw() ; can -> SaveAs( pdffile ) ;

     //----
      can -> cd(1) ;
      gPad -> SetLogy(1) ;
      can -> cd(2) ;
      gPad -> SetLogy(1) ;

      sprintf( pdffile, "outputfiles/kqcd-input-fakedata-logy.pdf" ) ;
      can -> Update() ; can -> Draw() ; can -> SaveAs( pdffile ) ;

     //----
      h_kqcd_input_data_lowdphi -> SetMaximum(20) ;
      h_kqcd_input_data_highdphi -> SetMaximum(20) ;

      can -> cd(1) ;
      gPad -> SetLogy(0) ;
      can -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/kqcd-input-fakedata-zoom1.pdf" ) ;
      can -> Update() ; can -> Draw() ; can -> SaveAs( pdffile ) ;

     //----
      h_kqcd_input_data_lowdphi -> SetMaximum(5) ;
      h_kqcd_input_data_highdphi -> SetMaximum(5) ;

      can -> cd(1) ;
      gPad -> SetLogy(0) ;
      can -> cd(2) ;
      gPad -> SetLogy(0) ;

      sprintf( pdffile, "outputfiles/kqcd-input-fakedata-zoom2.pdf" ) ;
      can -> Update() ; can -> Draw() ; can -> SaveAs( pdffile ) ;


      printf("\n\n Saving histograms to outputfiles/fakedata-input.root\n\n") ;
      TFile* tf_out = new TFile( "outputfiles/fakedata-input.root", "RECREATE" ) ;
      h_kqcd_input_data_lowdphi -> Write() ;
      h_kqcd_input_data_highdphi -> Write() ;
      tf_out -> Close() ;




   } // make_fakedata_input_files1


  //=======================================================================

   void  read_kqcd_pars( const char* kqcd_fitconfig_file ) {

      char pname[100] ;
      char pname2[100] ;

      ifstream ifs ;
      ifs.open( kqcd_fitconfig_file ) ;
      if ( !ifs.good() ) {
         printf("\n\n *** Problem with fitconfig file : %s\n\n", kqcd_fitconfig_file ) ; return ;
      }

      printf("\n QCD model parameters\n") ;
      int n_qcd_kht_pars = find_line_val( ifs, "N-QCD-Kht-pars" ) ;
      for ( int i=1; i<=n_qcd_kht_pars; i++ ) {
         sprintf( pname, "QCD-Kht%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_kht_val[i], qcd_kht_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_kht_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_kht_val[i] ) ;
         } else if ( qcd_kht_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_kht_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_kht_val[i], qcd_kht_err[i] ) ;
         }
      } // i

      int n_qcd_kmht_pars = find_line_val( ifs, "N-QCD-Kmht-pars" ) ;
      for ( int i=1; i<=n_qcd_kmht_pars; i++ ) {
         sprintf( pname, "QCD-Kmht%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_kmht_val[i], qcd_kmht_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_kmht_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_kmht_val[i] ) ;
         } else if ( qcd_kmht_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_kmht_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_kmht_val[i], qcd_kmht_err[i] ) ;
         }
      } // i

      int n_qcd_knjet_pars = find_line_val( ifs, "N-QCD-Knjet-pars" ) ;
      for ( int i=1; i<=n_qcd_knjet_pars; i++ ) {
         sprintf( pname, "QCD-Knjet%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_knjet_val[i], qcd_knjet_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_knjet_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_knjet_val[i] ) ;
         } else if ( qcd_knjet_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_knjet_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_knjet_val[i], qcd_knjet_err[i] ) ;
         }
      } // i


      int n_qcd_knb_pars = find_line_val( ifs, "N-QCD-Knb-pars" ) ;
      for ( int i=1; i<=n_qcd_knb_pars; i++ ) {
         sprintf( pname, "QCD-Knb%d", i ) ;
         if ( !find_line_val_err( ifs, pname, qcd_knb_val[i], qcd_knb_err[i] ) ) { printf("\n\n *** bad input\n\n") ; return ; }
         if ( qcd_knb_err[i] < 0. ) {
            printf("  %s : %5.3f  floating\n", pname, qcd_knb_val[i] ) ;
         } else if ( qcd_knb_err[i] == 0. ) {
            printf("  %s : %5.3f  fixed\n", pname, qcd_knb_val[i] ) ;
         } else {
            printf("  %s : %5.3f +/- %5.3f (constrained)\n", pname, qcd_knb_val[i], qcd_knb_err[i] ) ;
         }
      } // i

   } // read_kqcd_pars


  //======================================================================================================

  //=================================================================================

   float find_line_val( ifstream& ifs, const char* key ) {

      ifs.seekg(0) ;
      TString ts ;

      int lc(0) ;
      while ( ifs.good() ) {
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         /////////// printf(" %4d : %s\n", lc++, ts.Data() ) ;
         if ( strcmp( first->GetString().Data(), key ) == 0 ) {
            /////////// printf(" Found %s\n", key ) ;
            float val(-1) ;
            if ( tokens->GetEntries() > 1 ) {
               TObjString* second = (TObjString*) tokens->At(1) ;
               sscanf( second->GetString().Data(), "%f", &val ) ;
               return val ;
            } else {
               return -2 ;
            }
         }
      }

      return -1 ;

   } // find_line_val

  //=================================================================================

   bool find_line_val_err( ifstream& ifs, const char* key, float& val, float& err ) {

      val = 0. ;
      err = 0. ;

      ifs.seekg(0) ;
      TString ts ;

      int lc(0) ;
      while ( ifs.good() ) {
         ts.ReadLine( ifs ) ;
         TObjArray* tokens = ts.Tokenize(" ") ;
         TObjString* first = (TObjString*) tokens->At(0) ;
         /////////// printf(" %4d : %s\n", lc++, ts.Data() ) ;
         if ( strcmp( first->GetString().Data(), key ) == 0 ) {
            /////////// printf(" Found %s\n", key ) ;
            if ( tokens->GetEntries() >= 3 ) {
               TObjString* second = (TObjString*) tokens->At(1) ;
               sscanf( second->GetString().Data(), "%f", &val ) ;
               TObjString* third = (TObjString*) tokens->At(2) ;
               sscanf( third->GetString().Data(), "%f", &err ) ;
               return true ;
            }
         }
      }

      return false ;


   } // find_line_val_err

  //=================================================================================


