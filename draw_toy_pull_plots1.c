
#include "TString.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TStyle.h"

#include <fstream>

   float get_true_val( const char* pname, const char* toy_dir ) ;

  //--------

   void draw_toy_pull_plots1( const char* toy_dir = "outputfiles/kqcd-fakedata-toys" ) {

      gStyle -> SetLabelSize( 0.09, "x" ) ;
      gStyle -> SetLabelSize( 0.08, "y" ) ;
      gStyle -> SetNdivisions( 505, "x" ) ;
      gStyle -> SetNdivisions( 505, "y" ) ;
      gStyle -> SetPadBottomMargin( 0.20 ) ;
      gStyle -> SetTitleSize( 0.09, "x" ) ;
      gStyle -> SetStatW(0.30) ;
      gStyle -> SetStatH(0.24) ;

      bool perfect_closure(true) ;

      char command[10000] ;

      sprintf( command, "grep \"perfect closure = TRUE\" %s/make-fakedata-config.txt", toy_dir ) ;
      int status = gSystem -> Exec( command ) ;
      if ( status == 0 ) {
         printf("\n\n Toys prepared with perfect closure.\n\n") ;
      } else {
         printf("\n\n Toys NOT prepared with perfect closure.\n\n") ;
      }

      printf("\n\n Making list of toy fit result files.\n") ;
      sprintf( command, "ls -1 %s/qcdlhfit-results-ws-toy*/kqcd-parameter-fit-results.txt > toyfitfiles.txt", toy_dir ) ;
      gSystem -> Exec( command ) ;

      ifstream ifs ;
      ifs.open( "toyfitfiles.txt" ) ;
      int fi(0) ;

      TH1F* hp_val[100] ;
      TH1F* hp_err[100] ;
      TH1F* hp_pull[100] ;
      float true_val[100] ;
      for ( int i=0; i<100; i++ ) { hp_val[i] = 0x0 ; hp_err[i] = 0x0 ; hp_pull[i] = 0x0 ; }
      int n_pars(0) ;

      while ( ifs.good() ) {

         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( !ifs.good() ) break ;
         printf(" file %4d : %s\n", fi, ts.Data() ) ;

         ifstream ifs_fr ;
         ifs_fr.open( ts.Data() ) ;
         if ( !ifs_fr.good() ) {
            printf("\n\n *** Problem opening file : %s\n\n", ts.Data() ) ; return ;
         }

         int pi(0) ;
         while ( ifs_fr.good() ) {

            TString ts_fr ;
            ts_fr.ReadLine( ifs_fr ) ;
            if ( !ifs_fr.good() ) break ;

            bool constrained(false) ;
            if ( ts_fr.Index( "constrained" ) >= 0 ) { constrained = true ; }

            char pname[100] ;
            float val(-1.) ;
            float err(-1.) ;
            float mean(-1.) ;
            float sigma(-1.) ;

            if ( constrained ) {
               sscanf( ts_fr.Data(), " %s %f constrained by %f +/- %f", pname, &val, &mean, &sigma ) ;
               printf("  %30s is constrained.    fit %6.4f  mean %6.4f,  sigma %6.4f\n", pname, val, mean, sigma ) ;
            } else {
               sscanf( ts_fr.Data(), " %s %f +/- %f", pname, &val, &err ) ;
               printf("  %30s is unconstrained.  fit %6.4f   +/-   %6.4f\n", pname, val, err ) ;
            }

            if ( fi == 0 ) {

               char hname_val[100] ;
               char htitle_val[1000] ;
               char hname_err[100] ;
               char htitle_err[1000] ;
               char hname_pull[100] ;
               char htitle_pull[1000] ;
               float hmin_val, hmax_val ;
               float hmin_err, hmax_err ;

               if ( constrained ) {
                  sprintf( hname_val, "h_fit_val_%s", pname ) ;
                  sprintf( htitle_val, "%s, fit value, constrained, %.3f +/- %.3f", pname, mean, sigma ) ;
                  hmin_val = mean - 5*sigma ;
                  hmax_val = mean + 5*sigma ;
                  sprintf( hname_pull, "h_fit_pull_%s", pname ) ;
                  sprintf( htitle_pull, "%s, fit pull, constrained", pname ) ;
               } else {
                  sprintf( hname_val, "h_fit_val_%s", pname ) ;
                  sprintf( htitle_val, "%s, fit value, unconstrained", pname ) ;
                  hmin_val = val - 5*err ;
                  hmax_val = val + 7*err ;
                  sprintf( hname_err, "h_fit_err_%s", pname ) ;
                  sprintf( htitle_err, "%s, fit error, constrained", pname ) ;
                  hmin_err =  0. ;
                  hmax_err = 2.5*err ;
                  sprintf( hname_pull, "h_fit_pull_%s", pname ) ;
                  sprintf( htitle_pull, "%s, fit pull, unconstrained", pname ) ;
               }
               if ( hmin_val < 0 ) hmin_val = 0. ;
               hp_val[pi] = new TH1F( hname_val, htitle_val, 40, hmin_val, hmax_val ) ;
               hp_val[pi] -> SetXTitle( pname ) ;
               if ( !constrained ) { hp_err[pi] = new TH1F( hname_err, htitle_err, 40, hmin_err, hmax_err ) ; }
               hp_pull[pi] = new TH1F( hname_pull, htitle_pull, 40, -6., 6. ) ;
               char xtitle[1000] ;
               if ( constrained ) {
                  sprintf( xtitle, "%s (val-mean)/sigma", pname ) ;
               } else {
                  sprintf( xtitle, "%s pull", pname ) ;
               }
               hp_pull[pi] -> SetXTitle( xtitle ) ;

               true_val[pi] = get_true_val( pname, toy_dir ) ;

               n_pars ++ ;
            }

            hp_val[pi] -> Fill( val ) ;
            if ( !constrained ) hp_err[pi] -> Fill( err ) ;

            float pull (-9.) ;
            if ( !constrained ) {
               if ( err > 0 ) pull = (val-true_val[pi])/err ;
            } else {
               if ( sigma > 0 ) pull = (val-true_val[pi])/sigma ;
            }
            hp_pull[pi] -> Fill( pull ) ;

            pi++ ;

         }

         fi++ ;

      }


      char pdf_file[10000] ;

   //-----

      TCanvas* can_err = new TCanvas( "can_err", "Parameter fit errors", 900, 1300 ) ;
      can_err -> Divide(2,5) ;

      for ( int pi=0; pi<n_pars; pi++ ) {

         can_err -> cd(pi+1) ;

         if ( hp_err[pi] == 0x0 ) continue ;

         hp_err[pi] -> SetFillColor(11) ;

         hp_err[pi] -> Draw() ;
         hp_err[pi] -> Draw("axis same") ;

      } // pi

      sprintf( pdf_file, "%s/fit-results-err.pdf", toy_dir ) ;
      can_err -> SaveAs( pdf_file ) ;


   //-----

      TCanvas* can_pull = new TCanvas( "can_pull", "Parameter pulls", 900, 1300 ) ;
      can_pull -> Divide(2,5) ;

      for ( int pi=0; pi<n_pars; pi++ ) {

         can_pull -> cd(pi+1) ;

         if ( hp_pull[pi] == 0x0 ) continue ;

         hp_pull[pi] -> SetFillColor(11) ;

         hp_pull[pi] -> Draw() ;
         hp_pull[pi] -> Draw("axis same") ;

      } // pi

      sprintf( pdf_file, "%s/fit-results-pull.pdf", toy_dir ) ;
      can_pull -> SaveAs( pdf_file ) ;

    //-----

      TCanvas* can_val = new TCanvas( "can_val", "Parameter values", 900, 1300 ) ;
      can_val -> Divide(2,5) ;

      for ( int pi=0; pi<n_pars; pi++ ) {

         can_val -> cd(pi+1) ;

         if ( hp_val[pi] == 0x0 ) continue ;

         hp_val[pi] -> SetFillColor(11) ;

         hp_val[pi] -> Draw() ;
         hp_val[pi] -> Draw("axis same") ;

         TArrow* ta = new TArrow() ;
         ta -> SetLineColor(4) ;
         ta -> SetLineWidth(3) ;
         ta -> DrawArrow( true_val[pi], 0.7*(hp_val[pi] -> GetMaximum()), true_val[pi], 0., 0.02 ) ;

      } // pi

      sprintf( pdf_file, "%s/fit-results-val.pdf", toy_dir ) ;
      can_val -> SaveAs( pdf_file ) ;


   } // draw_toy_pull_plots1

  //==========================================================================

   float get_true_val( const char* pname, const char* toy_dir ) {

      char command[10000] ;
      sprintf( command, "ls -1 %s/kqcd-fitconfig*.txt", toy_dir ) ;
      TString fitconfig_file = gSystem -> GetFromPipe( command ) ;

      ifstream ifs ;
      ifs.open( fitconfig_file.Data() ) ;
      if ( !ifs.good() ) {
         printf("\n\n *** Problem finding true values.\n\n") ;
         gSystem -> Exit(-1) ;
      }

      char target_name[100] ;
      if ( strcmp( pname, "Kqcd_ht1" ) == 0 ) {
         sprintf( target_name, "QCD-Kht1" ) ;
      } else if ( strcmp( pname, "Kqcd_ht2" ) == 0 ) {
         sprintf( target_name, "QCD-Kht2" ) ;
      } else if ( strcmp( pname, "Kqcd_ht3" ) == 0 ) {
         sprintf( target_name, "QCD-Kht3" ) ;
      } else if ( strcmp( pname, "Kqcd_mht2" ) == 0 ) {
         sprintf( target_name, "QCD-Kmht2" ) ;
      } else if ( strcmp( pname, "Kqcd_mht3" ) == 0 ) {
         sprintf( target_name, "QCD-Kmht3" ) ;
      } else if ( strcmp( pname, "Kqcd_mht4" ) == 0 ) {
         sprintf( target_name, "QCD-Kmht4" ) ;
      } else if ( strcmp( pname, "Kqcd_njet2" ) == 0 ) {
         sprintf( target_name, "QCD-Knjet2" ) ;
      } else if ( strcmp( pname, "Kqcd_njet3" ) == 0 ) {
         sprintf( target_name, "QCD-Knjet3" ) ;
      } else if ( strcmp( pname, "Kqcd_njet4" ) == 0 ) {
         sprintf( target_name, "QCD-Knjet4" ) ;
      } else if ( strcmp( pname, "Kqcd_njet5" ) == 0 ) {
         sprintf( target_name, "QCD-Knjet5" ) ;
      } else {
         printf("\n\n *** Unknown pname: %s\n\n", pname ) ;
         gSystem -> Exit(-1) ;
      }

      float ret_val(-9.) ;
      while ( ifs.good() ) {
         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( !ifs.good() ) break ;
         if ( ts.Index( target_name ) >=0 ) {
            char tpname[100] ;
            float err ;
            sscanf( ts.Data(), "%s %f %f", tpname, &ret_val, &err ) ;
            printf("  Found true value of %s to be %f\n", pname, ret_val ) ;
         }
      }

      return ret_val ;

   } // get_true_val

  //==========================================================================



