


#include "TROOT.h"

#include "TText.h"
#include "TLatex.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TMinuit.h"
#include "TMatrixT.h"

#include "histio.c"

#include <iostream>
#include <fstream>

  using std::cout ;
  using std::endl ;


   double data_Rqcd[5][5][5][5] ;
   double data_Rqcd_err[5][5][5][5] ;

   double fit_Rqcd_HT[5] ;
   double fit_SFqcd_MET[5] ;
   double fit_SFqcd_nb[5] ;
   double fit_SFqcd_njet[5] ;

   const int nBinsMET(4) ;
   const int nBinsHT(3) ;
   const int nBinsBjets(1) ;
   const int nBinsNjets(4) ;

   double calc_fit_error( TMinuit* tm, int hbi, int mbi, int bbi, int nji ) ;


  //-----------

   void minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

      int idummy = npar ;
      double fdummy = gin[0] ;
      fdummy = 0. ;
      idummy = iflag ;

      f = 0. ;

      //--- unpack the stupid par vector.
      int parind(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         fit_Rqcd_HT[hbi] = par[parind] ;
         parind ++ ;
      } // hbi.
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         if ( mbi == 0 ) {
            fit_SFqcd_MET[mbi] = 1.0 ;
         } else {
            fit_SFqcd_MET[mbi] = par[parind] ;
            parind++ ;
         }
      } // mbi.
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         if ( bbi == 0 ) {
            fit_SFqcd_nb[bbi] = 1.0 ;
         } else {
            fit_SFqcd_nb[bbi] = par[parind] ;
            parind++ ;
         }
      } // bbi.
      for ( int nji=0; nji<nBinsNjets; nji++ ) {
         if ( nji == 0 ) {
            fit_SFqcd_njet[nji] = 1.0 ;
         } else {
            fit_SFqcd_njet[nji] = par[parind] ;
            parind++ ;
         }
      } // nji.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               for ( int nji=0; nji<nBinsNjets; nji++ ) {
                  if ( !( data_Rqcd[mbi][hbi][bbi][nji] > 0. && data_Rqcd_err[mbi][hbi][bbi][nji] > 0. ) ) { continue ; }
                  double delta = data_Rqcd[mbi][hbi][bbi][nji] - fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] * fit_SFqcd_njet[nji] ;
                  f += delta*delta / (data_Rqcd_err[mbi][hbi][bbi][nji] * data_Rqcd_err[mbi][hbi][bbi][nji] ) ;
               } // nji.
            } // bbi.
         } // hbi.
      } // mbi.

   } // minuit_fcn.

  //-----------

  //=============================================================================================================================

   void modelfit_nb_njet4( const char* infile = "outputfiles/minDeltaPhiN-ratio-histograms-cut4.00-nb-njet.root", bool use_dphin = true ) {


      gDirectory->Delete("h*") ;

      loadHist( infile ) ;

      char hname[1000] ;
      char htitle[1000] ;


      TH1F* h_ratio_nb_njet[4][4] ;
      for ( int nji=0; nji<nBinsNjets; nji++ ) {
         for ( int nbi=0; nbi<nBinsBjets; nbi++ ) {
            if ( use_dphin ) {
               sprintf( hname, "h_ratio_minDeltaPhiN_nb%d_njet%d", nbi, nji+1 ) ;
            } else {
               sprintf( hname, "h_ratio_minDeltaPhi_nb%d_njet%d", nbi, nji+1 ) ;
            }
            h_ratio_nb_njet[nbi][nji] = (TH1F*) gDirectory -> FindObject( hname ) ;
            if ( h_ratio_nb_njet[nbi][nji] == 0x0 ) {
               printf( "\n\n *** Histogram %s missing from file %s\n\n", hname, infile ) ;
               return ;
            }
         } // nbi
      } // nji


      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               for ( int nji=0; nji<nBinsNjets; nji++ ) {

               int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

               double val, err ;

               val = h_ratio_nb_njet[bbi][nji] -> GetBinContent( histbin ) ;
               err = h_ratio_nb_njet[bbi][nji] -> GetBinError( histbin ) ;

               data_Rqcd[mbi][hbi][bbi][nji] = val ;
               data_Rqcd_err[mbi][hbi][bbi][nji] = err ;

               } // nji.
            } // bbi.
         } // hbi.
      } // mbi.


      int n_minuit_pars = nBinsHT + nBinsMET-1 + nBinsBjets-1 + nBinsNjets-1 ;

      TMinuit *myMinuit = new TMinuit( n_minuit_pars ) ; // arg is # of parameters

      myMinuit->SetFCN( minuit_fcn ) ;


      Double_t arglist[10] ;
      Int_t ierflg = 0 ;

      arglist[0] = 1 ;
      myMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //--- do this for chi2 fit.

      int parind(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         sprintf( pname, "Rqcd_HT%d", hbi+1 ) ;
         myMinuit->mnparm( parind, pname, data_Rqcd[0][hbi][0][0], 0.03, 0., 2., ierflg ) ;
         parind++ ;
      } // hbi.
      for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
         char pname[1000] ;
         sprintf( pname, "SFqcd_MET%d", mbi+1 ) ;
         myMinuit->mnparm( parind, pname, 1.0, 0.10, 0., 4., ierflg ) ;
         parind++ ;
      } // mbi.
      for ( int bbi=1; bbi<nBinsBjets; bbi++ ) {
         char pname[1000] ;
         sprintf( pname, "SFqcd_nb%d", bbi+1 ) ;
         myMinuit->mnparm( parind, pname, 1.0, 0.10, 0., 4., ierflg ) ;
         parind++ ;
      } // mbi.
      for ( int nji=1; nji<nBinsNjets; nji++ ) {
         char pname[1000] ;
         sprintf( pname, "SFqcd_njet%d", nji+1 ) ;
         myMinuit->mnparm( parind, pname, 1.0, 0.10, 0., 4., ierflg ) ;
         parind++ ;
      } // nji.

      myMinuit->Migrad() ;
      myMinuit->mncomd("hesse",ierflg) ;


      printf("\n\n") ;
      parind = 0 ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "Rqcd_HT%d", hbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_Rqcd_HT[hbi] = val ;
         parind++ ;
      } // hbi.
      for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "SFqcd_MET%d", mbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_SFqcd_MET[mbi] = val ;
         parind++ ;
      } // mbi.
      for ( int bbi=1; bbi<nBinsBjets; bbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "SFqcd_nb%d", bbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_SFqcd_nb[bbi] = val ;
         parind++ ;
      } // mbi.
      for ( int nji=1; nji<nBinsNjets; nji++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "SFqcd_njet%d", nji+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_SFqcd_njet[nji] = val ;
         parind++ ;
      } // nji.
      printf("\n\n") ;


//            //   Double_t cov_mat[n_minuit_pars][n_minuit_pars] ;
//            //   myMinuit -> mnemat( &cov_mat[0][0], n_minuit_pars ) ;
//            //   printf("\n\n  ===== Covariance matrix:\n") ;
//            //   for ( int i=0; i<10; i++ ) {
//            //      printf( "  %2d : ", i+1 ) ;
//            //      for ( int j=0; j<10; j++ ) {
//            //         printf( "  %10.7f  ", cov_mat[i][j] ) ;
//            //      }
//            //      printf("\n") ;
//            //   }
//            //   printf("\n\n") ;




      TH1F* h_ratio_model[nBinsBjets][nBinsNjets] ;
      TH1F* h_ratio_model_hist[nBinsBjets][nBinsNjets] ;

      for ( int nji=0; nji<nBinsNjets; nji++ ) {
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
            sprintf( hname, "h_ratio_model_%db_njet%d", bbi+1, nji+1 ) ;
            h_ratio_model[bbi][nji] = (TH1F*) h_ratio_nb_njet[bbi][nji]->Clone( hname ) ;
            sprintf( hname, "h_ratio_model_%db_njet%d_hist", bbi+1, nji+1 ) ;
            h_ratio_model_hist[bbi][nji] = (TH1F*) h_ratio_nb_njet[bbi][nji]->Clone( hname ) ;
            h_ratio_model_hist[bbi][nji] -> SetLineWidth(2) ;
            h_ratio_model_hist[bbi][nji] -> SetLineColor(kRed+1) ;
            h_ratio_model[bbi][nji] -> SetFillColor(kRed-10) ;
            h_ratio_model[bbi][nji] -> SetMarkerStyle(0) ;
            h_ratio_model[bbi][nji] -> SetLabelSize( 0.06, "x" ) ;
            h_ratio_model[bbi][nji] -> SetLabelOffset( 0.017, "x" ) ;
         } // bbi.
      } // nji.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               for ( int nji=0; nji<nBinsNjets; nji++ ) {

                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

                  double model_val = fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] * fit_SFqcd_njet[nji] ;
                  double model_err = calc_fit_error( myMinuit, hbi, mbi, bbi, nji ) ;

                  h_ratio_model_hist[bbi][nji] -> SetBinContent( histbin, model_val ) ;
                  h_ratio_model[bbi][nji] -> SetBinContent( histbin, model_val ) ;
                  h_ratio_model[bbi][nji] -> SetBinError( histbin, model_err ) ;

               } // nji.
            } // bbi.
         } // hbi.
      } // mbi.


//             //  for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
//             //     printf("\n\n") ;
//             //     h_ratio_model[bbi] -> Print("all") ;
//             //     printf("\n\n") ;
//             //  } // bbi



      gStyle -> SetPadBottomMargin( 0.24 ) ;
      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1_modelfit_nb_njet2" ) ;
      if ( can1 == 0x0 ) can1 = new TCanvas( "can1_modelfit_nb_njet2", "QCD Model fit", 800, 800 ) ;
      can1 -> Clear() ;
      can1 -> Divide(nBinsNjets,nBinsBjets) ;
      ///can1 -> Divide(nBinsBjets,nBinsNjets) ;

      int can_index(1) ;

      for ( int nji=0; nji<nBinsNjets; nji++ ) {
         float hist_draw_max = 0.7 ;
         //if ( nji == 0 ) hist_draw_max = 0.70 ; //formerly 0.6
         //if ( nji == 1 ) hist_draw_max = 0.35 ; //formerly 0.25
         //if ( nji == 2 ) hist_draw_max = 0.35 ; //formerly 0.25

         h_ratio_model[0][nji] -> SetMaximum( hist_draw_max ) ;

         can1 -> cd(can_index) ;
         h_ratio_model[0][nji] -> Draw( "e2" ) ;
         h_ratio_model_hist[0][nji] -> Draw( "hist same" ) ;
         h_ratio_nb_njet[0][nji] -> Draw( "same" ) ;
         gPad -> SetGridy(1) ;
         can_index ++ ;

      }

     ////Dump txt file of event counts
     
     string inf ( infile ) ;
     unsigned long pos = inf.find_last_of(".root");
     string ret;
     if ( pos != std::string::npos ) {
       ret = inf.substr(0,pos-4);
     }
     else {
       ret = "./";
     }
     TString fitplot = ret;fitplot += "-fitplot.pdf"; 
     can1->SaveAs((TString)fitplot);


     ///
      ofstream ofs;
      string forofs  = ret;
      forofs += "-table.txt";
      ofs.open(forofs.c_str());
      char toofs[2000];

      TH1F* h_all_nb_njet[4][4] ;
      TH1F* h_pass_nb_njet[4][4] ;
      char njetbin_str[100] ;
      for ( int nji=0; nji<nBinsNjets; nji++ ) { 
         for ( int nbi=0; nbi<nBinsBjets; nbi++ ) { 
            sprintf( njetbin_str, "_njet%d", nji+1);      
      
            sprintf( hname, "h_qcd_yield_all_nb%d_%s%s", nbi, "minDeltaPhiN", njetbin_str ) ;
            h_all_nb_njet[nbi][nji] = (TH1F*) gDirectory -> FindObject ( hname );
            if ( h_all_nb_njet[nbi][nji] == 0x0 ) { 
               printf( "\n\n *** Histogram %s missing from file %s\n\n", hname, infile ) ; 
               return ;
            }   
            sprintf( hname, "h_qcd_yield_pass_nb%d_%s%s", nbi, "minDeltaPhiN", njetbin_str ) ;
            h_pass_nb_njet[nbi][nji] = (TH1F*) gDirectory -> FindObject ( hname );
            if ( h_all_nb_njet[nbi][nji] == 0x0 ) { 
               printf( "\n\n *** Histogram %s missing from file %s\n\n", hname, infile ) ; 
               return ;
            }   

         } // nbi
      } // nji

     //printf("Bin Label          ||      Nlow +/- error   :    Nhigh +/- error    : Nhi/Nlo +/- error  : RatioMC +/- err   : Fit Ratio +/- err\n\n");
     sprintf(toofs,"Bin Label          ||      Nlow +/- error   :    Nhigh +/- error    : Nhi/Nlo +/- error  : RatioMC +/- err   : Fit Ratio +/- err\n\n");
     ofs << toofs;

     for ( int bbi=0; bbi<nBinsBjets; bbi++ ) { 
       for ( int nji=0; nji<nBinsNjets; nji++ ) { 
         for ( int mbi=0; mbi<nBinsMET; mbi++ ) { 
           for ( int hbi=0; hbi<nBinsHT; hbi++ ) { 

                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ; 

                  double model_val = fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] * fit_SFqcd_njet[nji] ;
                  double model_err = calc_fit_error( myMinuit, hbi, mbi, bbi, nji ) ; 
                  float ratio = h_ratio_nb_njet[bbi][nji] -> GetBinContent(histbin);
                  float ratioerr = h_ratio_nb_njet[bbi][nji] -> GetBinError(histbin);
                
                  float nhi = h_pass_nb_njet[bbi][nji]->GetBinContent(histbin);
                  float nall = h_all_nb_njet[bbi][nji]->GetBinContent(histbin);
                  float nhierr = h_pass_nb_njet[bbi][nji]->GetBinError(histbin);
                  float nallerr = h_all_nb_njet[bbi][nji]->GetBinError(histbin);
                  float nlo = nall - nhi;
                  float nloerr = sqrt(nallerr*nallerr - nhierr*nhierr);
                
                  float myratio = 0;
                  float myratioerr = 0;
                  if ( nlo > 0. ) myratio = nhi/nlo;
                  if ( nlo > 0 && nhi > 0 ) myratioerr = myratio * sqrt( pow(nloerr/nlo,2) + pow(nhierr/nhi,2) ) ;

                  TString label = h_ratio_model[bbi][nji] -> GetXaxis() -> GetBinLabel(histbin);
                  label += " nB"; label += bbi;


                  //printf("%s || %9.4f +/- %7.3f : %9.4f +/- %7.3f : %7.4f +/- %6.3f : %6.4f +/- %6.3f : %6.4f +/- %6.3f \n",label.Data(),nlo,nloerr,nhi,nhierr,myratio,myratioerr,ratio,ratioerr,model_val,model_err);
                  sprintf(toofs,"%s || %9.4f +/- %7.3f : %9.4f +/- %7.3f : %7.4f +/- %6.3f : %6.4f +/- %6.3f : %6.4f +/- %6.3f \n",label.Data(),nlo,nloerr,nhi,nhierr,myratio,myratioerr,ratio,ratioerr,model_val,model_err);
                  ofs << toofs ;

               } // nji.
            } // bbi.
         } // hbi.
      } // mbi.


      ofs.close();
     ///


   } // modelfit_nb_njet4


  //==========================================================================================

   double calc_fit_error( TMinuit* tm, int hbi, int mbi, int bbi, int nji ) {

      //bool verb = true ;
      bool verb = false ;

      if ( verb ) {
         printf("\n\n ==================================================================================================\n") ;
         printf("     calc_fit_error: computing model error for HT%d, MET%D, NB%d, Njet%d\n\n", hbi+1, mbi+1, bbi, nji+1 ) ;
      }

      if ( tm == 0x0 ) return -1 ;

      int n_minuit_pars = tm -> GetNumPars() ;
      if ( n_minuit_pars <= 0 ) return -1 ;
      if (verb) printf(" Number of minuit parameters: %d\n", n_minuit_pars ) ;

      Double_t cov_mat[n_minuit_pars][n_minuit_pars] ;
      tm -> mnemat( &cov_mat[0][0], n_minuit_pars ) ;


     //--- Minuit parameter indices, counting from 0
      int ht_pind = hbi ;
      int met_pind = -1 ;
      if ( mbi > 0 ) met_pind = nBinsHT + mbi - 1 ;
      int nb_pind = -1 ;
      if ( bbi > 0 ) nb_pind = nBinsHT + (nBinsMET-1) + bbi - 1 ;
      int njet_pind = -1 ;
      if ( nji > 0 ) njet_pind = nBinsHT + (nBinsMET-1) + (nBinsBjets-1) + nji - 1 ;

      if ( verb ) printf("  Minuit parameter indices:  hbi=%d, pi=%d;   mbi=%d, pi=%d;   bbi=%d, pi=%d;   nji=%d, pi=%d\n",
        hbi, ht_pind, mbi, met_pind, bbi, nb_pind, nji, njet_pind ) ;

      double ht_par_val, ht_par_err ;
      tm -> GetParameter( ht_pind, ht_par_val, ht_par_err ) ;
      double met_par_val, met_par_err ;
      if ( mbi > 0 ) {
         tm -> GetParameter( met_pind, met_par_val, met_par_err ) ;
      } else {
         met_par_val = 1 ;  met_par_err = 0. ;
      }
      double nb_par_val, nb_par_err ;
      if ( bbi > 0 ) {
         tm -> GetParameter( nb_pind, nb_par_val, nb_par_err ) ;
      } else {
         nb_par_val = 1 ;  nb_par_err = 0. ;
      }
      double njet_par_val, njet_par_err ;
      if ( nji > 0 ) {
         tm -> GetParameter( njet_pind, njet_par_val, njet_par_err ) ;
      } else {
         njet_par_val = 1 ;  njet_par_err = 0. ;
      }

      if ( verb ) printf("   HT par: %5.3f +/- %5.3f\n", ht_par_val, ht_par_err ) ;
      if ( verb ) printf("  MET par: %5.3f +/- %5.3f\n", met_par_val, met_par_err ) ;
      if ( verb ) printf("   Nb par: %5.3f +/- %5.3f\n", nb_par_val, nb_par_err ) ;
      if ( verb ) printf(" Njet par: %5.3f +/- %5.3f\n", njet_par_val, njet_par_err ) ;




     //--- Calculate the three partial derivatives.
      double df_dhtpar   = met_par_val  *  nb_par_val * njet_par_val ;
      double df_dmetpar  =  ht_par_val  *  nb_par_val * njet_par_val ;
      double df_dnbpar   =  ht_par_val  * met_par_val * njet_par_val ;
      double df_dnjetpar =  ht_par_val  * met_par_val * nb_par_val   ;


     //--- Create the partial derivative vector.
      TMatrixT<double> pd_col_vec( n_minuit_pars, 1 ) ;
      TMatrixT<double> pd_row_vec( 1, n_minuit_pars ) ;
      for ( int i=0; i<n_minuit_pars; i++ ) {
         pd_col_vec(i,0) = 0. ;
         pd_row_vec(0,i) = 0. ;
      }
      pd_col_vec( ht_pind, 0 ) = df_dhtpar ;
      pd_row_vec( 0, ht_pind ) = df_dhtpar ;
      if ( mbi > 0 ) {
         pd_col_vec( met_pind, 0 ) = df_dmetpar ;
         pd_row_vec( 0, met_pind ) = df_dmetpar ;
      }
      if ( bbi > 0 ) {
         pd_col_vec( nb_pind, 0 ) = df_dnbpar ;
         pd_row_vec( 0, nb_pind ) = df_dnbpar ;
      }
      if ( nji > 0 ) {
         pd_col_vec( njet_pind, 0 ) = df_dnjetpar ;
         pd_row_vec( 0, njet_pind ) = df_dnjetpar ;
      }
      if ( verb ) {
         printf("   Partial derivative column vector:\n") ;
         pd_col_vec.Print() ;
      }
      if ( verb ) {
         printf("   Partial derivative row vector:\n") ;
         pd_row_vec.Print() ;
      }

      TMatrixT<double> cov_mat_tm( n_minuit_pars, n_minuit_pars ) ;
      for ( int i=0; i<n_minuit_pars; i++ ) {
         for ( int j=0; j<n_minuit_pars; j++ ) {
            cov_mat_tm( i, j ) = cov_mat[i][j] ;
         } // j
      } // i

      if ( verb ) {
         printf("\n\n  ===== Covariance matrix:\n") ;
         cov_mat_tm.Print() ;
      }


      TMatrixT<double> cov_times_pd_col( n_minuit_pars, 1 ) ;
      cov_times_pd_col.Mult( cov_mat_tm, pd_col_vec ) ;
      if ( verb ) {
         printf( "\n\n  ===== Cov mat * pd col vec:\n") ;
         cov_times_pd_col.Print() ;
      }

      TMatrixT<double> pd_row_times_prod( 1, 1 ) ;
      pd_row_times_prod.Mult( pd_row_vec, cov_times_pd_col ) ;
      if ( verb ) {
         printf( "\n\n  ===== pd row vec * prod:\n") ;
         pd_row_times_prod.Print() ;
      }


      double fit_error(0.) ;

      if ( pd_row_times_prod(0,0) > 0 ) {
         fit_error = sqrt( pd_row_times_prod(0,0) ) ;
      }

      if ( verb ) printf("\n Final answer for fit error is %5.3f\n\n", fit_error ) ;

      return fit_error ;

   } // calc_fit_error

  //==========================================================================================








