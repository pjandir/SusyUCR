
#include "TH1F.h"
#include "TAxis.h"

   void calc_0m_rms( TH1F* hp, double& rms_val, double& rms_err ) {

      rms_val = 0. ;
      rms_err = 0. ;

      if ( hp == 0x0 ) return ;

      double sum(0.) ;
      double integral(0.) ;

      int nbins = hp -> GetNbinsX() ;
      TAxis* xaxis = hp -> GetXaxis() ;

      for ( int bi=1; bi<nbins; bi++ ) { // don't include last bin, which has overflow added
         double bin_center = xaxis -> GetBinCenter( bi ) ;
         double sumw = hp -> GetBinContent( bi ) ;
         integral += sumw ;
         sum += sumw * pow( bin_center, 2 ) ;
      } // bi

      if ( integral <= 0 ) return ;

      double avesq = sum / integral ;
      rms_val = sqrt( avesq ) ;


      //-- Attempt to compute a reasonable stat error on this.
      double mssum(0.) ;
      double emssum(0.) ;
      for ( int bi=1; bi<nbins; bi++ ) { // don't include last bin, which has overflow added
         double bin_center = xaxis -> GetBinCenter( bi ) ;
         double sumw = hp -> GetBinContent( bi ) ;
         double sumwerr = hp -> GetBinError( bi ) ;
         mssum += (sumw * bin_center*bin_center / integral) ;
         emssum += pow( (bin_center*bin_center / integral)*sumwerr, 2 ) ;
      } // bi
      double ems = sqrt( emssum ) ;
      rms_err = ems/(2*rms_val) ;
      //////////// double ms = mssum ;
      //////////// double rms_check = sqrt( ms ) ;
      //////////// printf(" Check on RMS: compare %.3f with %.3f\n", rv, rms_check ) ;

      printf(" Zero mean RMS is %.3f +/- %.3f\n", rms_val, rms_err ) ;

      return ;

   }

