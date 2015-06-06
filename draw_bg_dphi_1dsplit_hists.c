
#include "TFile.h"
#include "TH1F.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "THStack.h"
#include "TSystem.h"
#include "TLine.h"

   TFile* p_hist_file ;
   float  mdp_cut ;

   float R_qcd_val ;
   float R_qcd_err ;
   char var_name[100] ;

  //---------
   void add_lostlep_samples( const char* hnamebase ) ;
   void add_nb_hists( const char* hname ) ;
   void add_njet_hists( const char* hname ) ;
   TH3* get_hist_3d( const char* hname ) ;
   TH1F* get_hist_1d( const char* hname ) ;
   void draw_stack( const char* hname_base, const char* hname_end, float xmin=0., float xmax=-1. ) ;
   void add_overflow_to_last_bin( TH1* hp ) ;
   float get_of_frac( TH1* hp, int xblow, int xbhigh ) ;
  //---------


   void draw_bg_dphi_1dsplit_hists( const char* arg_var_name = "mdpn",  float arg_mdp_cut = 6.0, const char* infile = "outputfiles/fill-bg-dphi-hists2.root" ) {

      gStyle -> SetOptTitle(0) ;
      sprintf( var_name, "%s", arg_var_name ) ;
      mdp_cut = arg_mdp_cut ;

      gDirectory -> cd("Rint:/") ;
      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_draw_bg_dphi_1dsplit_hists" ) ;
      if ( can == 0x0 ) {
         can = new TCanvas( "can_draw_bg_dphi_1dsplit_hists", "QCD 1D split", 700, 500 ) ;
      }
      can -> Clear() ;
      can -> SaveAs( "outputfiles/split1d-dummy.pdf" ) ;

      gDirectory -> Delete( "h*") ;
      char dirname[1000] ;
      sprintf( dirname, "%s:/", infile ) ;
      bool cdstat(false) ;
      cdstat = gDirectory -> cd( dirname ) ;
      if ( !cdstat ) {
         printf("\n\n === Opening file : %s\n\n", infile ) ;
         p_hist_file = new TFile( infile, "READ" ) ;
         gDirectory -> pwd() ;
      } else {
         printf("\n\n File %s already open.\n\n", infile ) ;
         p_hist_file = gFile ;
         if ( p_hist_file == 0x0 ) { printf("\n\n *** null file pointer.\n\n" ) ; return; }
      }



      int nsamples(3) ;
      char sname[3][20] ;
      int scolor[100] ;
      sprintf( sname[0], "znunu" ) ;    scolor[0] = kOrange - 3 ;
      sprintf( sname[1], "lostlep" ) ;  scolor[1] = kBlue - 4 ;
      sprintf( sname[2], "qcd" ) ;      scolor[2] = kRed + 3 ;


      printf("\n\n Adding lostlep samples.\n\n") ;
      for ( int nji=1; nji<=5; nji++ ) {
         for ( int nbi=0; nbi<=3; nbi++ ) {
            char hnamebase[1000] ;
            sprintf( hnamebase, "h_%svsmhtvsht_nb%d_nj%d_zl", var_name, nbi, nji ) ;
            add_lostlep_samples( hnamebase ) ;
         } // nbi
      } // nji
      printf("\n\n") ;

      for ( int si=0; si<nsamples; si++ ) {

         char hname[1000] ;
         for ( int nji=1; nji<=5; nji++ ) {
            sprintf( hname, "h_%svsmhtvsht_nb0_nj%d_zl_%s", var_name, nji, sname[si] ) ;
            add_nb_hists( hname ) ;
         }

         sprintf( hname, "h_%svsmhtvsht_nbsum_nj1_zl_%s", var_name, sname[si] ) ;
         add_njet_hists( hname ) ;

      } // si


      float htlow[4] ;
      float hthigh[4] ;
      htlow[1] = 500. ; hthigh[1] = 800. ;
      htlow[2] = 800. ; hthigh[2] = 1200. ;
      htlow[3] =1200. ; hthigh[3] = 20000. ;


      for ( int hbi=1; hbi<=3; hbi++ ) {

         for ( int si=0; si<nsamples; si++ ) {

            char hname3d[1000] ;
            char hnameproj[1000] ;

            sprintf( hname3d, "h_%svsmhtvsht_nbsum_njsum_zl_%s", var_name, sname[si] ) ;
            sprintf( hnameproj, "h_%s_nbsum_njsum_zl_%s_ht%d", var_name, sname[si], hbi ) ;

            TH3* hp3d = get_hist_3d( hname3d ) ;
            TH1D* hpproj = hp3d -> ProjectionZ( hnameproj, hbi+1, hbi+1, 2, 5 ) ;
            add_overflow_to_last_bin( hpproj ) ;

         } // si

         char hname_base[1000] ;
         char hname_end[1000] ;
         sprintf( hname_base, "h_%s_nbsum_njsum_zl", var_name ) ;
         sprintf( hname_end, "ht%d", hbi ) ;
         draw_stack( hname_base, hname_end ) ;

         TText* tlabel = new TText() ;
         char label[1000] ;
         sprintf( label, "HT%d [%.0f, %.0f], sum over MHT, Njet, and Nb bins", hbi, htlow[hbi], hthigh[hbi] ) ;
         tlabel -> DrawTextNDC( 0.10, 0.92, label ) ;

         char fname[10000] ;
         sprintf( fname, "outputfiles/split1d-ht%d-%s.pdf", hbi, var_name ) ;
         can -> SaveAs( fname ) ;

      } // hbi



      float R1_qcd_val ;
      float R1_qcd_err ;

      for ( int mbi=1; mbi<=4; mbi++ ) {

         for ( int si=0; si<nsamples; si++ ) {

            char hname3d[1000] ;
            char hnameproj[1000] ;

            sprintf( hname3d, "h_%svsmhtvsht_nbsum_njsum_zl_%s", var_name, sname[si] ) ;
            sprintf( hnameproj, "h_%s_nbsum_njsum_zl_%s_mht%d", var_name, sname[si], mbi ) ;

            TH3* hp3d = get_hist_3d( hname3d ) ;
            TH1D* hpproj = hp3d -> ProjectionZ( hnameproj, 2, 4, mbi+1, mbi+1 ) ;
            add_overflow_to_last_bin( hpproj ) ;

         } // si

         char hname_base[1000] ;
         char hname_end[1000] ;
         sprintf( hname_base, "h_%s_nbsum_njsum_zl", var_name ) ;
         sprintf( hname_end, "mht%d", mbi ) ;
         draw_stack( hname_base, hname_end ) ;

         TText* tlabel = new TText() ;
         char label[1000] ;
         char mhtbinstring[10] ;
         if ( mbi == 1 ) sprintf( mhtbinstring, "MHT1a [200,300]" ) ;
         if ( mbi == 2 ) sprintf( mhtbinstring, "MHT1b [300,500]" ) ;
         if ( mbi == 3 ) sprintf( mhtbinstring, "MHT2 [500,750]" ) ;
         if ( mbi == 4 ) sprintf( mhtbinstring, "MHT3 [750,+]" ) ;
         sprintf( label, "%s, sum over HT, Njet, and Nb bins", mhtbinstring ) ;
         tlabel -> DrawTextNDC( 0.10, 0.92, label ) ;

         if ( mbi==1 ) {
            R1_qcd_val = R_qcd_val ;
            R1_qcd_err = R_qcd_err ;
         } else {
            float rn_over_r1_val(0.) ;
            float rn_over_r1_err(0.) ;
            if ( R1_qcd_val > 0 && R_qcd_val > 0 ) {
               rn_over_r1_val = R_qcd_val / R1_qcd_val ;
               rn_over_r1_err = rn_over_r1_val * sqrt( pow( R_qcd_err/R_qcd_val, 2 ) + pow( R1_qcd_err/R1_qcd_val, 2 ) ) ;
            }
            TText* tt = new TText() ;
            tt -> SetTextSize( 0.045 ) ;
            char label[1000] ;
            sprintf( label, "R%d / R1 = %.3f +/- %.3f\n", mbi, rn_over_r1_val, rn_over_r1_err ) ;
            tt -> DrawTextNDC( 0.20, 0.62, label ) ;
         }

         char fname[10000] ;
         sprintf( fname, "outputfiles/split1d-mht%d-%s.pdf", mbi, var_name ) ;
         can -> SaveAs( fname ) ;


      } // mbi


      for ( int nji=1; nji<=5; nji++ ) {

         for ( int si=0; si<nsamples; si++ ) {

            char hname3d[1000] ;
            char hnameproj[1000] ;

            sprintf( hname3d, "h_%svsmhtvsht_nbsum_nj%d_zl_%s", var_name, nji, sname[si] ) ;
            sprintf( hnameproj, "h_%s_nbsum_nj%d_zl_%s", var_name, nji, sname[si] ) ;

            TH3* hp3d = get_hist_3d( hname3d ) ;
            TH1D* hpproj = hp3d -> ProjectionZ( hnameproj, 2, 4, 2, 5 ) ;
            add_overflow_to_last_bin( hpproj ) ;

         } // si

         char hname_base[1000] ;
         char hname_end[1000] ;
         sprintf( hname_base, "h_%s_nbsum_nj%d_zl", var_name, nji ) ;
         sprintf( hname_end, "" ) ;
         draw_stack( hname_base, hname_end ) ;

         TText* tlabel = new TText() ;
         char label[1000] ;
         char njbinstring[10] ;
         if ( nji == 1 ) sprintf( njbinstring, "Njet1a [4]" ) ;
         if ( nji == 2 ) sprintf( njbinstring, "Njet1b [5]" ) ;
         if ( nji == 3 ) sprintf( njbinstring, "Njet1c [6]" ) ;
         if ( nji == 4 ) sprintf( njbinstring, "Njet2 [7,8]" ) ;
         if ( nji == 5 ) sprintf( njbinstring, "Njet3 [9,+]" ) ;
         sprintf( label, "%s, sum over HT, MHT, and Nb bins", njbinstring ) ;
         tlabel -> DrawTextNDC( 0.10, 0.92, label ) ;

         if ( nji==1 ) {
            R1_qcd_val = R_qcd_val ;
            R1_qcd_err = R_qcd_err ;
         } else {
            float rn_over_r1_val(0.) ;
            float rn_over_r1_err(0.) ;
            if ( R1_qcd_val > 0 && R_qcd_val > 0 ) {
               rn_over_r1_val = R_qcd_val / R1_qcd_val ;
               rn_over_r1_err = rn_over_r1_val * sqrt( pow( R_qcd_err/R_qcd_val, 2 ) + pow( R1_qcd_err/R1_qcd_val, 2 ) ) ;
            }
            TText* tt = new TText() ;
            tt -> SetTextSize( 0.045 ) ;
            char label[1000] ;
            sprintf( label, "R%d / R1 = %.3f +/- %.3f\n", nji, rn_over_r1_val, rn_over_r1_err ) ;
            tt -> DrawTextNDC( 0.20, 0.62, label ) ;
         }

         char fname[10000] ;
         sprintf( fname, "outputfiles/split1d-njet%d-%s.pdf", nji, var_name ) ;
         can -> SaveAs( fname ) ;

      } // nji





   } // draw_bg_dphi_1dsplit_hists


  //===============================================================

   TH3* get_hist_3d( const char* hname ) {

      TH3* hp(0x0) ;
      hp = (TH3*) gDirectory -> Get( hname ) ;
      if ( hp == 0x0 ) {
         hp = (TH3*) p_hist_file -> Get( hname ) ;
         if ( hp == 0x0 ) { printf("\n\n **** get_hist_3d: Missing histogram: %s\n\n", hname ) ; gSystem->Exit(-1) ; }
         printf("  get_hist_3d: found %s in file.\n", hname ) ;
      } else {
         printf("  get_hist_3d: found %s in memory.\n", hname ) ;
      }

      return hp ;

   } // get_hist_3d

  //===============================================================

   TH1F* get_hist_1d( const char* hname ) {

      TH1F* hp(0x0) ;
      hp = (TH1F*) gDirectory -> Get( hname ) ;
      if ( hp == 0x0 ) {
         hp = (TH1F*) p_hist_file -> Get( hname ) ;
         if ( hp == 0x0 ) { printf("\n\n **** get_hist_1d: Missing histogram: %s\n\n", hname ) ; gSystem->Exit(-1) ; }
         printf("  get_hist_1d: found %s in file.\n", hname ) ;
      } else {
         printf("  get_hist_1d: found %s in memory.\n", hname ) ;
      }

      return hp ;

   } // get_hist_1d

  //===============================================================


   void draw_stack( const char* hname_base, const char* hname_end, float xmin, float xmax ) {

       char hpname[1000] ;

       char xaxislabel[100] ;
       if ( strcmp( var_name, "mdpn" ) == 0 ) {
          sprintf( xaxislabel, "minDeltaPhiN" ) ;
       } else {
          sprintf( xaxislabel, "minDeltaPhi" ) ;
       }

       if ( strlen( hname_end ) > 0 ) {
          sprintf( hpname, "%s_znunu_%s", hname_base, hname_end ) ;
       } else {
          sprintf( hpname, "%s_znunu", hname_base ) ;
       }
       TH1F* hp_znunu = get_hist_1d( hpname ) ;
       hp_znunu -> SetFillColor( kOrange - 5 ) ;

       if ( strlen( hname_end ) > 0 ) {
          sprintf( hpname, "%s_lostlep_%s", hname_base, hname_end ) ;
       } else {
          sprintf( hpname, "%s_lostlep", hname_base ) ;
       }
       TH1F* hp_lostlep = get_hist_1d( hpname ) ;
       hp_lostlep -> SetFillColor( kBlue - 6 ) ;

       if ( strlen( hname_end ) > 0 ) {
          sprintf( hpname, "%s_qcd_%s", hname_base, hname_end ) ;
       } else {
          sprintf( hpname, "%s_qcd", hname_base ) ;
       }
       TH1F* hp_qcd = get_hist_1d( hpname ) ;
       hp_qcd -> SetFillColor( kRed + 1 ) ;

       sprintf( hpname, "%s_stack_%s", hname_base, hname_end ) ;

       THStack* hs = new THStack( hpname, hpname ) ;
       hs -> Add( hp_znunu ) ;
       hs -> Add( hp_lostlep ) ;
       hs -> Add( hp_qcd ) ;

       hs -> Draw("hist") ;
       if ( xmax > xmin ) {
          TH1* hp = hs -> GetHistogram() ;
          hp -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
          hs -> Draw("hist") ;
          int iblow  = hp_znunu -> FindBin( xmax+1e-5 ) ;
          int ibhigh = hp_znunu -> GetNbinsX() ;
          float qcd_of_frac = get_of_frac( hp_qcd, iblow, ibhigh )  ;
          float lostlep_of_frac = get_of_frac( hp_lostlep, iblow, ibhigh )  ;
          float znunu_of_frac = get_of_frac( hp_znunu, iblow, ibhigh )  ;
          printf(" %s : overflow fracs:  QCD = %.3f ,   lostlep = %.3f ,   Znunu = %.3f\n",
             hpname, qcd_of_frac, lostlep_of_frac, znunu_of_frac ) ;
          TText* tt = new TText() ;
          tt -> SetTextAlign(33) ;
          char label[100] ;
          float tx = 0.84 ;
          tt -> DrawTextNDC( tx, 0.75, "Overflow frac" ) ;
          sprintf( label, "QCD %.2f", qcd_of_frac ) ;
          tt -> DrawTextNDC( tx, 0.70, label ) ;
          sprintf( label, "LL %.2f", lostlep_of_frac ) ;
          tt -> DrawTextNDC( tx, 0.65, label ) ;
          sprintf( label, "Zinv %.2f", znunu_of_frac ) ;
          tt -> DrawTextNDC( tx, 0.60, label ) ;
       }
       hs -> Draw("same") ;

       TLine* l = new TLine() ;
       l -> SetLineStyle(2) ;
       l -> DrawLine( mdp_cut, 0., mdp_cut, hs->GetMaximum() ) ;

       int cutbin = hp_qcd -> FindBin( mdp_cut+1e-5 ) ;
       int nbins = hp_qcd -> GetNbinsX() ;

       double qcd_low_val(0.) ;
       double qcd_low_err(0.) ;
       double qcd_high_val(0.) ;
       double qcd_high_err(0.) ;

       qcd_low_val  = hp_qcd -> IntegralAndError( 1, cutbin, qcd_low_err ) ;
       qcd_high_val = hp_qcd -> IntegralAndError( cutbin, nbins, qcd_high_err ) ;

       float qcd_ratio_val(0.) ;
       float qcd_ratio_err(0.) ;
       if ( qcd_low_val > 0 ) {
          qcd_ratio_val = qcd_high_val / qcd_low_val ;
          if ( qcd_high_val > 0 ) {
             qcd_ratio_err = qcd_ratio_val * sqrt( pow( qcd_low_err / qcd_low_val , 2 ) + pow( qcd_high_err / qcd_high_val , 2 ) ) ;
          }
       }

       TText* tt = new TText() ;
       tt -> SetTextSize( 0.045 ) ;
       char label[1000] ;

       sprintf( label, "Rqcd = (%.1f +/- %.1f ) / (%.1f +/- %.1f)", qcd_high_val, qcd_high_err, qcd_low_val, qcd_low_err ) ;
       tt -> DrawTextNDC( 0.20, 0.80, label ) ;

       sprintf( label, "Rqcd = %.3f +/- %.3f ", qcd_ratio_val, qcd_ratio_err ) ;
       tt -> DrawTextNDC( 0.20, 0.74, label ) ;

       sprintf( label, "QCD MC stat error only" ) ;
       tt -> DrawTextNDC( 0.20, 0.68, label ) ;

       R_qcd_val = qcd_ratio_val ;
       R_qcd_err = qcd_ratio_err ;

       tt -> SetTextAlign( 31 ) ;
       tt -> DrawTextNDC( 0.90, 0.01, xaxislabel ) ;


   } // draw_stack

  //===============================================================

   float get_of_frac( TH1* hp, int xblow, int xbhigh ) {
      float of  = hp -> Integral( xblow, xbhigh ) ;
      float all = hp -> Integral( 1, xbhigh ) ;
      float of_frac = 0. ;
      if ( all > 0 ) of_frac = of / all ;
      return of_frac ;
   } // get_of_frac

  //===============================================================

   void add_nb_hists( const char* hname ) {

      TH3* hp = get_hist_3d( hname ) ;

      TString hnameout = hp -> GetName() ;
      hnameout.ReplaceAll( "nb0", "nbsum" ) ;

      TString htitleout = hp -> GetTitle() ;
      htitleout.ReplaceAll( "nb0", "nbsum" ) ;

      TH3* hout = (TH3*) hp -> Clone( hnameout ) ;
      hout -> SetTitle( htitleout ) ;

      TString otherhname = hp -> GetName() ;
      TH3* otherhp ;

      otherhname.ReplaceAll( "nb0", "nb1" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;

      otherhname.ReplaceAll( "nb1", "nb2" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;

      otherhname.ReplaceAll( "nb2", "nb3" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;


   } // add_nb_hists

  //===============================================================

   void add_njet_hists( const char* hname ) {

      TH3* hp = get_hist_3d( hname ) ;

      TString hnameout = hp -> GetName() ;
      hnameout.ReplaceAll( "nj1", "njsum" ) ;

      TString htitleout = hp -> GetTitle() ;
      htitleout.ReplaceAll( "nj1", "njsum" ) ;

      TH3* hout = (TH3*) hp -> Clone( hnameout ) ;
      hout -> SetTitle( htitleout ) ;

      TString otherhname = hp -> GetName() ;
      TH3* otherhp ;

      otherhname.ReplaceAll( "nj1", "nj2" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;

      otherhname.ReplaceAll( "nj2", "nj3" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;

      otherhname.ReplaceAll( "nj3", "nj4" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;

      otherhname.ReplaceAll( "nj4", "nj5" ) ;
      otherhp = get_hist_3d( otherhname ) ;
      hout -> Add( otherhp ) ;


   } // add_njet_hists

  //===============================================================

   void add_lostlep_samples( const char* hnamebase ) {

       char hname[1000] ;

       TH3* h_ttbar ;
       TH3* h_wjets ;
       TH3* h_sngltop ;

       sprintf( hname, "%s_ttbar", hnamebase ) ;
       h_ttbar = get_hist_3d( hname ) ;

       sprintf( hname, "%s_wjets", hnamebase ) ;
       h_wjets = get_hist_3d( hname ) ;

       sprintf( hname, "%s_sngltop", hnamebase ) ;
       h_sngltop = get_hist_3d( hname ) ;

       sprintf( hname, "%s_lostlep", hnamebase ) ;
       TH3F* h_ll = (TH3F*) h_ttbar -> Clone( hname ) ;
       TString newtitle = h_ll -> GetTitle() ;
       newtitle.ReplaceAll( "ttbar", "lostlep" ) ;
       h_ll -> SetTitle( newtitle ) ;
       h_ll -> Add( h_wjets ) ;
       h_ll -> Add( h_sngltop ) ;

   } // add_lostlep_samples

  //===============================================================


   void add_overflow_to_last_bin( TH1* hp ) {

      if ( hp == 0x0 ) return ;

      int nb = hp -> GetNbinsX() ;

      double lbval = hp -> GetBinContent( nb ) ;
      double lberr = hp -> GetBinError( nb ) ;

      double ofval = hp -> GetBinContent( nb+1 ) ;
      double oferr = hp -> GetBinError( nb+1 ) ;

      double newval = lbval + ofval ;
      double newerr = sqrt( lberr*lberr + oferr*oferr ) ;

      hp -> SetBinContent( nb, newval ) ;
      hp -> SetBinError( nb, newerr ) ;

   } // add_overflow_to_last_bin


  //===============================================================




