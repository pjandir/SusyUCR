

#include "TSystem.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TText.h"
#include "THStack.h"
#include "TLine.h"
#include "TStyle.h"
#include "TPad.h"

#include "histio.c"


   //float mdp_cut[4] = { 0.35, 0.25, 0.20, 0.20 } ;
   ///float mdp_cut[4] = { 0.30, 0.30, 0.30, 0.30 } ;

   float mdp_cut ;

   TFile* p_hist_file ;

   TH1F* h_dphi_zl_qcd ;
   TH1F* h_dphi_zl_lostlep ;
   TH1F* h_dphi_zl_znunu ;
   TH1F* h_dphi_zl_sig ;

   TH1F* h_dphi_sl_qcd ;
   TH1F* h_dphi_sl_lostlep ;
   TH1F* h_dphi_sl_znunu ;
   TH1F* h_dphi_sl_sig ;

   TPad* text_pad ;
     char text_line[1000] ;
   TText* ttext ;

     float text_y_ndc(0.95) ;
     float text_x_ndc(0.0) ;
     float text_size(0.015) ;
     float text_dy_ndc(0.017) ;

   int sb_mht_ht_ind ;
   int sb_nj_ind ;

   int fb_mht_ind[10] ;
   int fb_ht_ind[10] ;
   int fb_nj_ind[10] ;
   int fb_nb_ind ;
   int n_fb_in_sb ;

   double n_zl_qcd_val[10] ;
   double n_zl_qcd_err[10] ;
   double n_zl_lostlep_val[10] ;
   double n_zl_lostlep_err[10] ;
   double n_zl_znunu_val[10] ;
   double n_zl_znunu_err[10] ;
   double n_zl_sig_val[10] ;
   double n_zl_sig_err[10] ;
   double n_zl_val[10] ;

   double n_ldp_qcd_val[10] ;
   double n_ldp_qcd_err[10] ;
   double n_ldp_lostlep_val[10] ;
   double n_ldp_lostlep_err[10] ;
   double n_ldp_znunu_val[10] ;
   double n_ldp_znunu_err[10] ;
   double n_ldp_val[10] ;

   double n_slldp_qcd_val[10] ;
   double n_slldp_qcd_err[10] ;
   double n_slldp_lostlep_val[10] ;
   double n_slldp_lostlep_err[10] ;
   double n_slldp_znunu_val[10] ;
   double n_slldp_znunu_err[10] ;
   double n_slldp_val[10] ;

   double n_sl_qcd_val[10] ;
   double n_sl_qcd_err[10] ;
   double n_sl_lostlep_val[10] ;
   double n_sl_lostlep_err[10] ;
   double n_sl_znunu_val[10] ;
   double n_sl_znunu_err[10] ;
   double n_sl_val[10] ;

   char sb_name[100] ;

   char signal_name[100] ;

   char varname[100] ;

  //---------
   void print_usage() ;
   TH1F* get_hist( const char* hname ) ;
   int hist_bin( int mht_bin, int ht_bin ) ;
   float get_of_frac( TH1* hp, int xblow, int xbhigh ) ;
   void fill_fb_index_arrays( ) ;
   void collect_dphi_hists() ;
   void calc_qcd() ;
   float calc_qcd_err() ;
   void fill_kqcd_arrays() ;
  //---------

   float Kqcd_ht_val[3] ;
   float Kqcd_ht_err[3] ;

   float Kqcd_mht_val[4] ;
   float Kqcd_mht_err[4] ;

   float Kqcd_nj_val[5] ;
   float Kqcd_nj_err[5] ;


  //===============================================================

   void draw_single_searchbin( int arg_sb_mht_ht_ind=-1, int arg_sb_nj_ind=1, int sb_nb_ind=0,
                               float arg_mdp_cut = 6.0,
                               const char* arg_varname = "mdpn",
                               float xmin = 0., float xmax = -1.,
                               const char* infile =  "outputfiles/fill-bg-dphi-hists2-t1bbbbH-postdraw.root",
                               const char* signame = "t1bbbbH" ) {

      sprintf( varname, "%s", arg_varname ) ;

      mdp_cut = arg_mdp_cut ;
      text_y_ndc = 0.95 ;

      sb_mht_ht_ind = arg_sb_mht_ht_ind ;
      sb_nj_ind = arg_sb_nj_ind ;

      fill_kqcd_arrays() ;

      if ( sb_mht_ht_ind < 1 ) {
         printf("\n\n") ;
         printf("  draw_single_searchbin( sb_mht_ht_ind, sb_nj_ind, sb_nb_ind, xmin, xmax, infile, signame )\n") ;
         printf("\n\n") ;
         return ;
      }

      sprintf( signal_name, "%s", signame ) ;

      fb_nb_ind = sb_nb_ind ;
      fill_fb_index_arrays( ) ;
      if ( n_fb_in_sb <= 0 ) { printf("\n\n *** bad sb indices.\n\n") ; return ; }




      gStyle -> SetOptStat(0) ;

      ///////gDirectory -> pwd() ;
      gDirectory -> cd("Rint:/") ;
      ///////gDirectory -> pwd() ;
      TCanvas* can = (TCanvas*) gDirectory -> FindObject( "can_draw_single_searchbin" ) ;
      if ( can == 0x0 ) {
         can = new TCanvas( "can_draw_single_searchbin", "Single search bin", 1900, 1250 ) ;
      }


      gDirectory -> Delete( "h*") ;
      p_hist_file = new TFile( infile, "READ" ) ;
      ///////gDirectory -> pwd() ;

      char hname[1000] ;

      collect_dphi_hists() ;



   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     can -> Clear() ;
     can -> Divide(3,4) ;
     can -> cd(1) ;



     int nbins = h_dphi_zl_qcd -> GetNbinsX() ;
     int cutbin = h_dphi_zl_qcd -> FindBin( mdp_cut+1e-5 ) ;

     double sb_n_zl_qcd_err(0.) ;
     float sb_n_zl_qcd_val = h_dphi_zl_qcd -> IntegralAndError( cutbin, nbins, sb_n_zl_qcd_err ) ;
     double sb_n_zl_lostlep_err(0.) ;
     float sb_n_zl_lostlep_val = h_dphi_zl_lostlep -> IntegralAndError( cutbin, nbins, sb_n_zl_lostlep_err ) ;
     double sb_n_zl_znunu_err(0.) ;
     float sb_n_zl_znunu_val = h_dphi_zl_znunu -> IntegralAndError( cutbin, nbins, sb_n_zl_znunu_err ) ;

     float sb_n_zl_val = sb_n_zl_qcd_val + sb_n_zl_lostlep_val + sb_n_zl_znunu_val ;

     double sb_n_ldp_qcd_err(0.) ;
     float sb_n_ldp_qcd_val = h_dphi_zl_qcd -> IntegralAndError( 1, cutbin-1, sb_n_ldp_qcd_err ) ;
     double sb_n_ldp_lostlep_err(0.) ;
     float sb_n_ldp_lostlep_val = h_dphi_zl_lostlep -> IntegralAndError( 1, cutbin-1, sb_n_ldp_lostlep_err ) ;
     double sb_n_ldp_znunu_err(0.) ;
     float sb_n_ldp_znunu_val = h_dphi_zl_znunu -> IntegralAndError( 1, cutbin-1, sb_n_ldp_znunu_err ) ;

     float sb_n_ldp_val = sb_n_ldp_qcd_val + sb_n_ldp_lostlep_val + sb_n_ldp_znunu_val ;


     double sb_n_sl_qcd_err(0.) ;
     float sb_n_sl_qcd_val = h_dphi_sl_qcd -> IntegralAndError( cutbin, nbins, sb_n_sl_qcd_err ) ;
     double sb_n_sl_lostlep_err(0.) ;
     float sb_n_sl_lostlep_val = h_dphi_sl_lostlep -> IntegralAndError( cutbin, nbins, sb_n_sl_lostlep_err ) ;
     double sb_n_sl_znunu_err(0.) ;
     float sb_n_sl_znunu_val = h_dphi_sl_znunu -> IntegralAndError( cutbin, nbins, sb_n_sl_znunu_err ) ;

     double sb_n_slldp_qcd_err(0.) ;
     float sb_n_slldp_qcd_val = h_dphi_sl_qcd -> IntegralAndError( 1, cutbin-1, sb_n_slldp_qcd_err ) ;
     double sb_n_slldp_lostlep_err(0.) ;
     float sb_n_slldp_lostlep_val = h_dphi_sl_lostlep -> IntegralAndError( 1, cutbin-1, sb_n_slldp_lostlep_err ) ;
     double sb_n_slldp_znunu_err(0.) ;
     float sb_n_slldp_znunu_val = h_dphi_sl_znunu -> IntegralAndError( 1, cutbin-1, sb_n_slldp_znunu_err ) ;








     TH1F* h_bincounts_qcd     = new TH1F( "h_bincounts_qcd"    , "Bin counts QCD"    , 9, 0.5, 9.5 ) ;
     TH1F* h_bincounts_lostlep = new TH1F( "h_bincounts_lostlep", "Bin counts lostlep", 9, 0.5, 9.5 ) ;
     TH1F* h_bincounts_znunu   = new TH1F( "h_bincounts_znunu"  , "Bin counts znunu"  , 9, 0.5, 9.5 ) ;

     h_bincounts_qcd -> SetFillColor( h_dphi_zl_qcd -> GetFillColor() ) ;
     h_bincounts_lostlep -> SetFillColor( h_dphi_zl_lostlep -> GetFillColor() ) ;
     h_bincounts_znunu -> SetFillColor( h_dphi_zl_znunu -> GetFillColor() ) ;


     h_bincounts_qcd -> SetBinContent( 2, sb_n_sl_qcd_val )    ; h_bincounts_qcd -> SetBinError( 2, sb_n_sl_qcd_err ) ;    h_bincounts_qcd -> GetXaxis() -> SetBinLabel( 2, "SL" ) ;
     h_bincounts_qcd -> SetBinContent( 4, sb_n_slldp_qcd_val ) ; h_bincounts_qcd -> SetBinError( 4, sb_n_slldp_qcd_err ) ; h_bincounts_qcd -> GetXaxis() -> SetBinLabel( 4, "SLLDP" ) ;
     h_bincounts_qcd -> SetBinContent( 6, sb_n_ldp_qcd_val )   ; h_bincounts_qcd -> SetBinError( 6, sb_n_ldp_qcd_err ) ;   h_bincounts_qcd -> GetXaxis() -> SetBinLabel( 6, "LDP" ) ;
     h_bincounts_qcd -> SetBinContent( 8, sb_n_zl_qcd_val )    ; h_bincounts_qcd -> SetBinError( 8, sb_n_zl_qcd_err ) ;    h_bincounts_qcd -> GetXaxis() -> SetBinLabel( 8, "ZL" ) ;

     h_bincounts_lostlep -> SetBinContent( 2, sb_n_sl_lostlep_val )    ; h_bincounts_lostlep -> SetBinError( 2, sb_n_sl_lostlep_err ) ;     h_bincounts_lostlep -> GetXaxis() -> SetBinLabel( 2, "SL" ) ;
     h_bincounts_lostlep -> SetBinContent( 4, sb_n_slldp_lostlep_val ) ; h_bincounts_lostlep -> SetBinError( 4, sb_n_slldp_lostlep_err ) ;  h_bincounts_lostlep -> GetXaxis() -> SetBinLabel( 4, "SLLDP" ) ;
     h_bincounts_lostlep -> SetBinContent( 6, sb_n_ldp_lostlep_val )   ; h_bincounts_lostlep -> SetBinError( 6, sb_n_ldp_lostlep_err ) ;    h_bincounts_lostlep -> GetXaxis() -> SetBinLabel( 6, "LDP" ) ;
     h_bincounts_lostlep -> SetBinContent( 8, sb_n_zl_lostlep_val )    ; h_bincounts_lostlep -> SetBinError( 8, sb_n_zl_lostlep_err ) ;     h_bincounts_lostlep -> GetXaxis() -> SetBinLabel( 8, "ZL" ) ;

     h_bincounts_znunu -> SetBinContent( 2, sb_n_sl_znunu_val )    ; h_bincounts_znunu -> SetBinError( 2, sb_n_sl_znunu_err ) ;      h_bincounts_znunu -> GetXaxis() -> SetBinLabel( 2, "SL" ) ;
     h_bincounts_znunu -> SetBinContent( 4, sb_n_slldp_znunu_val ) ; h_bincounts_znunu -> SetBinError( 4, sb_n_slldp_znunu_err ) ;   h_bincounts_znunu -> GetXaxis() -> SetBinLabel( 4, "SLLDP" ) ;
     h_bincounts_znunu -> SetBinContent( 6, sb_n_ldp_znunu_val )   ; h_bincounts_znunu -> SetBinError( 6, sb_n_ldp_znunu_err ) ;     h_bincounts_znunu -> GetXaxis() -> SetBinLabel( 6, "LDP" ) ;
     h_bincounts_znunu -> SetBinContent( 8, sb_n_zl_znunu_val )    ; h_bincounts_znunu -> SetBinError( 8, sb_n_zl_znunu_err ) ;      h_bincounts_znunu -> GetXaxis() -> SetBinLabel( 8, "ZL" ) ;

     // h_bincounts_qcd -> GetXaxis() -> LabelsOption( "v" ) ;
     // h_bincounts_lostlep -> GetXaxis() -> LabelsOption( "v" ) ;
     // h_bincounts_znunu -> GetXaxis() -> LabelsOption( "v" ) ;

     THStack* h_bincounts_stack = new THStack( "h_bincounts_stack", "Bincounts" ) ;
     h_bincounts_stack -> Add( h_bincounts_znunu ) ;
     h_bincounts_stack -> Add( h_bincounts_lostlep ) ;
     h_bincounts_stack -> Add( h_bincounts_qcd ) ;

     TH1F* h_bincounts_sum = (TH1F*) h_bincounts_znunu -> Clone( "h_bincounts_sum" ) ;
     h_bincounts_sum -> Add( h_bincounts_lostlep ) ;
     h_bincounts_sum -> Add( h_bincounts_qcd ) ;
     h_bincounts_sum -> SetLabelSize( 0.06, "x" ) ;


     can -> cd( 10 ) ;
     h_bincounts_sum -> Draw( "hist" ) ;
     h_bincounts_stack -> Draw( "hist same" ) ;
     h_bincounts_stack -> Draw( "same" ) ;



     TLine* line = new TLine() ;
     line -> SetLineStyle(2) ;


     sprintf( hname, "h_mdpvsmhtvsht_%s_zl_stack", sb_name ) ;
     THStack* hs_zl = new THStack( hname, hname ) ;
     hs_zl -> Add( h_dphi_zl_znunu ) ;
     hs_zl -> Add( h_dphi_zl_lostlep ) ;
     hs_zl -> Add( h_dphi_zl_qcd ) ;

     sprintf( hname, "h_mdpvsmhtvsht_%s_sl_stack", sb_name ) ;
     THStack* hs_sl = new THStack( hname, hname ) ;
     hs_sl -> Add( h_dphi_sl_znunu ) ;
     hs_sl -> Add( h_dphi_sl_lostlep ) ;
     hs_sl -> Add( h_dphi_sl_qcd ) ;

     can -> cd(1) ;
     hs_zl -> Draw("hist") ;
     if ( xmax > xmin ) {
        TH1* hp = hs_zl -> GetHistogram() ;
        hp -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
        hs_zl -> Draw("hist") ;
        int iblow  = h_dphi_zl_znunu -> FindBin( xmax+1e-5 ) ;
        int ibhigh = h_dphi_zl_znunu -> GetNbinsX() ;
        float qcd_of_frac = get_of_frac( h_dphi_zl_qcd, iblow, ibhigh )  ;
        float lostlep_of_frac = get_of_frac( h_dphi_zl_lostlep, iblow, ibhigh )  ;
        float znunu_of_frac = get_of_frac( h_dphi_zl_znunu, iblow, ibhigh )  ;
        printf(" %s : overflow fracs:  QCD = %.3f ,   lostlep = %.3f ,   Znunu = %.3f\n",
           hname, qcd_of_frac, lostlep_of_frac, znunu_of_frac ) ;
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
     hs_zl -> Draw("same") ;
     line -> DrawLine( mdp_cut, 0., mdp_cut, hs_zl -> GetMaximum() ) ;



     can -> cd(4) ;
     hs_sl -> Draw("hist") ;
     if ( xmax > xmin ) {
        TH1* hp = hs_sl -> GetHistogram() ;
        hp -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
        hs_sl -> Draw("hist") ;
        int iblow  = h_dphi_sl_znunu -> FindBin( xmax+1e-5 ) ;
        int ibhigh = h_dphi_sl_znunu -> GetNbinsX() ;
        float qcd_of_frac = get_of_frac( h_dphi_sl_qcd, iblow, ibhigh )  ;
        float lostlep_of_frac = get_of_frac( h_dphi_sl_lostlep, iblow, ibhigh )  ;
        float znunu_of_frac = get_of_frac( h_dphi_sl_znunu, iblow, ibhigh )  ;
        printf(" %s : overflow fracs:  QCD = %.3f ,   lostlep = %.3f ,   Znunu = %.3f\n",
           hname, qcd_of_frac, lostlep_of_frac, znunu_of_frac ) ;
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
     hs_sl -> Draw("same") ;
     line -> DrawLine( mdp_cut, 0., mdp_cut, hs_sl -> GetMaximum() ) ;



     can -> cd(7) ;
     h_dphi_zl_sig -> Draw("hist") ;
     if ( xmax > xmin ) {
        h_dphi_zl_sig -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
        int iblow  = h_dphi_zl_sig -> FindBin( xmax+1e-5 ) ;
        int ibhigh = h_dphi_zl_sig -> GetNbinsX() ;
        float sig_of_frac = get_of_frac( h_dphi_zl_sig, iblow, ibhigh )  ;
        TText* tt = new TText() ;
        tt -> SetTextAlign(33) ;
        char label[100] ;
        float tx = 0.84 ;
        tt -> DrawTextNDC( tx, 0.75, "Overflow frac" ) ;
        sprintf( label, "%s %.2f", signame, sig_of_frac ) ;
        tt -> DrawTextNDC( tx, 0.70, label ) ;
     }
     line -> DrawLine( mdp_cut, 0., mdp_cut, h_dphi_zl_sig -> GetMaximum() ) ;


  //--------

     can -> cd(0) ;
     text_pad = new TPad( "text_pad", "text pad", 0.35, 0.05, 1.0, 1.0 ) ;
     text_pad -> Draw() ;
     text_pad -> cd() ;


     ttext = new TText() ;
     ttext -> SetTextSize( text_size ) ;
     ttext -> SetTextFont( 102 ) ;

     calc_qcd() ;

      printf("\n\n") ;
      text_y_ndc -= text_dy_ndc ;
      sprintf( text_line, "====== Total Event counts for %s", sb_name ) ;
      printf( "%s\n", text_line ) ;
      text_y_ndc -= text_dy_ndc ;
      ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;
      sprintf( text_line, "   lostlep : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f",
            sb_n_sl_lostlep_val, sb_n_sl_lostlep_err,
            sb_n_slldp_lostlep_val, sb_n_slldp_lostlep_err,
            sb_n_ldp_lostlep_val, sb_n_ldp_lostlep_err,
            sb_n_zl_lostlep_val, sb_n_zl_lostlep_err ) ;
      printf( "%s\n", text_line ) ;
      text_y_ndc -= text_dy_ndc ;
      ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

      sprintf( text_line, "   QCD     : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f",
            sb_n_sl_qcd_val, sb_n_sl_qcd_err,
            sb_n_slldp_qcd_val, sb_n_slldp_qcd_err,
            sb_n_ldp_qcd_val, sb_n_ldp_qcd_err,
            sb_n_zl_qcd_val, sb_n_zl_qcd_err ) ;
      printf( "%s\n", text_line ) ;
      text_y_ndc -= text_dy_ndc ;
      ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

      sprintf( text_line, "   Znunu   : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f",
            sb_n_sl_znunu_val, sb_n_sl_znunu_err,
            sb_n_slldp_znunu_val, sb_n_slldp_znunu_err,
            sb_n_ldp_znunu_val, sb_n_ldp_znunu_err,
            sb_n_zl_znunu_val, sb_n_zl_znunu_err ) ;
      printf( "%s\n", text_line ) ;
      text_y_ndc -= text_dy_ndc ;
      ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

      text_y_ndc -= text_dy_ndc ;


      char save_file[10000] ;
      sprintf( save_file, "outputfiles/single-searchbin-%s-%.2f-MHTHT%d-Nj%d-Nb%d.png", varname, mdp_cut, sb_mht_ht_ind, sb_nj_ind, sb_nb_ind ) ;
      can -> SaveAs( save_file ) ;



   } // draw_single_searchbin


  //===============================================================

   void print_usage() {
      printf("\n\n\n\n") ;
      printf("  draw_single_searchbin( nb_bin, nj_bin, mht_bin, ht_bin, infile1, infile2 ) \n\n") ;
      printf("\n") ;
      printf("    nb_bin : 0,1,2,3\n") ;
      printf("    nj_bin : 1,2,3,4,5\n") ;
      printf("    mht_bin : 1,2,3,4\n") ;
      printf("    ht_bin : 1,2,3,4\n") ;
      printf("    infile1 : outputfiles/fill-bg-hists3b-t1bbbbH-postdraw.root  ,  for example\n") ;
      printf("    infile2 : outputfiles/fill-bg-dphi-hists2-postdraw.root  ,  for example\n") ;
      printf("\n\n\n\n") ;
   }

  //===============================================================

   TH1F* get_hist( const char* hname ) {

      ////// TH1F* hp = (TH1F*) gDirectory -> FindObject( hname ) ;
      TH1F* hp = (TH1F*) p_hist_file -> Get( hname ) ;
      if ( hp == 0x0 ) { printf("\n\n **** Missing histogram: %s\n\n", hname ) ; gSystem->Exit(-1) ; }

      return hp ;

   } // get_hist

  //===============================================================

   int hist_bin( int mht_bin, int ht_bin ) {

      int rv = 1 + 4*(mht_bin-1) + ht_bin ;

      return rv ;

   } // hist_bin

  //===============================================================

   float get_of_frac( TH1* hp, int xblow, int xbhigh ) {
      float of  = hp -> Integral( xblow, xbhigh ) ;
      float all = hp -> Integral( 1, xbhigh ) ;
      float of_frac = 0. ;
      if ( all > 0 ) of_frac = of / all ;
      return of_frac ;
   } // get_of_frac

  //===============================================================

   void fill_fb_index_arrays( ) {

      if ( sb_mht_ht_ind == 1 ) {

         sprintf( sb_name, "HT1_MHT1_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 2 ) {

         sprintf( sb_name, "HT2_MHT1_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 3 ) {

         sprintf( sb_name, "HT3_MHT1_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 1 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 2 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 4 ) {

         sprintf( sb_name, "HT12_MHT2_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 1 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 5 ) {

         sprintf( sb_name, "HT3_MHT2_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 3 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      } else if ( sb_mht_ht_ind == 6 ) {

         sprintf( sb_name, "HT23_MHT3_Nj%d_Nb%d", sb_nj_ind, fb_nb_ind ) ;
         int i=0;
         if ( sb_nj_ind == 1 ) {
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 1 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 2 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 3 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 3 ; i++ ;
         } else if ( sb_nj_ind == 2 ) {
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 4 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 4 ; i++ ;
         } else if ( sb_nj_ind == 3 ) {
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 2 ; fb_nj_ind[i] = 5 ; i++ ;
            fb_mht_ind[i] = 4 ; fb_ht_ind[i] = 3 ; fb_nj_ind[i] = 5 ; i++ ;
         }
         n_fb_in_sb = i ;

      }

   } // fill_fb_index_arrays

  //===============================================================

   void collect_dphi_hists() {

     char hname[1000] ;
     char newhname[1000] ;

     for ( int fbi=0; fbi<n_fb_in_sb; fbi++ ) {

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_qcd_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_qcd = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_qcd_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_qcd = (TH1F*) h_zl_qcd -> Clone( newhname ) ; } else { h_dphi_zl_qcd -> Add( h_zl_qcd ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_lostlep_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_lostlep = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_lostlep_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_lostlep = (TH1F*) h_zl_lostlep -> Clone( newhname ) ; } else { h_dphi_zl_lostlep -> Add( h_zl_lostlep ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_znunu_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_znunu = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_znunu_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_znunu = (TH1F*) h_zl_znunu -> Clone( newhname ) ; } else { h_dphi_zl_znunu -> Add( h_zl_znunu ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_zl_%s_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], signal_name, fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_zl_sig = get_hist( hname ) ;
        sprintf( newhname, "h_%s_zl_sig_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_zl_sig = (TH1F*) h_zl_sig -> Clone( newhname ) ; } else { h_dphi_zl_sig -> Add( h_zl_sig ) ; }


        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_sl_qcd_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_sl_qcd = get_hist( hname ) ;
        sprintf( newhname, "h_%s_sl_qcd_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_sl_qcd = (TH1F*) h_sl_qcd -> Clone( newhname ) ; } else { h_dphi_sl_qcd -> Add( h_sl_qcd ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_sl_lostlep_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_sl_lostlep = get_hist( hname ) ;
        sprintf( newhname, "h_%s_sl_lostlep_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_sl_lostlep = (TH1F*) h_sl_lostlep -> Clone( newhname ) ; } else { h_dphi_sl_lostlep -> Add( h_sl_lostlep ) ; }

        sprintf( hname, "h_%svsmhtvsht_nb%d_nj%d_sl_znunu_ht%d_mht%d", varname, fb_nb_ind, fb_nj_ind[fbi], fb_ht_ind[fbi], fb_mht_ind[fbi] ) ;
        TH1F* h_sl_znunu = get_hist( hname ) ;
        sprintf( newhname, "h_%s_sl_znunu_%s", varname, sb_name ) ;
        if ( fbi == 0 ) { h_dphi_sl_znunu = (TH1F*) h_sl_znunu -> Clone( newhname ) ; } else { h_dphi_sl_znunu -> Add( h_sl_znunu ) ; }


      //-------

        int nbins = h_zl_qcd -> GetNbinsX() ;
        int cutbin = h_zl_qcd -> FindBin( mdp_cut+1e-5 ) ;

        n_zl_qcd_val[fbi] = h_zl_qcd -> IntegralAndError( cutbin, nbins, n_zl_qcd_err[fbi] ) ;
        n_zl_lostlep_val[fbi] = h_zl_lostlep -> IntegralAndError( cutbin, nbins, n_zl_lostlep_err[fbi] ) ;
        n_zl_znunu_val[fbi] = h_zl_znunu -> IntegralAndError( cutbin, nbins, n_zl_znunu_err[fbi] ) ;

        n_ldp_qcd_val[fbi] = h_zl_qcd -> IntegralAndError( 1, cutbin-1, n_ldp_qcd_err[fbi] ) ;
        n_ldp_lostlep_val[fbi] = h_zl_lostlep -> IntegralAndError( 1, cutbin-1, n_ldp_lostlep_err[fbi] ) ;
        n_ldp_znunu_val[fbi] = h_zl_znunu -> IntegralAndError( 1, cutbin-1, n_ldp_znunu_err[fbi] ) ;

        n_sl_qcd_val[fbi] = h_sl_qcd -> IntegralAndError( cutbin, nbins, n_sl_qcd_err[fbi] ) ;
        n_sl_lostlep_val[fbi] = h_sl_lostlep -> IntegralAndError( cutbin, nbins, n_sl_lostlep_err[fbi] ) ;
        n_sl_znunu_val[fbi] = h_sl_znunu -> IntegralAndError( cutbin, nbins, n_sl_znunu_err[fbi] ) ;

        n_slldp_qcd_val[fbi] = h_sl_qcd -> IntegralAndError( 1, cutbin-1, n_slldp_qcd_err[fbi] ) ;
        n_slldp_lostlep_val[fbi] = h_sl_lostlep -> IntegralAndError( 1, cutbin-1, n_slldp_lostlep_err[fbi] ) ;
        n_slldp_znunu_val[fbi] = h_sl_znunu -> IntegralAndError( 1, cutbin-1, n_slldp_znunu_err[fbi] ) ;

        n_sl_val[fbi]  = n_sl_lostlep_val[fbi] + n_sl_qcd_val[fbi] + n_sl_znunu_val[fbi] ;
        n_ldp_val[fbi] = n_ldp_lostlep_val[fbi] + n_ldp_qcd_val[fbi] + n_ldp_znunu_val[fbi] ;
        n_zl_val[fbi]  = n_zl_lostlep_val[fbi] + n_zl_qcd_val[fbi] + n_zl_znunu_val[fbi] ;

     } // fbi

     char htitle[1000] ;

     sprintf( htitle, "ZL QCD %s", sb_name ) ;
     h_dphi_zl_qcd -> SetTitle( htitle ) ;
     sprintf( htitle, "ZL LostLep %s", sb_name ) ;
     h_dphi_zl_lostlep -> SetTitle( htitle ) ;
     sprintf( htitle, "ZL Znunu %s", sb_name ) ;
     h_dphi_zl_znunu -> SetTitle( htitle ) ;
     sprintf( htitle, "ZL %s %s", signal_name, sb_name ) ;
     h_dphi_zl_sig -> SetTitle( htitle ) ;

     sprintf( htitle, "SL QCD %s", sb_name ) ;
     h_dphi_sl_qcd -> SetTitle( htitle ) ;
     sprintf( htitle, "SL LostLep %s", sb_name ) ;
     h_dphi_sl_lostlep -> SetTitle( htitle ) ;
     sprintf( htitle, "SL Znunu %s", sb_name ) ;
     h_dphi_sl_znunu -> SetTitle( htitle ) ;


   } // collect_dphi_hists

  //===============================================================

   void calc_qcd() {

     float total_pred_qcd_zl_val(0.) ;

     for ( int fbi=0; fbi<n_fb_in_sb; fbi++ ) {

        int mht_bin = fb_mht_ind[fbi] ;
        int ht_bin = fb_ht_ind[fbi] ;
        int nj_bin = fb_nj_ind[fbi] ;

        printf("\n\n") ;
        sprintf( text_line, "====== Fine Bin: HT%d MHT%d Njet%d Nb%d", ht_bin, mht_bin, nj_bin, fb_nb_ind ) ;
        printf( "%s\n", text_line ) ;
        text_y_ndc -= text_dy_ndc ;
        ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

        sprintf( text_line, "   lostlep : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f\n",
              n_sl_lostlep_val[fbi], n_sl_lostlep_err[fbi],
              n_slldp_lostlep_val[fbi], n_slldp_lostlep_err[fbi],
              n_ldp_lostlep_val[fbi], n_ldp_lostlep_err[fbi],
              n_zl_lostlep_val[fbi], n_zl_lostlep_err[fbi] ) ;
        text_y_ndc -= text_dy_ndc ;
        ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

        sprintf( text_line, "   QCD     : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f\n",
              n_sl_qcd_val[fbi], n_sl_qcd_err[fbi],
              n_slldp_qcd_val[fbi], n_slldp_qcd_err[fbi],
              n_ldp_qcd_val[fbi], n_ldp_qcd_err[fbi],
              n_zl_qcd_val[fbi], n_zl_qcd_err[fbi] ) ;
        text_y_ndc -= text_dy_ndc ;
        ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

        sprintf( text_line, "   Znunu   : SL    %7.1f +/- %5.1f   |  SLLDP  %7.1f +/- %5.1f   |  LDP   %7.1f +/- %5.1f   |   ZL   %7.1f +/- %5.1f\n",
              n_sl_znunu_val[fbi], n_sl_znunu_err[fbi],
              n_slldp_znunu_val[fbi], n_slldp_znunu_err[fbi],
              n_ldp_znunu_val[fbi], n_ldp_znunu_err[fbi],
              n_zl_znunu_val[fbi], n_zl_znunu_err[fbi] ) ;
        text_y_ndc -= text_dy_ndc ;
        ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;


        float Rqcd_val = Kqcd_ht_val[ht_bin-1] * Kqcd_mht_val[mht_bin-1] * Kqcd_nj_val[nj_bin-1] ;
        float Rqcd_err = Rqcd_val * sqrt( pow( Kqcd_ht_err[ht_bin-1]/Kqcd_ht_val[ht_bin-1], 2 )
                                        + pow( Kqcd_mht_err[mht_bin-1]/Kqcd_mht_val[mht_bin-1], 2 )
                                        + pow( Kqcd_nj_err[nj_bin-1]/Kqcd_nj_val[nj_bin-1], 2 ) ) ;

        printf("\n") ;
        sprintf( text_line, "    Rqcd = Kqcd_ht * Kqcd_mht * Kqcd_nj = (%.3f +/- %.3f) * (%.3f +/- %.3f) * (%.3f +/- %.3f) = %.3f +/- %.3f",
              Kqcd_ht_val[ht_bin-1], Kqcd_ht_err[ht_bin-1],
              Kqcd_mht_val[mht_bin-1], Kqcd_mht_err[mht_bin-1],
              Kqcd_nj_val[nj_bin-1], Kqcd_nj_err[nj_bin-1],
              Rqcd_val, Rqcd_err ) ;
        text_y_ndc -= text_dy_ndc ;
        ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;


        float pred_qcd_zl_val = n_ldp_qcd_val[fbi] * Rqcd_val ;

        float pred_qcd_zl_err = 0. ;
        if ( n_ldp_qcd_val[fbi] > 0. ) {
           pred_qcd_zl_err = pred_qcd_zl_val * sqrt( pow( Rqcd_err / Rqcd_val, 2 ) + pow( sqrt(n_ldp_val[fbi])/n_ldp_qcd_val[fbi], 2 ) ) ;
        } else {
           pred_qcd_zl_err = pred_qcd_zl_val * sqrt( pow( Rqcd_err / Rqcd_val, 2 ) + n_ldp_val[fbi] ) ;
        }

        sprintf( text_line, "    Predicted ZL QCD : QCD_LDP * Rqcd_val = (%.2f +/- %.2f) * (%.3f +/- %.3f) =  %.2f +/- %.2f\n",
                  n_ldp_qcd_val[fbi], sqrt(n_ldp_val[fbi]),
                  Rqcd_val, Rqcd_err,
                  pred_qcd_zl_val, pred_qcd_zl_err ) ;
        text_y_ndc -= text_dy_ndc ;
        ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

        float error_ratio(1.) ;
        if ( n_zl_val[fbi] > 0 ) error_ratio = pred_qcd_zl_err / sqrt( n_zl_val[fbi] ) ;

        printf("   Error ratio ZL QCD / sqrt(N_ZL):  %.2f / %.2f = %.2f\n", pred_qcd_zl_err, sqrt( n_zl_val[fbi] ), error_ratio ) ;

        total_pred_qcd_zl_val += pred_qcd_zl_val ;

     } // fbi

    //-- Calculate the uncertainty, taking into account correlations between estimates for the fine bins.

     float total_pred_qcd_zl_err = calc_qcd_err() ;

     sprintf(text_line, "  Total predicted ZL QCD : %.2f +/- %.2f", total_pred_qcd_zl_val, total_pred_qcd_zl_err ) ;
     text_y_ndc -= text_dy_ndc ;
     ttext -> DrawTextNDC( text_x_ndc, text_y_ndc, text_line ) ;

   } // calc_qcd


  //===============================================================

   float calc_qcd_err() {

   //++++++++ The order matters for n_ldp array!  Check for consistency with fill_fb_index_arrays.
     float err2(0.) ;

     if ( sb_mht_ht_ind == 1 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[0] * Kqcd_nj_val[2] ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[0] * Kqcd_nj_val[2] ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_nj_val[2] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[0] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[3] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[0] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[4] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[0] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 2 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[1] * Kqcd_nj_val[2] ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[1] * Kqcd_nj_val[2] ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_nj_val[2] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[3] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[4] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[1] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 3 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[2] * Kqcd_nj_val[2] ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[2] * Kqcd_nj_val[2] ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_nj_val[2] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[0]  ,  2 ) ;

           //-- Kqcd_mht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[1]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[0] * Kqcd_ht_val[2]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[1] * Kqcd_ht_val[2] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 4 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[0] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[0] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[1] * Kqcd_nj_val[2] ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[2] * Kqcd_nj_val[2] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[2] * Kqcd_nj_val[2] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[2]   ,  2 ) ;

        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[3] ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[3] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;


           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[3]   ,  2 ) ;

        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[0] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[1] * Kqcd_nj_val[4] ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[4] ) * Kqcd_ht_err[0]   ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;


           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[1] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     } else if ( sb_mht_ht_ind == 5 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[2] * Kqcd_nj_val[2]  ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_nj_val[2] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[2]   ,  2 ) ;


        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[3]   ,  2 ) ;


        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[2]  ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[2] * Kqcd_ht_val[2] ) * Kqcd_nj_err[4]   ,  2 ) ;


        }

     } else if ( sb_mht_ht_ind == 6 ) {

        if ( sb_nj_ind == 1 ) {

           //-- Kqcd_mht4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_ht_val[1] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_ht_val[2] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_ht_val[1] * Kqcd_nj_val[2]
                          + n_ldp_qcd_val[5] * Kqcd_ht_val[2] * Kqcd_nj_val[2]  ) * Kqcd_mht_err[3]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[2] * Kqcd_mht_val[3] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[4] * Kqcd_mht_val[3] * Kqcd_nj_val[2] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_nj_val[0]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[3] * Kqcd_nj_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[3] * Kqcd_nj_val[2] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj1
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[0]   ,  2 ) ;

           //-- Kqcd_nj2
           err2 += pow( (   n_ldp_qcd_val[2] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[3] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[1]   ,  2 ) ;

           //-- Kqcd_nj3
           err2 += pow( (   n_ldp_qcd_val[4] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[5] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[2]   ,  2 ) ;



        } else if ( sb_nj_ind == 2 ) {

           //-- Kqcd_mht4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[3]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[3]  ) * Kqcd_mht_err[3]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_nj_val[3] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_nj_val[3] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[3]   ,  2 ) ;


        } else if ( sb_nj_ind == 3 ) {

           //-- Kqcd_mht4
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_ht_val[1] * Kqcd_nj_val[4]
                          + n_ldp_qcd_val[1] * Kqcd_ht_val[2] * Kqcd_nj_val[4]  ) * Kqcd_mht_err[3]  ,  2 ) ;

           //-- Kqcd_ht2
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_nj_val[4] ) * Kqcd_ht_err[1]   ,  2 ) ;

           //-- Kqcd_ht3
           err2 += pow( (   n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_nj_val[4] ) * Kqcd_ht_err[2]   ,  2 ) ;

           //-- Kqcd_nj5
           err2 += pow( (   n_ldp_qcd_val[0] * Kqcd_mht_val[3] * Kqcd_ht_val[1]
                          + n_ldp_qcd_val[1] * Kqcd_mht_val[3] * Kqcd_ht_val[2] ) * Kqcd_nj_err[4]   ,  2 ) ;

        }

     }

     //-- LDP stats
     for ( int fbi=0; fbi<n_fb_in_sb; fbi++ ) {

      //-- stat error on N_ldp
        if ( n_ldp_val[fbi] > 0 ) {
           err2 += n_ldp_val[fbi] * pow( Kqcd_mht_val[fbi] * Kqcd_ht_val[fbi] * Kqcd_nj_val[fbi] , 2 ) ;
        } else {
           err2 += pow( Kqcd_mht_val[fbi] * Kqcd_ht_val[fbi] * Kqcd_nj_val[fbi] , 2 ) ;
        }

      //-- stat error on lost-lep subtraction.
        err2 += n_ldp_lostlep_val[fbi] * pow( Kqcd_mht_val[fbi] * Kqcd_ht_val[fbi] * Kqcd_nj_val[fbi] , 2 ) ;

     } // fbi

     return sqrt(err2) ;

   } // calc_qcd_err

  //===============================================================

   void fill_kqcd_arrays() {

      if ( strcmp( varname, "mdpn" ) == 0 ) {

      //-- these are for mdpN > 4.0
   ///// Kqcd_ht_val[0] = 0.268 ;  Kqcd_ht_err[0] = 0.015 ;
   ///// Kqcd_ht_val[1] = 0.200 ;  Kqcd_ht_err[1] = 0.009 ;
   ///// Kqcd_ht_val[2] = 0.096 ;  Kqcd_ht_err[2] = 0.007 ;

   ///// Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
   ///// Kqcd_mht_val[1] = 1.34 ; Kqcd_mht_err[1] = 0.1 ;
   ///// Kqcd_mht_val[2] = 2.54 ; Kqcd_mht_err[2] = 0.9 ;
   ///// Kqcd_mht_val[3] = 3.5  ; Kqcd_mht_err[3] = 1.5 ;

   ///// Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
   ///// Kqcd_nj_val[1] = 0.92 ; Kqcd_nj_err[1] = 0.06 ;
   ///// Kqcd_nj_val[2] = 0.62 ; Kqcd_nj_err[2] = 0.07 ;
   ///// Kqcd_nj_val[3] = 0.53 ; Kqcd_nj_err[3] = 0.11 ;
   ///// Kqcd_nj_val[4] = 0.64 ; Kqcd_nj_err[4] = 0.41 ;

      //-- these are for mdpN > 6.0
         Kqcd_ht_val[0] = 0.140 ;  Kqcd_ht_err[0] = 0.012 ;
         Kqcd_ht_val[1] = 0.096 ;  Kqcd_ht_err[1] = 0.007 ;
         Kqcd_ht_val[2] = 0.038 ;  Kqcd_ht_err[2] = 0.005 ;

         Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
         Kqcd_mht_val[1] = 1.29 ; Kqcd_mht_err[1] = 0.2 ;
         Kqcd_mht_val[2] = 2.09 ; Kqcd_mht_err[2] = 1.1 ;
         Kqcd_mht_val[3] = 3.3  ; Kqcd_mht_err[3] = 1.5 ;

         Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
         Kqcd_nj_val[1] = 0.69 ; Kqcd_nj_err[1] = 0.08 ;
         Kqcd_nj_val[2] = 0.69 ; Kqcd_nj_err[2] = 0.13 ;
         Kqcd_nj_val[3] = 0.50 ; Kqcd_nj_err[3] = 0.23 ;
         Kqcd_nj_val[4] = 1.20 ; Kqcd_nj_err[4] = 0.90 ;

      //-- these are for mdpN > 6.0, with aggressive isotrk veto
  ///    Kqcd_ht_val[0] = 0.128 ;  Kqcd_ht_err[0] = 0.012 ;
  ///    Kqcd_ht_val[1] = 0.089 ;  Kqcd_ht_err[1] = 0.007 ;
  ///    Kqcd_ht_val[2] = 0.040 ;  Kqcd_ht_err[2] = 0.005 ;

  ///    Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
  ///    Kqcd_mht_val[1] = 1.00 ; Kqcd_mht_err[1] = 0.3 ;
  ///    Kqcd_mht_val[2] = 2.10 ; Kqcd_mht_err[2] = 2.0 ;
  ///    Kqcd_mht_val[3] = 3.0  ; Kqcd_mht_err[3] = 2.0 ;

  ///    Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
  ///    Kqcd_nj_val[1] = 1.17 ; Kqcd_nj_err[1] = 0.26 ;
  ///    Kqcd_nj_val[2] = 0.30 ; Kqcd_nj_err[2] = 0.38 ;
  ///    Kqcd_nj_val[3] = 0.31 ; Kqcd_nj_err[3] = 0.52 ;
  ///    Kqcd_nj_val[4] = 0.31 ; Kqcd_nj_err[4] = 0.70 ;

      } else {

 ///  //-- these are for plain mdp > 0.3
 ///     Kqcd_ht_val[0] = 0.068 ;  Kqcd_ht_err[0] = 0.010 ;
 ///     Kqcd_ht_val[1] = 0.074 ;  Kqcd_ht_err[1] = 0.010 ;
 ///     Kqcd_ht_val[2] = 0.105 ;  Kqcd_ht_err[2] = 0.020 ;

 ///     Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
 ///     Kqcd_mht_val[1] = 0.28 ; Kqcd_mht_err[1] = 0.1 ;
 ///     Kqcd_mht_val[2] = 0.33 ; Kqcd_mht_err[2] = 0.3 ;
 ///     Kqcd_mht_val[3] = 0.61 ; Kqcd_mht_err[3] = 0.7 ;

 ///     Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
 ///     Kqcd_nj_val[1] = 1.38 ; Kqcd_nj_err[1] = 0.4 ;
 ///     Kqcd_nj_val[2] = 1.66 ; Kqcd_nj_err[2] = 1.0 ;
 ///     Kqcd_nj_val[3] = 1.8  ; Kqcd_nj_err[3] = 1.0 ;
 ///     Kqcd_nj_val[4] = 3.0  ; Kqcd_nj_err[4] = 2.5 ;

      //-- these are for plain mdp > 0.4
         Kqcd_ht_val[0] = 0.046 ;  Kqcd_ht_err[0] = 0.007 ;
         Kqcd_ht_val[1] = 0.041 ;  Kqcd_ht_err[1] = 0.004 ;
         Kqcd_ht_val[2] = 0.045 ;  Kqcd_ht_err[2] = 0.005 ;

         Kqcd_mht_val[0] = 1.   ; Kqcd_mht_err[0] = 0. ;
         Kqcd_mht_val[1] = 0.36 ; Kqcd_mht_err[1] = 0.11 ;
         Kqcd_mht_val[2] = 0.43 ; Kqcd_mht_err[2] = 0.56 ;
         Kqcd_mht_val[3] = 0.86 ; Kqcd_mht_err[3] = 2.0 ;

         Kqcd_nj_val[0] = 1.   ; Kqcd_nj_err[0] = 0. ;
         Kqcd_nj_val[1] = 1.89 ; Kqcd_nj_err[1] = 0.3 ;
         Kqcd_nj_val[2] = 2.02 ; Kqcd_nj_err[2] = 0.4 ;
         Kqcd_nj_val[3] = 2.66 ; Kqcd_nj_err[3] = 0.5 ;
         Kqcd_nj_val[4] = 3.5  ; Kqcd_nj_err[4] = 2.5 ;

      }

   } // fill_kqcd_arrays

  //===============================================================






