

#include "TChain.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TSystem.h"

#include "histio.c"


  //--------------------------------------------------------------------------------------------------------

   int nhtbins = 4 ;   // first bin is not used in analysis.
   double htbins[5] = { 0., 500., 800., 1200., 20000. } ;  // first bin is not used in analysis.

   int nmhtbins = 5 ;  // first bin is not used in analysis.
   double mhtbins[6] = { 0., 200., 300., 500., 750., 20000. } ;  // first bin is not used in analysis.

   int nnbbins = 4 ;
   double nbbins[5] = { -0.5, 0.5, 1.5, 2.5, 20. } ;

   int nnjetbins = 6 ; // first bin is not used in analysis.
   double njetbins[7] = { 0., 3.5, 4.5, 5.5, 6.5, 8.5, 20. } ;  // first bin is not used in analysis.

  //--------------------------------------------------------------------------------------------------------

   int n_merged_bins ;
   int merged_bin_index[8][8][8][8] ; // nbi, nji, mbi, hbi
   char merged_bin_label[300][20] ;

  //--------------------------------------------------------------------------------------------------------

   int n_samples ;
   char sname[100][100] ;
   int scolor[100] ;
   bool isll[100] ;
   bool issignal[100] ;
   char signal_name[100] ;

   double xsec_over_ngen ;

   void add_overflow_to_last_bin( TH1* hp ) ;
   void draw_ht_vs_mht_canvas( const char* hname, const char* plotname, const char* sample_name, const char* config_string, int fill_color, float xmin=0., float xmax=-1., bool is_signal=false ) ;
   void draw_ht_vs_mht_stack_canvas( const char* hname, const char* plotname, const char* config_string, float xmin=0., float xmax=-1. ) ;
   void draw_ht_vs_mht_canvas_searchbins( const char* hname, const char* plotname, const char* sample_name, const char* config_string, float xmin=0., float xmax=-1. ) ;
   void draw_ht_vs_mht_stack_canvas_searchbins( const char* hname, const char* plotname, const char* config_string, float xmin=0., float xmax=-1. ) ;
   void draw_stack( const char* hname_base, const char* hname_end, float xmin=0., float xmax=-1. ) ;
   void add_lostlep_samples( const char* hnamebase ) ;
   void add_nb_hists( const char* hname ) ;
   void add_njet_hists( const char* hname ) ;
   float get_of_frac( TH1* hp, int xblow, int xbhigh ) ;
   TH3* get_hist_3d( const char* hname ) ;
   TH1* get_hist_1d( const char* hname ) ;


   TCanvas* can ;

  //--------------------------------------------------------------------------------------------------------

   void draw_bgsig_dphi_hists2( const char* infile_bg = "outputfiles/fill-bg-dphi-hists2.root",
                                const char* infile_sig = "outputfiles/fill-sig-dphi-hists2.root",
                                const char* signame = "t1bbbbH" ) {

      if ( strcmp( signame, "t1bbbbH" )==0 ) {
         xsec_over_ngen = 0.01419 / 105149 ;
      } else if ( strcmp( signame, "t1bbbbC" )==0 ) {
         xsec_over_ngen = 0.3254  /  97134 ;
      } else if ( strcmp( signame, "t1ttttH" )==0 ) {
         xsec_over_ngen = 0.01419 / 105679 ;
      } else if ( strcmp( signame, "t1ttttC" )==0 ) {
         xsec_over_ngen = 0.08564 / 100322 ;
      } else if ( strcmp( signame, "t1qqqqH" )==0 ) {
         xsec_over_ngen = 0.02530 / 102891 ;
      } else if ( strcmp( signame, "t1qqqqC" )==0 ) {
         xsec_over_ngen = 0.3254  /  96681 ;
      } else {
         printf("\n\n *** Unknown signal name : %s\n\n", signame ) ; return ;
      }


      sprintf( signal_name, "%s", signame ) ;

      gStyle -> SetPadBottomMargin( 0.20 ) ;
      gStyle -> SetOptStat(0) ;

      //setup_merged_bins() ;

      gDirectory -> Delete( "h*" ) ;
      printf("\n\n Reading in BG 3D mdp vs MHT vs HT histograms from %s\n\n", infile_bg ) ;
      loadHist( infile_bg ) ;
      printf("\n\n Reading in Signal 3D mdp vs MHT vs HT histograms from %s\n\n", infile_sig ) ;
      loadHist( infile_sig ) ;

    //--------
      int n_selections = 2 ;
      char selname[2][100] ;
      sprintf( selname[0], "zl" ) ;
      sprintf( selname[1], "sl" ) ;
    //--------

      int si(0) ;
      sprintf( sname[si], "qcd" )     ; isll[si] = false ; scolor[si] = kRed + 3 ; issignal[si] = false ; si++ ;
      sprintf( sname[si], "ttbar" )   ; isll[si] = true  ; scolor[si] = kBlue - 4 ; issignal[si] = false ; si++ ;
      sprintf( sname[si], "wjets" )   ; isll[si] = true  ; scolor[si] = kGreen + 2 ; issignal[si] = false ; si++ ;
      sprintf( sname[si], "sngltop" ) ; isll[si] = true  ; scolor[si] = kYellow - 7 ; issignal[si] = false ; si++ ;
      sprintf( sname[si], "znunu" )   ; isll[si] = false ; scolor[si] = kOrange - 3 ; issignal[si] = false ; si++ ;
      sprintf( sname[si], "%s", signame )   ; isll[si] = false ; scolor[si] = kMagenta - 3 ; issignal[si] = true ; si++ ;
      n_samples = si ;

      can = (TCanvas*) gDirectory -> FindObject( "can_draw_bgsig_dphi_hists" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_draw_bgsig_dphi_hists", "BG hists", 900, 1100 ) ;
      can -> Clear() ;

   //-------------------------------

      for ( int si=0; si<n_samples; si++ ) {

         for ( int selind=0; selind<n_selections; selind++ ) {

            for ( int nbi=0; nbi<nnbbins; nbi++ ) {

               for ( int nji=1; nji<nnjetbins; nji++ ) {

/////////////     char configstr[1000] ;
/////////////     sprintf( configstr, "nb%d_nj%d_%s", nbi, nji, selname[selind] ) ;

/////////////     char hname[1000] ;

/////////////     sprintf( hname, "h_mdpvsmhtvsht_%s_%s", configstr, sname[si] ) ;
/////////////     draw_ht_vs_mht_canvas( hname, "dphi", sname[si], configstr, scolor[si] ) ;

/////////////     sprintf( hname, "h_mdpnvsmhtvsht_%s_%s", configstr, sname[si] ) ;
/////////////     draw_ht_vs_mht_canvas( hname, "dphin", sname[si], configstr, scolor[si] ) ;


               } // nji

            } // nbi.

         } // selind.

      } // si


   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //   Split by Njet and Nb


      {
         for ( int selind=0; selind<n_selections; selind++ ) {
         //for ( int selind=0; selind<1; selind++ ) {

            for ( int nbi=0; nbi<nnbbins; nbi++ ) {
            //for ( int nbi=0; nbi<1; nbi++ ) {

               for ( int nji=1; nji<nnjetbins; nji++ ) {
               //for ( int nji=1; nji<2; nji++ ) {

                  char configstr[1000] ;
                  sprintf( configstr, "nb%d_nj%d_%s", nbi, nji, selname[selind] ) ;

                  char hnamebase[1000] ;
                  char hname[1000] ;

               //------

                  sprintf( hnamebase, "h_mdpvsmhtvsht_%s", configstr ) ;

                  add_lostlep_samples( hnamebase ) ;
                  sprintf( hname, "%s_lostlep", hnamebase ) ;
                  draw_ht_vs_mht_canvas( hname, "dphi", "lostlep", configstr, scolor[1]-2 ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "lostlep", configstr, scolor[1]-2 ) ;

                  sprintf( hname, "%s_qcd", hnamebase ) ;
                  draw_ht_vs_mht_canvas( hname, "dphi", "qcd", configstr, scolor[0]-2 ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "qcd", configstr, scolor[0]-2 ) ;

                  sprintf( hname, "%s_znunu", hnamebase ) ;
                  draw_ht_vs_mht_canvas( hname, "dphi", "znunu", configstr, scolor[4]-2 ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "znunu", configstr, scolor[4]-2 ) ;

                  sprintf( hname, "%s_%s", hnamebase, signame ) ;
                  draw_ht_vs_mht_canvas( hname, "dphi", signame, configstr, scolor[5]-3, 0., -1, true ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", signame, configstr, scolor[5]-3 ) ;

                  draw_ht_vs_mht_stack_canvas( hnamebase, "dphi", configstr ) ;
                  draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphi", configstr, 0., 1.0 ) ;
                  draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphi-searchbins", configstr ) ;
                  draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphi-searchbins", configstr, 0., 1.0 ) ;


               //------

                  sprintf( hnamebase, "h_mdpnvsmhtvsht_%s", configstr ) ;

                  add_lostlep_samples( hnamebase ) ;
                  sprintf( hname, "%s_lostlep", hnamebase ) ;
                  draw_ht_vs_mht_canvas( hname, "dphin", "lostlep", configstr, scolor[1]-2 ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "lostlep", configstr, scolor[1]-2 ) ;

                  sprintf( hname, "%s_qcd", hnamebase ) ;
                  draw_ht_vs_mht_canvas( hname, "dphin", "qcd", configstr, scolor[0]-2 ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "qcd", configstr, scolor[0]-2 ) ;

                  sprintf( hname, "%s_znunu", hnamebase ) ;
                  draw_ht_vs_mht_canvas( hname, "dphin", "znunu", configstr, scolor[4]-2 ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "znunu", configstr, scolor[4]-2 ) ;

                  sprintf( hname, "%s_%s", hnamebase, signame ) ;
                  draw_ht_vs_mht_canvas( hname, "dphin", signame, configstr, scolor[5]-3, 0., -1, true ) ;
                  draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", signame, configstr, scolor[5]-3 ) ;

                  draw_ht_vs_mht_stack_canvas( hnamebase, "dphin", configstr ) ;
                  draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphin", configstr, 0., 20 ) ;
                  draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphin-searchbins", configstr ) ;
                  draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphin-searchbins", configstr, 0., 20 ) ;


               } // nji

            } // nbi.

         } // selind.
      }



   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //   Sum over Nb


//    {
//       for ( int selind=0; selind<n_selections; selind++ ) {

//          for ( int nji=1; nji<nnjetbins; nji++ ) {
//          //for ( int nji=1; nji<2; nji++ ) {

//             char configstr[1000] ;
//             char sumconfigstr[1000] ;
//             sprintf( configstr, "nb0_nj%d_%s", nji, selname[selind] ) ;

//             char hnamebase[1000] ;
//             char hname[1000] ;


//          //------

//             sprintf( hnamebase, "h_mdpvsmhtvsht_%s", configstr ) ;

//             sprintf( hname, "%s_lostlep", hnamebase ) ;
//             add_nb_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nbsum_nj%d_%s_lostlep", nji, selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_nj%d_%s", nji, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "lostlep", sumconfigstr, scolor[1]-2 ) ;

//             sprintf( hname, "%s_qcd", hnamebase ) ;
//             add_nb_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nbsum_nj%d_%s_qcd", nji, selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_nj%d_%s", nji, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "qcd", sumconfigstr, scolor[0]-2 ) ;

//             sprintf( hname, "%s_znunu", hnamebase ) ;
//             add_nb_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nbsum_nj%d_%s_znunu", nji, selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_nj%d_%s", nji, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "znunu", sumconfigstr, scolor[4]-2 ) ;

//             draw_ht_vs_mht_stack_canvas( hnamebase, "dphi", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphi", sumconfigstr, 0., 1.0 ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphi-searchbins", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphi-searchbins", sumconfigstr, 0., 1.0 ) ;



//          //------

//             sprintf( hnamebase, "h_mdpnvsmhtvsht_%s", configstr ) ;

//             sprintf( hname, "%s_lostlep", hnamebase ) ;
//             add_nb_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nbsum_nj%d_%s_lostlep", nji, selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_nj%d_%s", nji, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "lostlep", sumconfigstr, scolor[1]-2 ) ;

//             sprintf( hname, "%s_qcd", hnamebase ) ;
//             add_nb_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nbsum_nj%d_%s_qcd", nji, selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_nj%d_%s", nji, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "qcd", sumconfigstr, scolor[0]-2 ) ;

//             sprintf( hname, "%s_znunu", hnamebase ) ;
//             add_nb_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nbsum_nj%d_%s_znunu", nji, selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_nj%d_%s", nji, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "znunu", sumconfigstr, scolor[4]-2 ) ;

//             draw_ht_vs_mht_stack_canvas( hnamebase, "dphin", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphin", sumconfigstr, 0., 20 ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphin-searchbins", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphin-searchbins", sumconfigstr, 0., 20 ) ;


//          } // nji

//       } // selind.
//    }


// //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// //   Sum over Njet

//    {
//       for ( int selind=0; selind<n_selections; selind++ ) {

//          for ( int nbi=0; nbi<nnbbins; nbi++ ) {

//             char configstr[1000] ;
//             char sumconfigstr[1000] ;
//             sprintf( configstr, "nb%d_nj1_%s", nbi, selname[selind] ) ;

//             char hnamebase[1000] ;
//             char hname[1000] ;

//          //------

//             sprintf( hnamebase, "h_mdpvsmhtvsht_%s", configstr ) ;

//             sprintf( hname, "%s_lostlep", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nb%d_njsum_%s_lostlep", nbi, selname[selind] ) ;
//             sprintf( sumconfigstr, "nb%d_njsum_%s", nbi, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "lostlep", sumconfigstr, scolor[1]-2 ) ;

//             sprintf( hname, "%s_qcd", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nb%d_njsum_%s_qcd", nbi, selname[selind] ) ;
//             sprintf( sumconfigstr, "nb%d_njsum_%s", nbi, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "qcd", sumconfigstr, scolor[0]-2 ) ;

//             sprintf( hname, "%s_znunu", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nb%d_njsum_%s_znunu", nbi, selname[selind] ) ;
//             sprintf( sumconfigstr, "nb%d_njsum_%s", nbi, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "znunu", sumconfigstr, scolor[4]-2 ) ;

//             draw_ht_vs_mht_stack_canvas( hnamebase, "dphi", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphi", sumconfigstr, 0., 1.0 ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphi-searchbins", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphi-searchbins", sumconfigstr, 0., 1.0 ) ;


//          //------

//             sprintf( hnamebase, "h_mdpnvsmhtvsht_%s", configstr ) ;

//             sprintf( hname, "%s_lostlep", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nb%d_njsum_%s_lostlep", nbi, selname[selind] ) ;
//             sprintf( sumconfigstr, "nb%d_njsum_%s", nbi, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "lostlep", sumconfigstr, scolor[1]-2 ) ;

//             sprintf( hname, "%s_qcd", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nb%d_njsum_%s_qcd", nbi, selname[selind] ) ;
//             sprintf( sumconfigstr, "nb%d_njsum_%s", nbi, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "qcd", sumconfigstr, scolor[0]-2 ) ;

//             sprintf( hname, "%s_znunu", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nb%d_njsum_%s_znunu", nbi, selname[selind] ) ;
//             sprintf( sumconfigstr, "nb%d_njsum_%s", nbi, selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "znunu", sumconfigstr, scolor[4]-2 ) ;

//             draw_ht_vs_mht_stack_canvas( hnamebase, "dphin", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphin", sumconfigstr, 0., 20 ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphin-searchbins", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphin-searchbins", sumconfigstr, 0., 20 ) ;

//          //------


//          } // nji

//       } // selind.
//    }


// //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// //   Sum over Nb and Njet

//    {
//       for ( int selind=0; selind<n_selections; selind++ ) {

//             char configstr[1000] ;
//             char sumconfigstr[1000] ;
//             sprintf( configstr, "nbsum_nj1_%s", selname[selind] ) ;

//             char hnamebase[1000] ;
//             char hname[1000] ;

//          //------

//             sprintf( hnamebase, "h_mdpvsmhtvsht_%s", configstr ) ;

//             sprintf( hname, "%s_lostlep", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nbsum_njsum_%s_lostlep", selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_njsum_%s", selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas( hname, "zoom-dphi", "lostlep", sumconfigstr, scolor[1]-2, 0., 1.0 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "zoom-dphi-searchbins", "lostlep", sumconfigstr, 0., 1.0 ) ;

//             sprintf( hname, "%s_qcd", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nbsum_njsum_%s_qcd", selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_njsum_%s", selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas( hname, "zoom-dphi", "qcd", sumconfigstr, scolor[0]-2, 0., 1.0 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "zoom-dphi-searchbins", "qcd", sumconfigstr, 0., 1.0 ) ;

//             sprintf( hname, "%s_znunu", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpvsmhtvsht_nbsum_njsum_%s_znunu", selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_njsum_%s", selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphi", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphi-searchbins", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas( hname, "zoom-dphi", "znunu", sumconfigstr, scolor[4]-2, 0., 1.0 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "zoom-dphi-searchbins", "znunu", sumconfigstr, 0., 1.0 ) ;

//             draw_ht_vs_mht_stack_canvas( hnamebase, "dphi", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphi", sumconfigstr, 0., 1.0 ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphi-searchbins", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphi-searchbins", sumconfigstr, 0., 1.0 ) ;

//          //------

//             sprintf( hnamebase, "h_mdpnvsmhtvsht_%s", configstr ) ;

//             sprintf( hname, "%s_lostlep", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nbsum_njsum_%s_lostlep", selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_njsum_%s", selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "lostlep", sumconfigstr, scolor[1]-2 ) ;
//             draw_ht_vs_mht_canvas( hname, "zoom-dphin", "lostlep", sumconfigstr, scolor[1]-2, 0., 20. ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "zoom-dphin-searchbins", "lostlep", sumconfigstr, 0., 20. ) ;

//             sprintf( hname, "%s_qcd", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nbsum_njsum_%s_qcd", selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_njsum_%s", selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "qcd", sumconfigstr, scolor[0]-2 ) ;
//             draw_ht_vs_mht_canvas( hname, "zoom-dphin", "qcd", sumconfigstr, scolor[0]-2, 0., 20. ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "zoom-dphin-searchbins", "qcd", sumconfigstr, 0., 20. ) ;

//             sprintf( hname, "%s_znunu", hnamebase ) ;
//             add_njet_hists( hname ) ;
//             sprintf( hname, "h_mdpnvsmhtvsht_nbsum_njsum_%s_znunu", selname[selind] ) ;
//             sprintf( sumconfigstr, "nbsum_njsum_%s", selname[selind] ) ;
//             draw_ht_vs_mht_canvas( hname, "dphin", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "dphin-searchbins", "znunu", sumconfigstr, scolor[4]-2 ) ;
//             draw_ht_vs_mht_canvas( hname, "zoom-dphin", "znunu", sumconfigstr, scolor[4]-2, 0., 20. ) ;
//             draw_ht_vs_mht_canvas_searchbins( hname, "zoom-dphin-searchbins", "znunu", sumconfigstr, 0., 20. ) ;

//             draw_ht_vs_mht_stack_canvas( hnamebase, "dphin", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas( hnamebase, "zoom-dphin", sumconfigstr, 0., 20 ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "dphin-searchbins", sumconfigstr ) ;
//             draw_ht_vs_mht_stack_canvas_searchbins( hnamebase, "zoom-dphin-searchbins", sumconfigstr, 0., 20 ) ;

//          //------

//       } // selind.
//    }


      TString savefilets( infile_bg ) ;
      char newend[100] ;
      sprintf( newend, "-%s-postdraw.root", signame ) ;
      savefilets.ReplaceAll( ".root", newend ) ;
      printf("\n\n Saving all histograms as %s\n", savefilets.Data() ) ;
      saveHist( savefilets.Data(), "h*" ) ;


   } // draw_bgsig_dphi_hists2

  //====================================================================================================================

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

  //====================================================================================================================

   void draw_ht_vs_mht_canvas( const char* hname, const char* plotname, const char* sample_name, const char* config_string, int fill_color, float xmin, float xmax, bool is_signal ) {

      TH3* h3p ;
      char savename[1000] ;
      char label[100] ;
      int ci(1) ;

      TText* label_tt = new TText() ;
      label_tt -> SetTextAlign( 33 ) ;

      printf("\n draw_ht_vs_mht_canvas: Looking for %s\n", hname ) ;
      h3p = get_hist_3d( hname ) ;

      can -> Clear() ;
      can -> Divide( 3, 4 ) ;
      ci = 1 ;
      for ( int mhti=5; mhti>=2; mhti-- ) {
         for ( int hti=2; hti<=4; hti++ ) {
            char hpname[1000] ;
            sprintf( hpname, "%s_ht%d_mht%d", hname, hti-1, mhti-1 ) ;
            TH1D* hp = h3p->ProjectionZ( hpname, hti, hti,  mhti, mhti ) ;
            if ( is_signal ) {
               printf("\n scaling %s by %g\n", hname, xsec_over_ngen ) ;
               hp -> Scale( xsec_over_ngen ) ;
            }
            hp -> SetFillColor( fill_color ) ;
            add_overflow_to_last_bin( hp ) ;
            can -> cd(ci) ;
            if ( xmax > xmin ) hp -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
            hp -> Draw("hist") ;
            hp -> Draw("same") ;
            sprintf( label, "Fine bins: HT%d, MHT%d", hti-1, mhti-1 ) ;
            label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;
            ci++ ;
         } // hti
      } // mhti

      sprintf( savename, "outputfiles/%s-%s_%s.pdf", plotname, sample_name, config_string ) ;
      can -> SaveAs( savename ) ;

   } // draw_ht_vs_mht_canvas



  //====================================================================================================================

   void draw_ht_vs_mht_stack_canvas( const char* hname, const char* plotname, const char* config_string, float xmin, float xmax ) {

      TH3* h3p ;
      char savename[1000] ;
      char label[100] ;
      int ci(1) ;

      TText* label_tt = new TText() ;
      label_tt -> SetTextAlign( 33 ) ;


      printf("\n draw_ht_vs_mht_stack_canvas: Looking for %s\n", hname ) ;

      can -> Clear() ;
      can -> Divide( 3, 4 ) ;
      ci = 1 ;
      for ( int mhti=5; mhti>=2; mhti-- ) {
         for ( int hti=2; hti<=4; hti++ ) {

            can -> cd(ci) ;
            char hname_end[1000] ;
            sprintf( hname_end, "ht%d_mht%d", hti-1, mhti-1 ) ;
            draw_stack( hname, hname_end, xmin, xmax ) ;

            sprintf( label, "Fine bins: HT%d, MHT%d", hti-1, mhti-1 ) ;
            label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;
            ci++ ;

         } // hti
      } // mhti

      sprintf( savename, "outputfiles/%s-stack_%s.pdf", plotname, config_string ) ;
      can -> SaveAs( savename ) ;

   } // draw_ht_vs_mht_canvas



  //====================================================================================================================

  //--- to be called after draw_ht_vs_mht_canvas
   void draw_ht_vs_mht_canvas_searchbins( const char* hname, const char* plotname, const char* sample_name, const char* config_string, float xmin, float xmax ) {

      char savename[1000] ;
      char label[100] ;
      char hpname[1000] ;
      char hnewname[1000] ;
      TH1* hp1 ;
      TH1* hp2 ;
      TH1* hpc ;
      int ci(1) ;

      TText* label_tt = new TText() ;
      label_tt -> SetTextAlign( 33 ) ;

      can -> Clear() ;
      can -> Divide( 3, 4 ) ;

    //----------
      ci = 10 ;
      sprintf( hpname, "%s_ht1_mht1", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht1_mht1a", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT1, MHT1a" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci = 11 ;
      sprintf( hpname, "%s_ht2_mht1", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht2_mht1a", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT2, MHT1a" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci = 12 ;
      sprintf( hpname, "%s_ht3_mht1", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht3_mht1a", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT3, MHT1a" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;





    //----------
      ci =  7 ;
      sprintf( hpname, "%s_ht1_mht2", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht1_mht1b", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT1, MHT1b" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci =  8 ;
      sprintf( hpname, "%s_ht2_mht2", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht2_mht1b", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT2, MHT1b" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci =  9 ;
      sprintf( hpname, "%s_ht3_mht2", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht3_mht1b", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT3, MHT1b" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;





    //----------
      ci =  5 ;
      sprintf( hpname, "%s_ht1_mht3", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hpname, "%s_ht2_mht3", hname ) ;
      hp2 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht12_mht2", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      hpc -> Add( hp2 ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT12, MHT2" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci =  6 ;
      sprintf( hpname, "%s_ht3_mht3", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht3_mht2", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT3, MHT2" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;





    //----------
      ci =  3 ;
      sprintf( hpname, "%s_ht2_mht4", hname ) ;
      hp1 = get_hist_1d( hpname ) ;
      sprintf( hpname, "%s_ht3_mht4", hname ) ;
      hp2 = get_hist_1d( hpname ) ;
      sprintf( hnewname, "%s_searchbin_ht23_mht3", hname ) ;
      hpc = (TH1*) hp1 -> Clone( hnewname ) ;
      hpc -> Add( hp2 ) ;
      can -> cd(ci) ;
      if ( xmax > xmin ) hpc -> GetXaxis() -> SetRangeUser( xmin, xmax ) ;
      hpc -> Draw("hist") ;
      hpc -> Draw("same") ;
      sprintf( label, "Search bins: HT23, MHT3" ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;


      sprintf( savename, "outputfiles/%s-%s_%s.pdf", plotname, sample_name, config_string ) ;
      can -> SaveAs( savename ) ;

   } // draw_ht_vs_mht_canvas_searchbins



  //====================================================================================================================

  //--- to be called after draw_ht_vs_mht_canvas
   void draw_ht_vs_mht_stack_canvas_searchbins( const char* hname, const char* plotname, const char* config_string, float xmin, float xmax ) {

      char savename[1000] ;
      char label[100] ;
      char hpname[1000] ;
      char hsname[1000] ;
      char sbin_name[1000] ;
      TH1* hp_znunu ;
      TH1* hp_lostlep ;
      TH1* hp_qcd ;
      THStack* hs ;
      int ci(1) ;

      TText* label_tt = new TText() ;
      label_tt -> SetTextAlign( 33 ) ;

      can -> Clear() ;
      can -> Divide( 3, 4 ) ;

    //----------
      ci = 10 ;
      sprintf( sbin_name, "searchbin_ht1_mht1a" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;


    //----------
      ci = 11 ;
      sprintf( sbin_name, "searchbin_ht2_mht1a" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci = 12 ;
      sprintf( sbin_name, "searchbin_ht3_mht1a" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;





    //----------
      ci =  7 ;
      sprintf( sbin_name, "searchbin_ht1_mht1b" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci =  8 ;
      sprintf( sbin_name, "searchbin_ht2_mht1b" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci =  9 ;
      sprintf( sbin_name, "searchbin_ht3_mht1b" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;




    //----------
      ci =  5 ;
      sprintf( sbin_name, "searchbin_ht12_mht2" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;

    //----------
      ci =  6 ;
      sprintf( sbin_name, "searchbin_ht3_mht2" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;




    //----------
      ci =  3 ;
      sprintf( sbin_name, "searchbin_ht23_mht3" ) ;
      can -> cd(ci) ;
      draw_stack( hname, sbin_name, xmin, xmax ) ;
      sprintf( label, "Search bins: %s", sbin_name ) ;
      label_tt -> DrawTextNDC( 0.80, 0.85, label ) ;



      sprintf( savename, "outputfiles/%s-stack_%s.pdf", plotname, config_string ) ;
      can -> SaveAs( savename ) ;

   } // draw_ht_vs_mht_canvas_searchbins



  //====================================================================================================================


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
       //////////h_ll -> Sumw2() ;
       TString newtitle = h_ll -> GetTitle() ;
       newtitle.ReplaceAll( "ttbar", "lostlep" ) ;
       h_ll -> SetTitle( newtitle ) ;
       h_ll -> Add( h_wjets ) ;
       h_ll -> Add( h_sngltop ) ;

   } // add_lostlep_samples


  //====================================================================================================================

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

  //====================================================================================================================

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

  //====================================================================================================================

   TH3* get_hist_3d( const char* hname ) {

      TH3* rp = (TH3*) gDirectory -> FindObject( hname ) ;
      if ( rp == 0x0 ) {
         printf("\n\n *** Missing histogram: %s\n\n", hname ) ;
         gSystem -> Exit(-1) ;
      }
      return rp ;

   } // get_hist_3d


  //====================================================================================================================

   TH1* get_hist_1d( const char* hname ) {

      TH1* rp = (TH1*) gDirectory -> FindObject( hname ) ;
      if ( rp == 0x0 ) {
         printf("\n\n *** Missing histogram: %s\n\n", hname ) ;
         gSystem -> Exit(-1) ;
      }
      return rp ;

   } // get_hist_1d


  //====================================================================================================================

   void draw_stack( const char* hname_base, const char* hname_end, float xmin, float xmax ) {

       char hpname[1000] ;

       sprintf( hpname, "%s_znunu_%s", hname_base, hname_end ) ;
       TH1* hp_znunu = get_hist_1d( hpname ) ;

       sprintf( hpname, "%s_lostlep_%s", hname_base, hname_end ) ;
       TH1* hp_lostlep = get_hist_1d( hpname ) ;

       sprintf( hpname, "%s_qcd_%s", hname_base, hname_end ) ;
       TH1* hp_qcd = get_hist_1d( hpname ) ;

       sprintf( hpname, "%s_%s_%s", hname_base, signal_name, hname_end ) ;
       TH1* hp_sig = get_hist_1d( hpname ) ;

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
          float sig_of_frac = get_of_frac( hp_sig, iblow, ibhigh )  ;
          printf(" %s : overflow fracs:  QCD = %.3f ,   lostlep = %.3f ,   Znunu = %.3f ,   %s = %.3f\n",
             hpname, qcd_of_frac, lostlep_of_frac, znunu_of_frac, signal_name, sig_of_frac ) ;
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
          sprintf( label, "%s %.2f", signal_name, sig_of_frac ) ;
          tt -> DrawTextNDC( tx, 0.55, label ) ;
       }
       hs -> Draw("same") ;

       hp_sig -> SetLineColor( kMagenta -3 ) ;
       hp_sig -> SetLineWidth( 2 ) ;
       hp_sig -> SetFillColor( 0 ) ;
       hp_sig -> Draw( "hist same" ) ;


   } // draw_stack

  //====================================================================================================================


   float get_of_frac( TH1* hp, int xblow, int xbhigh ) {
      float of  = hp -> Integral( xblow, xbhigh ) ;
      float all = hp -> Integral( 1, xbhigh ) ;
      float of_frac = 0. ;
      if ( all > 0 ) of_frac = of / all ;
      return of_frac ;
   } // get_of_frac

  //====================================================================================================================


