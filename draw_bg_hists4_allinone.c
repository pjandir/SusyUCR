

#include "TChain.h"
#include "TH1F.h"
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

#include <fstream>

  //--------------------------------------------------------------------------------------------------------

   int nhtbins = 4 ;   // first bin is not used in analysis.
   double htbins[5] = { 0., 500., 800., 1200., 20000. } ;  // first bin is not used in analysis.

   int nmhtbins ;
   double mhtbins[10] ;

   int nnbbins = 4 ;
   double nbbins[5] = { -0.5, 0.5, 1.5, 2.5, 20. } ;

   int nnjetbins = 6 ; // first bin is not used in analysis.
   double njetbins[7] = { 0., 3.5, 4.5, 5.5, 6.5, 8.5, 20. } ;  // first bin is not used in analysis.

  //--------------------------------------------------------------------------------------------------------

   int n_samples ;
   char sname[100][100] ;
   int scolor[100] ;
   bool isll[100] ;

   const int N_FB_MAX(1000) ;
   int n_fine_bins(0) ;
   char fb_name[N_FB_MAX][100] ;
   int fb_sbi[N_FB_MAX] ;
   int fb_nji[N_FB_MAX] ;
   int fb_nbi[N_FB_MAX] ;
   int fb_mhti[N_FB_MAX] ;
   int fb_hti[N_FB_MAX] ;

   int n_search_bins ;
   const int MAX_SB(1000) ;
   char sb_name[MAX_SB][100] ;



   TCanvas* can ;

   void draw_stack( const char* selname ) ;
   void make_legend(  ) ;
   void setup_binmap( const char* lhbuilder_file ) ;
   bool  find_line( ifstream& ifs, const char* key ) ;
   float find_line_val( ifstream& ifs, const char* key ) ;

  //--------------------------------------------------------------------------------------------------------

   void draw_bg_hists4_allinone( const char* infile = "outputfiles/fill-bg-hists-tm4.root",
                                 const char* lhbuilder_file = "fb-sb-bin-map.txt"
                                  ) {

      setup_binmap( lhbuilder_file ) ;

      nmhtbins = 5 ;
      mhtbins[0] = 0. ;
      mhtbins[1] = 200. ;
      mhtbins[2] = 300. ;
      mhtbins[3] = 500. ;
      mhtbins[4] = 750. ;
      mhtbins[5] = 20000. ;

      gStyle -> SetPadBottomMargin( 0.40 ) ;
      gStyle -> SetOptStat(0) ;


      gDirectory -> Delete( "h*" ) ;
      printf("\n\n Reading in 3D histograms from %s\n\n", infile ) ;
      loadHist( infile ) ;

      int n_selections = 4 ;
      char selname[4][100] ;
      sprintf( selname[0], "zl" ) ;
      sprintf( selname[1], "sl" ) ;
      sprintf( selname[2], "ldp" ) ;
      sprintf( selname[3], "slldp" ) ;

      int si(0) ;
      sprintf( sname[si], "qcd" )     ; isll[si] = false ; scolor[si] = kRed + 3 ; si++ ;
      sprintf( sname[si], "ttbar" )   ; isll[si] = true  ; scolor[si] = kBlue - 4 ; si++ ;
      sprintf( sname[si], "wjets" )   ; isll[si] = true  ; scolor[si] = kGreen + 2 ; si++ ;
      sprintf( sname[si], "sngltop" ) ; isll[si] = true  ; scolor[si] = kYellow - 7 ; si++ ;
      sprintf( sname[si], "znunu" )   ; isll[si] = false ; scolor[si] = kOrange - 3 ; si++ ;
      n_samples = si ;

      can = (TCanvas*) gDirectory -> FindObject( "can_draw_bg_hists4_allinone" ) ;
      if ( can == 0x0 ) can = new TCanvas( "can_draw_bg_hists4_allinone", "BG hists", 1500, 750 ) ;
      can -> Clear() ;

      for ( int selind=0; selind<n_selections; selind++ ) {

         draw_stack( selname[selind] ) ;

      } // selind.

      make_legend( ) ;

      TString savefilets( infile ) ;
      char newend[100] ;
      sprintf( newend, "-allinone-postdraw.root" ) ;
      savefilets.ReplaceAll( ".root", newend ) ;
      printf("\n\n Saving all histograms as %s\n", savefilets.Data() ) ;
      saveHist( savefilets.Data(), "h*" ) ;


   } // draw_bg_hists4_allinone

  //====================================================================================================================

   void draw_stack( const char* selname ) {

      printf("\n\n ============== Drawing stack for %s\n\n", selname ) ;

      char hname[1000] ;
      char htitle[1000] ;

      sprintf( hname, "h_stack_%s", selname ) ;
      THStack* h_stack = new THStack( hname, selname ) ;

      TH1F* h_bgsum(0x0) ;
      TH1F* h_llbgsum(0x0) ;

      TH1F* h1d_samples[100] ;

      for ( int si=0; si<n_samples; si++ ) {

         printf("          Sample %s\n", sname[si] ) ;

         char hname1d[1000] ;
         char htitle1d[1000] ;
         sprintf( hname1d, "h_1d_allinone_%s_%s", selname, sname[si] ) ;
         sprintf( htitle1d, "%s, %s", selname, sname[si] ) ;
         TH1F* h1p = new TH1F( hname1d, htitle1d, n_fine_bins-20, 0.5, n_fine_bins-20+0.5 ) ;

         int hbi(1) ;
         for ( int fbi=0; fbi<n_fine_bins; fbi++ ) {

            if ( fb_mhti[fbi] == 4 && fb_hti[fbi] == 1 ) continue ;

            char hname3d[1000] ;
            sprintf( hname3d, "h_njvsmhtvsht_nb%d_%s_%s", fb_nbi[fbi], selname, sname[si] ) ;
            TH3F* h3p = (TH3F*) gDirectory -> FindObject( hname3d ) ;
            if ( h3p == 0x0 ) { printf("\n\n *** Missing histogram %s\n\n", hname3d ) ; return ; }
            int h3xbi = fb_hti[fbi]+1 ;
            int h3ybi = fb_mhti[fbi]+1 ;
            int h3zbi = fb_nji[fbi]+1 ;
            float val = h3p -> GetBinContent( h3xbi, h3ybi, h3zbi ) ;
            float err2 = pow( h3p -> GetBinError( h3xbi, h3ybi, h3zbi ), 2. ) ;

            char binlabel[1000] ;
            sprintf( binlabel, "%s %3d", fb_name[fbi], hbi ) ;
            h1p -> GetXaxis() -> SetBinLabel( hbi, binlabel ) ;

            h1p -> SetBinContent( hbi, val ) ;
            h1p -> SetBinError( hbi, sqrt(err2) ) ;

            hbi++ ;

         } // fbi


         printf("\n\n") ;
         h1p -> Print("all") ;
         printf("\n\n") ;

         h1p -> GetXaxis() -> LabelsOption( "v" ) ; ;


         h1d_samples[si] = h1p ;
         h1p -> SetFillColor( scolor[si] ) ;
         h_stack -> Add( h1p ) ;

         if ( si == 0 ) {
            sprintf( hname1d, "h_1d_allinone_%s_bgsum", selname ) ;
            h_bgsum = (TH1F*) h1p -> Clone( hname1d ) ;
            sprintf( htitle1d, "%s, bgsum", selname ) ;
            h_bgsum -> SetTitle( htitle1d ) ;
         } else {
            h_bgsum -> Add( h1p ) ;
         }

         if ( isll[si] ) {
            if ( h_llbgsum == 0x0 ) {
               sprintf( hname1d, "h_1d_allinone_%s_llbgsum", selname ) ;
               h_llbgsum = (TH1F*) h1p -> Clone( hname1d ) ;
               sprintf( htitle1d, "%s, llbgsum", selname ) ;
               h_llbgsum -> SetTitle( htitle1d ) ;
            } else {
               h_llbgsum -> Add( h1p ) ;
            }
         }

      } // si

      sprintf( hname, "h_frac_stack_%s", selname ) ;
      sprintf( htitle, "sample fractions, %s", selname ) ;
      THStack* h_frac_stack = new THStack( hname, htitle ) ;

      TH1F* h_frac_sum(0x0) ;

      for ( int si=0; si<n_samples; si++ ) {
         sprintf( hname, "%s_fr", h1d_samples[si] -> GetName() ) ;
         TH1F* h1d_fr = (TH1F*) h1d_samples[si] -> Clone( hname ) ;
         char htitle[1000] ;
         sprintf( htitle, "%s fraction", h1d_samples[si] -> GetTitle() ) ;
         h1d_fr -> SetTitle( htitle ) ;
         for ( int bi=1; bi<=h1d_samples[0] -> GetNbinsX(); bi++ ) {
            float sum_val = h_bgsum -> GetBinContent( bi ) ;
            if ( sum_val <= 0 ) continue ;
            float sum_err = h_bgsum -> GetBinError( bi ) ;
            float comp_val = h1d_samples[si] -> GetBinContent( bi ) ;
            float comp_err = h1d_samples[si] -> GetBinError( bi ) ;
            float fr_val = comp_val / sum_val ;
            float rem_val = sum_val - comp_val ;
            float rem_err = sqrt( pow(sum_err,2) - pow(comp_err,2) ) ;
            float fr_err = 0. ;
            if ( rem_val > 0. && comp_val > 0. ) fr_err = sqrt( pow( (rem_val*comp_err) / (sum_val*sum_val), 2 ) + pow( (comp_val*rem_err) / (sum_val*sum_val), 2 ) ) ;
            h1d_fr -> SetBinContent( bi, fr_val ) ;
            h1d_fr -> SetBinError( bi, fr_err ) ;
         } // bi
         h_frac_stack -> Add( h1d_fr ) ;
         if ( si == 0 ) {
            sprintf( hname, "h_bgfracsum_%s", selname ) ;
            h_frac_sum = (TH1F*) h1d_fr -> Clone( hname ) ;
            sprintf( htitle, "sample fractions, %s", selname ) ;
            h_frac_sum -> SetTitle( htitle ) ;
         } else {
            h_frac_sum -> Add( h1d_fr ) ;
         }
      } // si

      char fname[10000] ;

      h_bgsum -> SetMinimum( 0.1 ) ;
      h_frac_sum -> SetMaximum( 1.3 ) ;

      h_bgsum -> Draw( ) ;
      h_stack -> Draw( "hist same" ) ;
      h_bgsum -> Draw( "same" ) ;

      gPad -> SetLogy(0) ;
      sprintf( fname, "outputfiles/bgplot-1d-allinone-%s-liny.pdf", selname ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

      gPad -> SetLogy(1) ;
      sprintf( fname, "outputfiles/bgplot-1d-allinone-%s-logy.pdf", selname ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

      h_bgsum -> SetMaximum(25.) ;
      h_bgsum -> Draw( ) ;
      h_stack -> Draw( "hist same" ) ;
      h_bgsum -> Draw( "same" ) ;

      gPad -> SetLogy(0) ;
      sprintf( fname, "outputfiles/bgplot-1d-allinone-%s-zoom.pdf", selname ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;


      h_bgsum -> SetMaximum(5.) ;
      h_bgsum -> Draw( ) ;
      h_stack -> Draw( "hist same" ) ;
      h_bgsum -> Draw( "same" ) ;

      gPad -> SetLogy(0) ;
      sprintf( fname, "outputfiles/bgplot-1d-allinone-%s-zoom2.pdf", selname ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

      h_frac_sum -> Draw( "hist" ) ;
      h_frac_stack -> Draw( "hist same" ) ;
      h_frac_stack -> Draw( "same" ) ;
      sprintf( fname, "outputfiles/bgplot-1d-allinone-%s-frac.pdf", selname ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;


   } // draw_stack


  //====================================================================================================================


  //====================================================================================================================

   void make_legend( ) {

      char hname[1000] ;

      TH1F* h_legend_dummy = new TH1F( "h_legend_dummy", "", 2, 0., 1. ) ;

      TLegend* legend = new TLegend( 0.1, 0.2, 0.9, 0.9 ) ;

      for ( int si=n_samples-1; si>=0; si-- ) {

         sprintf( hname, "h_1d_allinone_zl_%s", sname[si] ) ;
         printf("   hname: %s\n", hname ) ;

         TH1F* h1d = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( h1d == 0x0 ) { printf("\n\n *** Missing histogram %s\n\n", hname ) ; return ; }

         legend -> AddEntry( h1d, sname[si] ) ;

      } // si

      h_legend_dummy -> Draw( "AH" ) ;
      legend -> Draw() ;

      char fname[10000] ;
      sprintf( fname, "outputfiles/bgplot-1d-allinone-legend.pdf" ) ;
      can -> Update() ; can -> Draw() ;
      can -> SaveAs( fname ) ;

   } // make_legend

  //====================================================================================================================

   void setup_binmap( const char* lhbuilder_file ) {

      ifstream ifs ;
      ifs.open( lhbuilder_file ) ;
      if ( !ifs.good() ) {
         printf("\n\n *** Problem with input file : %s\n\n", lhbuilder_file ) ; gSystem->Exit(0) ;
      }

      n_fine_bins = find_line_val( ifs, "N-fine-bins" ) ;
      printf("\n Number of fine bins is %d\n", n_fine_bins ) ;
      for ( int fbi=0; fbi<n_fine_bins; fbi++ ) {
         TString line ;
         line.ReadLine( ifs ) ;
         sscanf( line.Data(), "%s", fb_name[fbi] ) ;
         int nji, nbi, mhti, hti ;
         sscanf( fb_name[fbi], "FB-Njet%d-Nb%d-MHT%d-HT%d", &nji, &nbi, &mhti, &hti ) ;
         fb_nji[fbi] = nji ;
         fb_nbi[fbi] = nbi ;
         fb_mhti[fbi] = mhti ;
         fb_hti[fbi] = hti ;
         printf("  Fine bin %3d : %s,  nji=%d, nbi=%d, mhti=%d, hti=%d\n",
          fbi, fb_name[fbi], fb_nji[fbi], fb_nbi[fbi], fb_mhti[fbi], fb_hti[fbi] ) ;
      } // fbi.

      n_search_bins = find_line_val( ifs, "N-search-bins" ) ;
      printf("\n Number of search bins is %d\n", n_search_bins ) ;
      for ( int sbi=0; sbi<n_search_bins; sbi++ ) {
         TString line ;
         line.ReadLine( ifs ) ;
         sprintf( sb_name[sbi], "%s", line.Data() ) ;
      } // sbi

      if ( !find_line( ifs, "Fine-bin-search-bin-map" ) ) {
         printf("\n\n *** can't find map.\n\n") ; return ;
      }

      for ( int fbi=0; fbi<n_fine_bins; fbi++ ) {
         TString line ;
         line.ReadLine( ifs ) ;
         char fbname[100] ;
         char sbname[100] ;
         sscanf( line, "%s %s", fbname, sbname ) ;
         if ( strcmp( fbname, fb_name[fbi] ) != 0 ) { printf("\n\n *** Inconsistency %s %s\n", fbname, fb_name[fbi] ) ; return ; }
         if ( strcmp( sbname, "X" ) == 0 ) {
            fb_sbi[fbi] = -1 ;
         } else {
            for ( int sbi=0; sbi<n_search_bins; sbi++ ) {
               bool found(false) ;
               if ( strcmp( sbname, sb_name[sbi] ) == 0 ) {
                  fb_sbi[fbi] = sbi ;
                  found = true ;
                  break ;
               }
               if ( !found ) {
                  fb_sbi[fbi] = -1 ;
               }
            } // sbi
         }
         printf( " %3d %30s - %3d ", fbi, fb_name[fbi], fb_sbi[fbi] ) ;
         if ( fb_sbi[fbi] >= 0 ) {
            printf( "%30s\n", sb_name[fb_sbi[fbi]] ) ;
         } else {
            printf( "not used\n") ;
         }
      } // fbi.


   } // setup_binmap

  //====================================================================================================================

   bool find_line( ifstream& ifs, const char* key ) {

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
            return true ;
         }
      }

      return false ;


   } // find_line


  //====================================================================================================================

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

  //====================================================================================================================










