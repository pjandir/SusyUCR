


   void draw_per_bin_lh_results( const char* infile = "outputfiles/per-bin-lh-analysis-t1bbbbH.root",
                                 float ymax_signif = 2.5, float ymax_evts = 7. ) {

      gStyle -> SetPadBottomMargin( 0.45 ) ;
      gStyle -> SetPadRightMargin( 0.03 ) ;
      gStyle -> SetOptTitle(0) ;

      TCanvas* can = new TCanvas( "can_draw_per_bin_lh_results", "Per bin results", 1300, 1300 ) ;
      can -> Divide(1,2) ;

      TLine* line = new TLine() ;
      float ymax ;
      TText* ttext = new TText() ;
      ttext -> SetTextSize( 0.035 ) ;

      char text[1000] ;


      TFile f( infile ) ;

      f.ls() ;

    //-----------
      can -> cd(1) ;

      TH1F* h_signif = f.Get( "h_signif" ) ;
      if ( h_signif == 0x0 ) return ;
      h_signif -> SetMaximum( ymax_signif ) ;

      h_signif -> DrawCopy() ;
      gPad -> SetGridy(1) ;

      ymax = ymax_signif ;
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

      float sum2(0.) ;
      for ( int bi=1; bi<=h_signif->GetNbinsX(); bi++ ) {
         float val = h_signif -> GetBinContent( bi ) ;
         sum2 += val*val ;
      } //

      printf("\n\n Combined significance: %6.3f\n\n", sqrt(sum2) ) ;

      sprintf( text, "Combined significance: %6.3f", sqrt(sum2) ) ;
      ttext -> SetTextAlign(11) ;
      ttext -> SetTextColor( kRed ) ;
      ttext -> DrawTextNDC( 0.75, 0.95, text ) ;



    //-----------
      can -> cd(2) ;

      TH1F* h_sig = f.Get( "h_sig" ) ;
      if ( h_sig == 0x0 ) return ;
      h_sig -> SetMaximum( ymax_evts ) ;

      h_sig -> DrawCopy() ;
      gPad -> SetGridy(1) ;

      TGraph* g_bg = f.Get( "g_bg" ) ;
      g_bg -> Draw("P") ;

      TGraph* g_bg_stat_only = f.Get( "g_bg_stat_only" ) ;
      g_bg_stat_only -> Draw("P") ;

      ymax = ymax_evts ;
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


      TString outfile( infile ) ;
      outfile.ReplaceAll( ".root", "-results.pdf" ) ;
      can -> SaveAs( outfile ) ;


   }

