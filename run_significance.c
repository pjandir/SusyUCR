

void run_significance( bool check_signif_with_all_bins = false ) {

   gROOT -> LoadMacro( "significance_by_sb.c+" ) ;

   significance_by_sb( "outputfiles/ws-t1bbbbH.root", "outputfiles/significance-per-bin-t1bbbbH.pdf" ) ;
   significance_by_sb( "outputfiles/ws-t1bbbbC.root", "outputfiles/significance-per-bin-t1bbbbC.pdf" ) ;
   significance_by_sb( "outputfiles/ws-t1ttttH.root", "outputfiles/significance-per-bin-t1ttttH.pdf" ) ;
   significance_by_sb( "outputfiles/ws-t1ttttC.root", "outputfiles/significance-per-bin-t1ttttC.pdf" ) ;
   significance_by_sb( "outputfiles/ws-t1qqqqH.root", "outputfiles/significance-per-bin-t1qqqqH.pdf" ) ;
   significance_by_sb( "outputfiles/ws-t1qqqqC.root", "outputfiles/significance-per-bin-t1qqqqC.pdf" ) ;

}



