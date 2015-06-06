

void run_significance( bool check_signif_with_all_bins = false, float ymax = 3.5 ) {


   gROOT -> LoadMacro( "significance_by_sb.c+" ) ;

   bool fix_nuisance_pars ;
   bool fix_bg_mu_pars ;

   fix_nuisance_pars = true ;
   fix_bg_mu_pars = true ;
   significance_by_sb( "outputfiles/ws-t1bbbbH.root", "outputfiles/significance-per-bin-t1bbbbH-NPfixed-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1bbbbC.root", "outputfiles/significance-per-bin-t1bbbbC-NPfixed-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1ttttH.root", "outputfiles/significance-per-bin-t1ttttH-NPfixed-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1ttttC.root", "outputfiles/significance-per-bin-t1ttttC-NPfixed-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1qqqqH.root", "outputfiles/significance-per-bin-t1qqqqH-NPfixed-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1qqqqC.root", "outputfiles/significance-per-bin-t1qqqqC-NPfixed-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;

// fix_nuisance_pars = false ;
// fix_bg_mu_pars = true ;
// significance_by_sb( "outputfiles/ws-t1bbbbH.root", "outputfiles/significance-per-bin-t1bbbbH-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1bbbbC.root", "outputfiles/significance-per-bin-t1bbbbC-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1ttttH.root", "outputfiles/significance-per-bin-t1ttttH-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1ttttC.root", "outputfiles/significance-per-bin-t1ttttC-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1qqqqH.root", "outputfiles/significance-per-bin-t1qqqqH-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1qqqqC.root", "outputfiles/significance-per-bin-t1qqqqC-bgMufixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;

// fix_nuisance_pars = true ;
// fix_bg_mu_pars = false ;
// significance_by_sb( "outputfiles/ws-t1bbbbH.root", "outputfiles/significance-per-bin-t1bbbbH-NPfixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1bbbbC.root", "outputfiles/significance-per-bin-t1bbbbC-NPfixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1ttttH.root", "outputfiles/significance-per-bin-t1ttttH-NPfixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1ttttC.root", "outputfiles/significance-per-bin-t1ttttC-NPfixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1qqqqH.root", "outputfiles/significance-per-bin-t1qqqqH-NPfixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
// significance_by_sb( "outputfiles/ws-t1qqqqC.root", "outputfiles/significance-per-bin-t1qqqqC-NPfixed.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;

   fix_nuisance_pars = false ;
   fix_bg_mu_pars = false ;
   significance_by_sb( "outputfiles/ws-t1bbbbH.root", "outputfiles/significance-per-bin-t1bbbbH.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1bbbbC.root", "outputfiles/significance-per-bin-t1bbbbC.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1ttttH.root", "outputfiles/significance-per-bin-t1ttttH.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1ttttC.root", "outputfiles/significance-per-bin-t1ttttC.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1qqqqH.root", "outputfiles/significance-per-bin-t1qqqqH.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;
   significance_by_sb( "outputfiles/ws-t1qqqqC.root", "outputfiles/significance-per-bin-t1qqqqC.pdf", ymax, check_signif_with_all_bins, fix_nuisance_pars, fix_bg_mu_pars ) ;

}



