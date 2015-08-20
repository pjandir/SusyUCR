
#include "fill_sig_hists4.c"

void do_step1b() {

   bool sum_over_nb ;

   sum_over_nb = false ;
   fill_sig_hists4( 3000, sum_over_nb ) ;

   gDirectory -> Delete("h*") ;

   sum_over_nb = true ;
   fill_sig_hists4( 3000, sum_over_nb ) ;

}

