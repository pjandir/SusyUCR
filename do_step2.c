

#include "draw_bgsig_hists4.c"

void do_step2() {

   draw_bgsig_hists4( "outputfiles/fill-bg-hists4-nbsum.root", "outputfiles/fill-sig-hists4-nbsum.root", "t1bbbbH" ) ;

   gDirectory -> Delete( "h*" ) ;

   draw_bgsig_hists4( "outputfiles/fill-bg-hists4.root", "outputfiles/fill-sig-hists4.root", "t1bbbbH" ) ;

}




