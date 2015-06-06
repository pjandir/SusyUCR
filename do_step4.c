
//
//  Do NOT load this one with the compile option.  That is, do this
//      .L do_step4.c
//  not this
//      .L do_step4.c+
//

#include "wrapper_build_ra2b_ws4b.c"

   void do_step4() {

      wrapper_build_ra2b_ws4b( "outputfiles/lhbuilder-input-v4b-nonperfect-closure-ss0-nbsum-t1bbbbH-qcdlhfit.txt",
                               "outputfiles/ws-v4b-t1bbbbH-nonperfect-closure-nbsum-qcdlhfit.root",
                               0., true, true, 0., true ) ;

      wrapper_build_ra2b_ws4b( "outputfiles/lhbuilder-input-v4b-perfect-closure-ss0-nbsum-t1bbbbH-qcdlhfit.txt",
                               "outputfiles/ws-v4b-t1bbbbH-perfect-closure-nbsum-qcdlhfit.root",
                               0., true, true, 0., true ) ;

      wrapper_build_ra2b_ws4b( "outputfiles/lhbuilder-input-v4b-nonperfect-closure-ss1-fullfit-t1bbbbH.txt",
                               "outputfiles/ws-v4b-t1bbbbH-nonperfect-closure-ss1-fullfit.root",
                               0., true, true, 0., true ) ;

      wrapper_build_ra2b_ws4b( "outputfiles/lhbuilder-input-v4b-perfect-closure-ss1-fullfit-t1bbbbH.txt",
                               "outputfiles/ws-v4b-t1bbbbH-perfect-closure-ss1-fullfit.root",
                               0., true, true, 0., true ) ;

   }

