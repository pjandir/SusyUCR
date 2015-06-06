
#include "make_lhbuilder_input4b.c"

   void do_step3() {

      bool perfect_closure ;
      float true_sig_strength ;
      bool setup_qcdlhfit ;

      perfect_closure = false ;
      true_sig_strength = 0. ;
      setup_qcdlhfit = true ;
      make_lhbuilder_input4b( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists4-nbsum-t1bbbbH-postdraw.root", "t1bbbbH",
        "outputfiles/lhbuilder-input-v4b-nonperfect-closure-ss0-nbsum",
        setup_qcdlhfit ) ;

      perfect_closure = true ;
      true_sig_strength = 0. ;
      setup_qcdlhfit = true ;
      make_lhbuilder_input4b( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists4-nbsum-t1bbbbH-postdraw.root", "t1bbbbH",
        "outputfiles/lhbuilder-input-v4b-perfect-closure-ss0-nbsum",
        setup_qcdlhfit ) ;

      perfect_closure = false ;
      true_sig_strength = 1. ;
      setup_qcdlhfit = false ;
      make_lhbuilder_input4b( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists4-t1bbbbH-postdraw.root", "t1bbbbH",
        "outputfiles/lhbuilder-input-v4b-nonperfect-closure-ss1-fullfit",
        setup_qcdlhfit ) ;

      perfect_closure = true ;
      true_sig_strength = 1. ;
      setup_qcdlhfit = false ;
      make_lhbuilder_input4b( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists4-t1bbbbH-postdraw.root", "t1bbbbH",
        "outputfiles/lhbuilder-input-v4b-perfect-closure-ss1-fullfit",
        setup_qcdlhfit ) ;


   }

