

void setup_lhbuilder_inputs( bool perfect_closure = true,
                             float true_sig_strength = 1.0,
                             const char* outfilebase = "outputfiles/lhbuilder-input" ) {

      gROOT -> LoadMacro("make_lhbuilder_input2.c") ;

    make_lhbuilder_input2( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists2-t1bbbbH-postdraw.root", "t1bbbbH", outfilebase ) ;
    make_lhbuilder_input2( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists2-t1bbbbC-postdraw.root", "t1bbbbC", outfilebase ) ;
    make_lhbuilder_input2( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists2-t1ttttH-postdraw.root", "t1ttttH", outfilebase ) ;
    make_lhbuilder_input2( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists2-t1ttttC-postdraw.root", "t1ttttC", outfilebase ) ;
    make_lhbuilder_input2( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists2-t1qqqqH-postdraw.root", "t1qqqqH", outfilebase ) ;
    make_lhbuilder_input2( perfect_closure, true_sig_strength, "outputfiles/fill-bg-hists2-t1qqqqC-postdraw.root", "t1qqqqC", outfilebase ) ;

} // setup_lhbuilder_inputs

