
#include "make_fakedata_kqcd_input_file1.c"

#include "TSystem.h"


   void make_toy_fakedata_kqcd_input_files1(
          int ntoys = 100,
          int first_seed = 1,
          bool perfect_qcd_closure = true,
          const char* kqcd_output_dir = "outputfiles/kqcd-fakedata-toys",
          const char* kqcd_input_lostlep_file = "outputfiles/kqcd-input-lostlep.txt",
          const char* kqcd_input_hadtau_file = "outputfiles/kqcd-input-hadtau.txt",
          const char* kqcd_input_znunu_file = "outputfiles/kqcd-input-znunu.txt",
          const char* kqcd_input_qcd_file = "outputfiles/kqcd-input-qcdmc.txt",
          const char* kqcd_fitconfig_file = "outputfiles/kqcd-fitconfig-with-constraints1.txt"
          ) {

      char command[10000] ;


      sprintf( command, "mkdir -p %s", kqcd_output_dir ) ;
      gSystem -> Exec( command ) ;

      sprintf( command, "cp -p %s %s/", kqcd_fitconfig_file, kqcd_output_dir ) ;
      gSystem -> Exec( command ) ;
      sprintf( command, "cp -p %s %s/", kqcd_input_lostlep_file, kqcd_output_dir ) ;
      gSystem -> Exec( command ) ;
      sprintf( command, "cp -p %s %s/", kqcd_input_hadtau_file, kqcd_output_dir ) ;
      gSystem -> Exec( command ) ;
      sprintf( command, "cp -p %s %s/", kqcd_input_znunu_file, kqcd_output_dir ) ;
      gSystem -> Exec( command ) ;
      sprintf( command, "cp -p %s %s/", kqcd_input_qcd_file, kqcd_output_dir ) ;
      gSystem -> Exec( command ) ;

      FILE* ofp ;
      char output_config_file[10000] ;
      sprintf( output_config_file, "%s/make-fakedata-config.txt", kqcd_output_dir ) ;
      if ( (ofp=fopen( output_config_file, "w" ))==NULL ) {
         printf("\n\n *** Can't open output file %s\n\n", output_config_file ) ;
      }
      fprintf( ofp, "Ntoys %d\n", ntoys ) ;
      fprintf( ofp, "first seed %d\n", first_seed ) ;
      if ( perfect_qcd_closure ) {
         fprintf( ofp, "perfect closure = TRUE\n") ;
      } else {
         fprintf( ofp, "perfect closure = FALSE\n") ;
      }
      fclose(ofp) ;


      for ( int ti=1; ti<=ntoys; ti++ ) {

         bool randomize_nobs(true) ;
         int  random_seed = first_seed + ti-1 ;
         bool make_plots(false) ;

         char kqcd_output_fakedata_file[10000] ;
         sprintf( kqcd_output_fakedata_file, "%s/kqcd-input-fakedata-toy%04d.txt", kqcd_output_dir, ti ) ;
         make_fakedata_kqcd_input_file1( perfect_qcd_closure, randomize_nobs, random_seed, make_plots,
            kqcd_output_fakedata_file,
            kqcd_input_lostlep_file,
            kqcd_input_hadtau_file,
            kqcd_input_znunu_file,
            kqcd_input_qcd_file,
            kqcd_fitconfig_file ) ;


      } // ti


   } // make_toy_fakedata_kqcd_input_files1


