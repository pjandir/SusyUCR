
#include "run_qcdlhfit1.c"

#include <fstream>

#include "TSystem.h"

   void run_toy_qcdlhfits1( const char* io_dir = "outputfiles/kqcd-fakedata-toys",
                       float fixed_sig_strength = 0.,
                       bool make_all_plots = false
   ) {

      char command[10000] ;

      printf("\n\n Making list of toy files.\n") ;
      sprintf( command, "ls -1 %s/ws*.root > toywsfiles.txt", io_dir ) ;
      gSystem -> Exec( command ) ;

      ifstream ifs ;
      ifs.open( "toywsfiles.txt" ) ;
      int fi(0) ;

      while ( ifs.good() ) {

         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( !ifs.good() ) break ;
         printf(" file %4d : %s\n", fi, ts.Data() ) ;

         printf("  %s \n", ts.Data() ) ;

         run_qcdlhfit1( ts.Data(), fixed_sig_strength, make_all_plots ) ;

         fi ++ ;
      }


   } // run_toy_qcdlhfits1

