
#include "wrapper_build_qcdlhfit_ws1.c"

#include "TSystem.h"

   void make_toy_fakedata_kqcd_ws_files1( const char* io_dir = "outputfiles/kqcd-fakedata-toys",
                            const char* fname_fitconfig = "outputfiles/kqcd-fitconfig-with-constraints1.txt",
                            const char* fname_lostlep = "outputfiles/kqcd-input-lostlep.txt",
                            const char* fname_hadtau  = "outputfiles/kqcd-input-hadtau.txt",
                            const char* fname_znunu   = "outputfiles/kqcd-input-znunu.txt"
         ) {

      bool  skip_testfit = true ;
      bool  skip_modelconfig = true ;

      char command[10000] ;

      printf("\n\n Making list of toy files.\n") ;
      sprintf( command, "ls -1 %s/*data*.txt > toyfiles.txt", io_dir ) ;
      gSystem -> Exec( command ) ;

      ifstream ifs ;
      ifs.open( "toyfiles.txt" ) ;
      int fi(0) ;

      while ( ifs.good() ) {

         TString ts ;
         ts.ReadLine( ifs ) ;
         if ( !ifs.good() ) break ;
         printf(" file %4d : %s\n", fi, ts.Data() ) ;

         TString outfile = ts ;
         outfile.ReplaceAll( "kqcd-input-fakedata", "ws" ) ;
         outfile.ReplaceAll( ".txt", ".root" ) ;

         printf("  %s --- %s\n", ts.Data(), outfile.Data() ) ;

         wrapper_build_qcdlhfit_ws1( outfile.Data(),
            fname_fitconfig,
            ts.Data(),
            skip_testfit,
            skip_modelconfig,
            fname_lostlep,
            fname_hadtau,
            fname_znunu
            ) ;

         fi ++ ;
      }



   } // make_toy_fakedata_kqcd_ws_files1


