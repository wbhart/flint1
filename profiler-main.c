/******************************************************************************

 Command-line profiling utility

 (C) 2007 William Hart and David Harvey

******************************************************************************/

#include "flint.h"
#include "profiler-main.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


prof2d_main_t prof2d_active_main = NULL;
prof2d_exec_t prof2d_active_exec = NULL;


void prof2d_exec(unsigned long x, unsigned long y)
{
   // eventually this will try a bunch of different iteration counts,
   // let's just run it 100 times for now
   
   init_clock(0);
   start_clock(0);
   prof2d_active_exec(x, y, 100);
   stop_clock(0);
   printf("%d %d %lf\n", x, y, get_clock(0) / 100);
}


void do_target(int index)
{
   printf("\n");
   printf("=========================\n");
   printf("print some header here which contains e.g. current timestamp,\n");
   printf("name of module, name of profile target, target description\n");
   printf("e.g. \"%s\", current machine name\n", prof2d_target_string[index]);
   printf("(from an environment variable probably), etc.\n");
   printf("=========================\n");

   prof2d_active_main = prof2d_target_main[index];
   prof2d_active_exec = prof2d_target_exec[index];

   if (prof2d_active_main != NULL)
   {
      prof2d_active_main();
   }
}


int main(int argc, char* argv[])
{
   if (argc == 1)
   {
      // no args supplied; work in interactive mode
      printf("FLINT profiling utility for module \"%s\"\n\n", prof_module_name);
      
      for (int i = 0; i < prof2d_target_count; i++)
         printf("[ %2d ]: %s\n", i, prof2d_target_name[i]);
      
      printf("\n");
      printf("Select one of the above: ");
      
      unsigned choice;
      scanf("%d", &choice);
      
      if (choice >= prof2d_target_count)
      {
         printf("\n\nUmm yeah try again.\n");
         return 0;
      }
      
      do_target(choice);

   }
   else
   {
      // if user supplies profile target on command, run it directly
      for (int i = 0; i < prof2d_target_count; i++)
      {
         if (!strcmp(prof2d_target_name[i], argv[1]))
         {
            do_target(i);
            return 0;
         }
      }
      printf("\n\nUnrecognised target name\n");
   }

   return 0;
}



// end of file ****************************************************************
