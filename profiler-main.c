/******************************************************************************

 Command-line profiling utility

 (C) 2007 William Hart and David Harvey

******************************************************************************/

#include "flint.h"
#include "profiler-main.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


prof2d_Driver_t prof2d_active_Driver = NULL;
prof2d_DriverString_t prof2d_active_DriverString = NULL;
prof2d_Sampler_t prof2d_active_Sampler = NULL;


double do_single_run(prof2d_Sampler_t sampler, unsigned long x,
                     unsigned long y, unsigned long count)
{
   init_clock(0);
   start_clock(0);
   sampler(x, y, count);
   stop_clock(0);
   return get_clock(0);
}


// Timing runs need to last at least this many microseconds to be counted:
#define DURATION_THRESHOLD 200000
// Microseconds per timing run that the prof2d_sample function aims for:
#define DURATION_TARGET 300000


void prof2d_sample(unsigned long x, unsigned long y)
{
   // the number of timings that were at least DURATION_THRESHOLD microseconds:
   unsigned long good_count = 0;
   double max_time, min_time;

   // first try one loop
   unsigned long num_trials = 1;
   double last_time = do_single_run(prof2d_active_Sampler, x, y, 1);

   // loop until we have enough good times
   while (1)
   {
      double per_trial = last_time / num_trials;
      
      // if the last recorded time was long enough, record it
      if (last_time > DURATION_THRESHOLD)
      {
         if (good_count)
         {
            if (per_trial > max_time)
               max_time = per_trial;
            if (per_trial < min_time)
               min_time = per_trial;
         }
         else
            max_time = min_time = per_trial;

         if (++good_count == 5)
         {
            // we've got enough data
            break;
         }
      }

      // adjust num_trials so that the elapsed time gravitates towards
      // DURATION_TARGET; num_trials can be changed by a factor of
      // at most 25%, and must be at least 1
      if (last_time < 0.0001)
         last_time = 0.0001;
      double adjust_ratio = DURATION_TARGET / last_time;
      if (adjust_ratio > 1.25)
         adjust_ratio = 1.25;
      if (adjust_ratio < 0.75)
         adjust_ratio = 0.75;
      num_trials = (unsigned long) ceil(adjust_ratio * num_trials);
      // just to be safe:
      if (num_trials == 0)
         num_trials = 1;

      // run another trial
      last_time = do_single_run(prof2d_active_Sampler, x, y, num_trials);
   }

   // print results
   printf("%d %d %lf %lf\n", x, y, min_time, max_time);
}


void prof2d_set_sampler(prof2d_Sampler_t sampler)
{
   prof2d_active_Sampler = sampler;
}


void do_target(int index, int argc, char* argv[])
{
   printf("\n");
   printf("=========================\n");
   printf("print some header here which contains e.g. current timestamp,\n");
   printf("name of module, name of profile target, target description\n");
   printf("e.g. \"%s\", current machine name\n", prof2d_DriverString_list[index](argc, argv));
   printf("(from an environment variable probably), etc.\n");
   printf("=========================\n");

   prof2d_active_Driver = prof2d_Driver_list[index];

   if (prof2d_active_Driver != NULL)
   {
      prof2d_active_Driver(argc, argv);
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
      
      do_target(choice, argc, argv);

   }
   else
   {
      // if user supplies profile target on command, run it directly
      for (int i = 0; i < prof2d_target_count; i++)
      {
         if (!strcmp(prof2d_target_name[i], argv[1]))
         {
            do_target(i, argc, argv);
            return 0;
         }
      }
      printf("\n\nUnrecognised target name\n");
   }

   return 0;
}



// end of file ****************************************************************
