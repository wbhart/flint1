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
#include <time.h>
#include <gmp.h>


prof2d_Driver_t prof2d_active_Driver = NULL;
prof2d_DriverString_t prof2d_active_DriverString = NULL;
prof2d_Sampler_t prof2d_active_Sampler = NULL;


#define MACHINE_NAME_MAXLEN 1000
char machine_name[MACHINE_NAME_MAXLEN + 1];

#define PROFILE_PARAMS_MAXLEN 1000
char profile_params[PROFILE_PARAMS_MAXLEN + 1];


gmp_randstate_t profiler_main_randstate;


void profiler_random_limbs(unsigned long* output, unsigned long count)
{
   for (unsigned long i = 0; i < count; i++)
      output[i] = gmp_urandomb_ui(profiler_main_randstate,
                                  FLINT_BITS);
}


void prof2d_set_sampler(prof2d_Sampler_t sampler)
{
   prof2d_active_Sampler = sampler;
}


void prof2d_start()
{
   start_clock(0);
}

void prof2d_stop()
{
   stop_clock(0);
}


typedef struct
{
   unsigned long x, y;
   prof2d_Sampler_t sampler;
} prof2d_arg_t;


void prof2d_arg_exec(void* arg, unsigned long count)
{
   prof2d_arg_t* z = (prof2d_arg_t*) arg;
   
   z->sampler(z->x, z->y, count);
}


void prof2d_sample(unsigned long x, unsigned long y)
{
   double min_time, max_time;
   prof2d_arg_t arg;
   
   arg.x = x;
   arg.y = y;
   arg.sampler = prof2d_active_Sampler;
   
   prof_repeat(&min_time, &max_time, prof2d_arg_exec, &arg);

   printf("%d\t%d\t%.3le\t%.3le\n", x, y, min_time, max_time);
   fflush(stdout);
}


void do_target(int index, char* params)
{
   printf("FLINT profile output\n\n");

   time_t now;
   time(&now);
   printf("TIMESTAMP: %s", ctime(&now));
   printf("MACHINE: %s\n\n", machine_name);

   printf("MODULE: %s\n", prof_module_name);
   printf("TARGET: %s\n", prof2d_target_name[index]);
   printf("PARAMETERS: %s\n", strlen(params) ? params : "(none)");

   printf("\n");
   if (prof2d_DriverString_list[index])
      printf("DESCRIPTION:\n%s\n", prof2d_DriverString_list[index](params));
   
   printf("\n");
   printf("============================================== begin data \n");

   prof2d_active_Driver = prof2d_Driver_list[index];

   if (prof2d_active_Driver != NULL)
   {
      prof2d_active_Driver(params);
   }
}


// returns -1 if target name not found
int lookup_target_name(char* name)
{
   for (int i = 0; i < prof2d_target_count; i++)
   {
      if (!strcmp(prof2d_target_name[i], name))
         return i;
   }
   return -1;
}


void error(char* message)
{
   printf("Error: %s\n", message);
   exit(0);
}


void help()
{
   printf("FLINT profiling utility for module \"%s\"\n", prof_module_name);
   printf("\n");
   printf("options:\n");
   printf(" -h           Show this help screen\n");
   printf(" -t <target>  Target to run. Overrides environment variable\n");
   printf("              FLINT_PROFILE_TARGET.\n");
   printf(" -p <params>  Parameters to pass to target's Driver function.\n");
   printf("              Overrides environment variable FLINT_PROFILE_PARAMS.\n");
   printf("\n");
   printf("Targets in this profiling module are:\n");
   for (int i = 0; i < prof2d_target_count; i++)
      printf("  %s\n", prof2d_target_name[i]);
}


int main(int argc, char* argv[])
{
   gmp_randinit_default(profiler_main_randstate);

   // get name of current machine from environment variable
   char* machine_name_env = getenv("FLINT_MACHINE_NAME");
   if (machine_name_env)
      strncpy(machine_name, machine_name_env, MACHINE_NAME_MAXLEN);
   else
      strcpy(machine_name, "unknown");
      
   int selected_target = -1;
   char* profile_target_env = getenv("FLINT_PROFILE_TARGET");
   if (profile_target_env)
      selected_target = lookup_target_name(profile_target_env);
      
   char* profile_params_env = getenv("FLINT_PROFILE_PARAMS");
   if (profile_params_env)
      strncpy(profile_params, profile_params_env, PROFILE_PARAMS_MAXLEN);
   

   int num_args = 0;
   char** args = NULL;

   // scan command line options
   for (int n = 1; n < argc; n++)
   {
      if (!strcmp(argv[n], "-h"))
      {
         help();
         return 0;
      }
      else if (!strcmp(argv[n], "-t"))
      {
         // grab target name
         if (++n == argc)
            error("missing target name.");
         selected_target = lookup_target_name(argv[n]);
         if (selected_target == -1)
            error("unknown target name.");
      }
      else if (!strcmp(argv[n], "-p"))
      {
         if (++n == argc)
            error("missing target parameters.");
         strncpy(profile_params, argv[n], PROFILE_PARAMS_MAXLEN);
      }
      else
         error("unrecognised option.");
   }
   
   if (selected_target == -1)
   {
      help();
      return 0;
   }

   do_target(selected_target, profile_params);

   return 0;
}



// end of file ****************************************************************
