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


prof_Driver_t prof_active_Driver = NULL;
prof_DriverString_t prof_active_DriverString = NULL;
prof_DriverDefaultParams_t prof_active_DriverDefaultParams = NULL;

prof1d_Sampler_t prof1d_active_Sampler = NULL;
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

void prof1d_set_sampler(prof1d_Sampler_t sampler)
{
   prof1d_active_Sampler = sampler;
}


void prof_start()
{
   start_clock(0);
}

void prof_stop()
{
   stop_clock(0);
}


typedef struct
{
   unsigned long x;
   void* arg;
   prof1d_Sampler_t sampler;
} prof1d_info_t;


typedef struct
{
   unsigned long x, y;
   void* arg;
   prof2d_Sampler_t sampler;
} prof2d_info_t;


void prof1d_info_exec(void* info, unsigned long count)
{
   prof1d_info_t* z = (prof1d_info_t*) info;
   
   z->sampler(z->x, z->arg, count);
}


void prof2d_info_exec(void* info, unsigned long count)
{
   prof2d_info_t* z = (prof2d_info_t*) info;
   
   z->sampler(z->x, z->y, z->arg, count);
}


void prof1d_sample(unsigned long x, void* arg)
{
   double min_time, max_time;
   prof1d_info_t info;
   
   info.x = x;
   info.arg = arg;
   info.sampler = prof1d_active_Sampler;
   
   prof_repeat(&min_time, &max_time, prof1d_info_exec, &info);

   printf("%d\t%.3le\t%.3le\n", x, min_time, max_time);
   fflush(stdout);
}


void prof2d_sample(unsigned long x, unsigned long y, void* arg)
{
   double min_time, max_time;
   prof2d_info_t info;
   
   info.x = x;
   info.y = y;
   info.arg = arg;
   info.sampler = prof2d_active_Sampler;
   
   prof_repeat(&min_time, &max_time, prof2d_info_exec, &info);

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
   printf("TARGET: %s\n", prof_target_name[index]);
   printf("PARAMETERS: %s\n", params);

   printf("\n");
   if (prof_DriverString_list[index])
      printf("DESCRIPTION:\n%s\n", prof_DriverString_list[index](params));
   
   printf("\n");
   printf("============================================== begin data \n");

   prof_active_Driver = prof_Driver_list[index];

   if (prof_active_Driver != NULL)
      prof_active_Driver(params);
}


// returns -1 if target name not found
int lookup_target_name(char* name)
{
   for (int i = 0; i < prof_target_count; i++)
   {
      if (!strcmp(prof_target_name[i], name))
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
   for (int i = 0; i < prof_target_count; i++)
      printf("  %s\n", prof_target_name[i]);
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

   // get default parameters from appropriate function if no
   // parameters supplied on command line
   if (!strlen(profile_params))
   {
      prof_DriverDefaultParams_t f =
                  prof_DriverDefaultParams_list[selected_target];
      if (f)
         strcpy(profile_params, f());
   }

   do_target(selected_target, profile_params);

   return 0;
}



// end of file ****************************************************************
