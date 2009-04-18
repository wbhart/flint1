/*
   tune_main.c:  main() routine for tuning program
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "support.h"
#include "zn_poly_internal.h"
#include "profiler.h"


char* header = 
"/*\n"
"   NOTE: do not edit this file! It is auto-generated by the \"tune\" program.\n"
"   (Run \"make tune\" and then \"./tune > tuning.c\" to regenerate it.)\n"
"*/\n"
"\n"
"/*\n"
"   tuning.c:  global tuning values\n"
"\n"
"   Copyright (C) 2007, 2008, David Harvey\n"
"\n"
"   This file is part of the zn_poly library (version 0.9).\n"
"\n"
"   This program is free software: you can redistribute it and/or modify\n"
"   it under the terms of the GNU General Public License as published by\n"
"   the Free Software Foundation, either version 2 of the License, or\n"
"   (at your option) version 3 of the License.\n"
"\n"
"   This program is distributed in the hope that it will be useful,\n"
"   but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"   GNU General Public License for more details.\n"
"\n"
"   You should have received a copy of the GNU General Public License\n"
"   along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
"\n"
"*/\n"
"\n"
"#include \"zn_poly_internal.h\"\n"
"\n";


char* footer =
"};\n"
"\n"
"// end of file ****************************************************************\n";


int
main (int argc, char* argv[])
{
   fprintf (stderr, "zn_poly tuning program\n");
   fprintf (stderr, "(use -v flag for verbose output)\n\n");

   gmp_randinit_default (randstate);

#if ZNP_HAVE_CYCLE_COUNTER

   calibrate_cycle_scale_factor ();

   int verbose = 0;

   if (argc == 2 && !strcmp (argv[1], "-v"))
      verbose = 1;

   // call various routines to generate tuning data
   tune_mpn_smp_kara (stderr, verbose);
   tune_mpn_mulmid_fallback (stderr, verbose);
   tune_mul_KS (stderr, 0, verbose);
   tune_mul_KS (stderr, 1, verbose);
   tune_mulmid_KS (stderr, verbose);
   tune_nuss (stderr, 0, verbose);
   tune_nuss (stderr, 1, verbose);
   tune_mul (stderr, 0, verbose);
   tune_mul (stderr, 1, verbose);
   tune_mulmid (stderr, verbose);

#else

   fprintf (stderr, "\n"
            "Warning: no cycle counting code available on this system,\n"
            "using default tuning values.\n\n");

#endif
   
   unsigned bits;
   size_t x;

   // generate tuning.c file
   printf (header);

   x = ZNP_mpn_smp_kara_thresh;
   printf ("size_t ZNP_mpn_smp_kara_thresh = ");
   printf (x == SIZE_MAX ? "SIZE_MAX;\n" : "%lu;\n", x);

   x = ZNP_mpn_mulmid_fallback_thresh;
   printf ("size_t ZNP_mpn_mulmid_fallback_thresh = ");
   printf (x == SIZE_MAX ? "SIZE_MAX;\n" : "%lu;\n", x);

   printf ("\n");

   printf ("tuning_info_t tuning_info[] = \n");
   printf ("{\n");
   printf ("   {  // bits = 0\n");
   printf ("   },\n");
   printf ("   {  // bits = 1\n");
   printf ("   },\n");
   
   for (bits = 2; bits <= ULONG_BITS; bits++)
   {
      printf ("   {  // bits = %u\n", bits);
      
      x = tuning_info[bits].mul_KS2_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS1 -> KS2 multiplication threshold\n");

      x = tuning_info[bits].mul_KS4_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS2 -> KS4 multiplication threshold\n");

      x = tuning_info[bits].mul_fft_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS4 -> FFT multiplication threshold\n");
      
      x = tuning_info[bits].sqr_KS2_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS1 -> KS2 squaring threshold\n");

      x = tuning_info[bits].sqr_KS4_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS2 -> KS4 squaring threshold\n");

      x = tuning_info[bits].sqr_fft_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS4 -> FFT squaring threshold\n");
      
      x = tuning_info[bits].mulmid_KS2_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS1 -> KS2 middle product threshold\n");

      x = tuning_info[bits].mulmid_KS4_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS2 -> KS4 middle product threshold\n");

      x = tuning_info[bits].mulmid_fft_thresh;
      printf (x == SIZE_MAX ? "   SIZE_MAX," : "      %5lu,", x);
      printf ("   // KS4 -> FFT middle product threshold\n");
      
      printf ("      %5lu,   // nussbaumer multiplication threshold\n",
              tuning_info[bits].nuss_mul_thresh);
      printf ("      %5lu    // nussbaumer squaring threshold\n",
              tuning_info[bits].nuss_sqr_thresh);

      printf ("   },\n");
   }
   
   printf (footer);
   
   gmp_randclear (randstate);
   
   return 0;
}


// end of file ****************************************************************
