/*
Script for 2d profiles of MAGMA.

This file is part of the FLINT project. It generates output in the same
format as profiler-main.c.

Please make a COPY of this file before using it; don't overwrite this template
file in the FLINT repository!

Usage: run magma with the -b flag to prevent the start up banner, i.e.

   magma -b magma-profile.m > output.prof

(C) 2007 David Harvey, GPL

*/

/************************************************************************

 The first section of the file is a template you need to fill in to
 run a particular profile.

************************************************************************/

target_name := "target name here";
target_description := "some description here";


// Timing runs need to last at least this many microseconds to be counted:
DURATION_THRESHOLD := 200000;
// Microseconds per timing run that the prof2d_sample function aims for:
DURATION_TARGET := 300000;
    

forward prof2d_sample;

/*
This function should run "count" iterations at position (x, y),
and return the total time in seconds, using the Cputime() function.
*/
function sampler(x, y, count)

   // add setup code here

   time1 := Cputime();
   
   // add stuff to be timed here

   time2 := Cputime();

   return time2 - time1;
end function;


/*
This function should loop over appropriate combinations of (x, y),
and call prof2d_sample(x, y) for each one.
*/
procedure driver()

   // here is an example that calls prof2d_sample for a range of x and y
   
   for x := 1 to 10 do
      for y := 1 to 10 do
         prof2d_sample(x, y);
      end for;
   end for;
    
end procedure;



/************************************************************************

 This last section is the generic profiling code. Just leave this
 stuff alone.

************************************************************************/

/*
Formats in scientific notation with 3 decimal places
*/
function format_sci(x)
   L := Floor(Log(10, x));
   x := x / 10^L;
   s := Sprintf("%.3oe", x);
   if L lt 0 then
   s := s cat "-";
   else
   s := s cat "+";
   end if;

   s := s cat Sprintf("%o", Floor(Abs(L / 10)));
   s := s cat Sprintf("%o", (Abs(L) mod 10));

   return s;
end function;


procedure prof2d_sample(x, y)
   // number of timings that were at least DURATION_THRESHOLD microseconds:
   good_count := 0;

   // first try just a few loops
   num_trials := 4;
   last_time := sampler(x, y, num_trials) * 1000000.0;

   max_time := 0;
   min_time := 0;

   // loop until we have enough good times
   while true do
      per_trial := last_time / num_trials;

      // if the last recorded time was long enough, record it
      if last_time gt DURATION_THRESHOLD then
         if good_count gt 0 then
            max_time := Max(max_time, per_trial);
            min_time := Min(min_time, per_trial);
         else
            max_time := per_trial;
            min_time := per_trial;
         end if;

         good_count := good_count + 1;
         if good_count eq 5 then
            // we've got enough data
            // print it out and return
            print Sprintf("%o\t%o\t%o\t%o", x, y, format_sci(min_time), format_sci(max_time));
            return;
         end if;
      end if;

      // adjust num_trials so that the elapsed time gravitates towards
      // DURATION_TARGET; num_trials can be changed by a factor of
      // at most 25%, and must be at least 1
      if last_time lt 0.0001 then
         last_time := 0.0001;
      end if;
      adjust_ratio := 1.0 * DURATION_TARGET / last_time;
      if adjust_ratio gt 1.25 then
         adjust_ratio := 1.25;
      end if;
      if adjust_ratio lt 0.75 then
         adjust_ratio := 0.75;
      end if;
      num_trials := Ceiling(adjust_ratio * num_trials);
      // just to be safe:
      if num_trials eq 0 then
         num_trials := 1;
      end if;

      // run another trial
      last_time := sampler(x, y, num_trials) * 1000000.0;
   end while;
end procedure;


procedure print_header()
    print "FLINT profile output";
    print "";
    print "TIMESTAMP: (todo: write code to generate timestamp)";
    print "MACHINE: (todo: write code to get machine from environment var)";

    print "";
    print "MODULE: magma";
    print "TARGET:", target_name;
    print "";
    print "DESCRIPTION:";
    print target_description;
    print "";
    print "============================================== begin data";
    
end procedure;


print_header();
driver();

quit;

// ------------- end of file ------------------------------------
