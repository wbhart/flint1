/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/

/*
Profiling MAGMA polynomial multiplication in Z[x].

Usage: run magma with the -b flag to prevent the start up banner, i.e.

   magma -b magma-profile.m > output.prof

(C) 2007 David Harvey + Bill Hart, GPL

*/

target_name := "PolyMul";
target_description := "MAGMA polynomial multiplication in Z[x] over various lengths and bitsizes, NON-NEGATIVE coefficients only";


max := 4000000;   // maximum total bitsize of input polys
ratio := 1.2;      // ratio between consecutive lengths/bitsizes


// Timing runs need to last at least this many microseconds to be counted:
DURATION_THRESHOLD := 200000;
// Microseconds per timing run that the prof2d_sample function aims for:
DURATION_TARGET := 300000;


forward prof2d_sample;

/*
This function should run count iterations at position (x, y),
and return the total time in seconds, using the Cputime() function.
*/
function sampler(length, bits, count)

    // first time random poly generation + multiplication

    countmod := 4;
    if count gt 1000 then countmod := 100; end if;
    if count gt 100 then countmod := 10; end if;
    time1 := Cputime();
    for i := 1 to count do
      if (i-1) mod countmod eq 0 then
         a:=Polynomial([RandomBits(bits): x in [1..length]]);
         b:=Polynomial([RandomBits(bits): x in [1..length]]);
      end if;
      c:=a*b;
    end for;

    time2 := Cputime();

    // now time just the random poly generation

    for i := 1 to count do
      if (i-1) mod countmod eq 0 then
       a:=Polynomial([RandomBits(bits): x in [1..length]]);
       b:=Polynomial([RandomBits(bits): x in [1..length]]);
       end if;
    end for;

    time3 := Cputime();
    return (time2 - time1) - (time3 - time2);
end function;


/*
This function should loop over appropriate combinations of (x, y),
and call prof2d_sample(x, y) for each one.
*/
procedure driver()

   max_iter := Ceiling(Log(max) / Log(ratio));

   last_length := 0;
   for i := 0 to max_iter do
      length := Floor(ratio^i);
      if length ne last_length then
         last_length := length;

         for bits := 1 to 31 do
            m := Ceiling(Log(length)) + 2*bits;
            if m le 62 then
                  prof2d_sample(length, bits);
            end if;
         end for;
      end if;
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
