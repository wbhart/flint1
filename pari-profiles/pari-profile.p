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

/***************************************************************** 

pari-profile.p - code for timing PARI polynomial multiplication in Z[x]
 
over various lengths and bitsizes. Based on Magma profiling code 
by David Harvey.

Copyright (C) 2007, Tomasz Lechowski

Some corrections (C) 2007, Bill Allombert


*****************************************************************/

target_name="PolyMul";
target_description="PARI polynomial multiplication in Z[x] over various
 lengths and bitsizes, NON-NEGATIVE coefficients only";

Max=16000000;
ratio=1.2;
\p4;

DURATION_THRESHOLD=200000;
DURATION_TARGET=300000;


sampler(Length, bits, count)={

countmod=4;
if(count>1000,countmod=100);
if(count>100,countmod=10);
gettime();

for(i=1,count,
	if(((i-1)%countmod)==0,
	a=vector(Length);
	b=vector(Length);
		for(j=1,Length,
		a[j]=random(2^(bits));
		b[j]=random(2^(bits)));
        a=Pol(a);
        b=Pol(b));
c=a*b);

time1=gettime();

for(i=1,count,
	if(((i-1)%countmod)==0,
	a=vector(Length); 	b=vector(Length);
		for(j=1,Length,
		a[j]=random(2^(bits));
		b[j]=random(2^(bits)));
        a=Pol(a);
        b=Pol(b));
);

time2=gettime();

return((time1-time2)/1000);


};






format_sci(x)={
    L=floor(log(x)/log(10));
    x=x/(10^L);
    s=Str(x);
    
    if(L<0,
    s=concat(s,"e-"),
    s=concat(s,"e+")
    );

    s=concat(s,Str(floor(abs(L / 10))));
    s=concat(s,Str((abs(L))%10));

    return(s);
};

prof2d_sample(x, y)={

    good_count=0;

    num_trials=4;
    last_time=sampler(x,y,num_trials)*1000000.0;

    max_time=0;
    min_time=0;

    true=1;
    while(true,
      per_trial=last_time/num_trials;
      if(last_time>DURATION_THRESHOLD,
          if(good_count>0,
            max_time=max(max_time,per_trial);
            min_time=min(min_time, per_trial),
            max_time=per_trial;
            min_time=per_trial);

          good_count=good_count+1;
          if(good_count==5,
            print(Str(x,"	", y,"	",format_sci(min_time),"	",format_sci(max_time)));
            write(Pariprof,Str(x,"	", y,"	",format_sci(min_time),"	",format_sci(max_time)));
	    return;
          )
      );

     
      if(last_time<0.0001,
          last_time=0.0001
      );
      adjust_ratio= 1.0*DURATION_TARGET/last_time;
      if(adjust_ratio>1.25,
          adjust_ratio=1.25
      );
      if(adjust_ratio<0.75,
          adjust_ratio=0.75
      );
      num_trials=ceil(adjust_ratio*num_trials);
      if(num_trials==0,
          num_trials=1
      );

      last_time=sampler(x,y,num_trials)*1000000.0;
    )
};

print_header()={
    print("FLINT profile output");
    write(Pariprof,"FLINT profile output");
    print("");
    write(Pariprof,"   ");
    print("TIMESTAMP: ");
    write(Pariprof,"TIMESTAMP:");
    print("MACHINE: ");
    write(Pariprof,"MACHINE:");

    print("");
    write(Pariprof,"   ");
    print("MODULE:PARI");
    write(Pariprof,"MODULE:PARI");
    print("TARGET:", target_name);
    write(Pariprof, "TARGET:", Str(target_name));
    print("");
    write(Pariprof,"  ");
    print("DESCRIPTION:");
    write(Pariprof,"DESCRIPTION:");
    print(target_description);
    write(Pariprof, Str(target_description));
    print("");
    write(Pariprof,"  ");
    print("============================================== begin data");
    write(Pariprof,"============================================== begin data");
   
};

driver()={

max_iter=ceil(log(Max)/log(ratio));

last_length=0;
for(i=0,max_iter,
	Length=floor((ratio)^i);
	if(Length!=last_length,
		last_length=Length;
		last_bits=0;
		for(j=0,max_iter,
			bits=floor(ratio^j);
			if(bits!=last_bits,
				last_bits=bits;
				if(Length*bits<Max,
					prof2d_sample(Length, bits);
				)
			)
		)
	)
)

};       



print_header();
driver()
