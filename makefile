ifdef PLEASE_DO_NOT_ERASE_DAVIDS_64_BIT_FLAG_AGAIN_THANKS_VERY_MUCH
	BITFLAG = -m64
else
	BITFLAG =
      TUNE = -mtune=opteron -march=opteron 
endif

CC = gcc -std=c99 $(BITFLAG)

ifndef FLINT_GMP_INCLUDE_DIR
	FLINT_GMP_INCLUDE_DIR = "/home/dmharvey/gmp/install/include"
endif

ifndef FLINT_GMP_LIB_DIR
	FLINT_GMP_LIB_DIR = "/home/dmharvey/gmp/install/lib"
endif

LIBS = -L$(FLINT_GMP_LIB_DIR) -lgmp -lpthread -lm
INCS =  -I"/usr/include" -I$(FLINT_GMP_INCLUDE_DIR)

CFLAGS = $(INCS) -funroll-loops -fexpensive-optimizations $(TUNE) -O3

RM = rm -f

.PHONY: all all-before all-after clean clean-custom driver test

all: all-before ssmul all-after

clean: clean-custom
	${RM} $(FLINTOBJ) $(PROFOBJ) $(TESTOBJ) $(DRIVEROBJ) $(BIN) $(FBIN)

HEADERS = Z.h Z_mpn.h ZmodF.h ZmodFpoly.h Zpoly.h fmpz_poly.h \
			 extras.h memory-manager.h flint.h longlong.h longlong_wrapper.h \
			 mpn_extras.h profiler-main.h profiler.h


mpn_extras.o: mpn_extras.c $(HEADERS)
	$(CC) -c mpn_extras.c -o mpn_extras.o $(CFLAGS)

Z.o: Z.c $(HEADERS)
	$(CC) -c Z.c -o Z.o $(CFLAGS)

memory-manager.o: memory-manager.c $(HEADERS)
	$(CC) -c memory-manager.c -o memory-manager.o $(CFLAGS)



profiler.o: profiler.c $(HEADERS)
	$(CC) -c profiler.c -o profiler.o $(CFLAGS)
	
profiler-main.o: profiler-main.c $(HEADERS)
	$(CC) -c profiler-main.c -o profiler-main.o $(CFLAGS)



Z_mpn.o: Z_mpn.c $(HEADERS)
	$(CC) -c Z_mpn.c -o Z_mpn.o $(CFLAGS)



Zpoly.o: Zpoly.c $(HEADERS)
	$(CC) -c Zpoly.c -o Zpoly.o $(CFLAGS)
	
Zpoly-test.o: Zpoly-test.c $(HEADERS)
	$(CC) -c Zpoly-test.c -o Zpoly-test.o $(CFLAGS)
	
Zpoly-test: Zpoly-test.o Zpoly.o memory-manager.o
	$(CC) Zpoly.o Zpoly-test.o memory-manager.o -o Zpoly-test $(CFLAGS) $(LIBS)
	


fmpz_poly.o: fmpz_poly.c $(HEADERS)
	$(CC) -c fmpz_poly.c -o fmpz_poly.o $(CFLAGS)
	
fmpz_poly-test.o: fmpz_poly-test.c $(HEADERS)
	$(CC) -c fmpz_poly-test.c -o fmpz_poly-test.o $(CFLAGS)
	
fmpz_poly-test: fmpz_poly-test.o ZmodFpoly.o ZmodF.o ZmodF_mul.o Z_mpn.o mpn_extras.o Zpoly.o fmpz_poly.o Zpoly.h memory-manager.o 
	$(CC) fmpz_poly.o ZmodFpoly.o ZmodF.o ZmodF_mul.o Z_mpn.o mpn_extras.o Zpoly.o fmpz_poly-test.o memory-manager.o -o fmpz_poly-test $(CFLAGS) $(LIBS)

fmpz_poly-profile-tables.o: fmpz_poly-profile.c $(HEADERS)
	python make-profile-tables.py fmpz_poly
	$(CC) -c fmpz_poly-profile-tables.c -o fmpz_poly-profile-tables.o $(CFLAGS)
	rm fmpz_poly-profile-tables.c
	
fmpz_poly-profile.o: fmpz_poly-profile.c $(HEADERS)
	$(CC) -c fmpz_poly-profile.c -o fmpz_poly-profile.o $(CFLAGS)

fmpz_poly-profile: Z_mpn.o ZmodFpoly.o memory-manager.o fmpz_poly-profile.o fmpz_poly-profile-tables.o fmpz_poly.o profiler-main.o profiler.o ZmodF.o ZmodF_mul.o fmpz_poly.o Zpoly.o mpn_extras.o
	$(CC) -o fmpz_poly-profile Z_mpn.o fmpz_poly-profile.o fmpz_poly-profile-tables.o profiler.o profiler-main.o fmpz_poly.o Zpoly.o memory-manager.o mpn_extras.o ZmodFpoly.o ZmodF.o ZmodF_mul.o $(CFLAGS) $(LIBS)



ZmodF.o: ZmodF.c $(HEADERS)
	$(CC) -c ZmodF.c -o ZmodF.o $(CFLAGS)
	
ZmodF_mul.o: ZmodF_mul.c $(HEADERS)
	$(CC) -c ZmodF_mul.c -o ZmodF_mul.o $(CFLAGS)

ZmodF-test.o: ZmodF-test.c $(HEADERS)
	$(CC) -c ZmodF-test.c -o ZmodF-test.o $(CFLAGS)

ZmodF-test: ZmodF-test.o ZmodF.o ZmodF_mul.o ZmodFpoly.o fmpz_poly.o memory-manager.o mpn_extras.o Z_mpn.o
	$(CC) ZmodF-test.o ZmodF.o ZmodF_mul.o ZmodFpoly.o fmpz_poly.o memory-manager.o mpn_extras.o Z_mpn.o -o ZmodF-test $(CFLAGS) $(LIBS)

ZmodF_mul-profile-tables.o: ZmodF_mul-profile.c $(HEADERS)
	python make-profile-tables.py ZmodF_mul
	$(CC) -c ZmodF_mul-profile-tables.c -o ZmodF_mul-profile-tables.o $(CFLAGS)
	rm ZmodF_mul-profile-tables.c

ZmodF_mul-profile.o: ZmodF_mul-profile.c $(HEADERS)
	$(CC) -c ZmodF_mul-profile.c -o ZmodF_mul-profile.o $(CFLAGS)

ZmodF_mul-profile: ZmodF_mul-profile.o ZmodF.o ZmodF_mul.o ZmodFpoly.o fmpz_poly.o memory-manager.o mpn_extras.o Z_mpn.o profiler-main.o ZmodF_mul-profile-tables.o profiler.o
	$(CC) -o ZmodF_mul-profile ZmodF_mul-profile.o ZmodF.o ZmodF_mul.o ZmodFpoly.o fmpz_poly.o memory-manager.o mpn_extras.o Z_mpn.o profiler-main.o ZmodF_mul-profile-tables.o profiler.o $(CFLAGS) $(LIBS)



ZmodFpoly.o: ZmodFpoly.c $(HEADERS)
	$(CC) -c ZmodFpoly.c -o ZmodFpoly.o $(CFLAGS)
	
ZmodFpoly-test.o: ZmodFpoly-test.c $(HEADERS)
	$(CC) -c ZmodFpoly-test.c -o ZmodFpoly-test.o $(CFLAGS)

ZmodFpoly-test: Z_mpn.o fmpz_poly.o Zpoly.o memory-manager.o mpn_extras.o ZmodFpoly-test.o ZmodFpoly.o ZmodF.o ZmodF_mul.o
	$(CC) Z_mpn.o fmpz_poly.o Zpoly.o memory-manager.o mpn_extras.o ZmodFpoly-test.o ZmodFpoly.o ZmodF.o ZmodF_mul.o -o ZmodFpoly-test $(CFLAGS) $(LIBS)

ZmodFpoly-profile-tables.o: ZmodFpoly-profile.c $(HEADERS)
	python make-profile-tables.py ZmodFpoly
	$(CC) -c ZmodFpoly-profile-tables.c -o ZmodFpoly-profile-tables.o $(CFLAGS)
	rm ZmodFpoly-profile-tables.c

ZmodFpoly-profile.o: ZmodFpoly-profile.c $(HEADERS)
	$(CC) -c ZmodFpoly-profile.c -o ZmodFpoly-profile.o $(CFLAGS)

ZmodFpoly-profile: Z_mpn.o memory-manager.o ZmodFpoly-profile.o ZmodFpoly-profile-tables.o ZmodFpoly.o profiler-main.o profiler.o ZmodF.o fmpz_poly.o Zpoly.o mpn_extras.o ZmodF_mul.o
	$(CC) -o ZmodFpoly-profile Z_mpn.o ZmodFpoly-profile.o ZmodFpoly-profile-tables.o profiler.o profiler-main.o fmpz_poly.o Zpoly.o memory-manager.o mpn_extras.o ZmodFpoly.o ZmodF.o ZmodF_mul.o $(CFLAGS) $(LIBS)
	

delta_qexp.o: delta_qexp.c $(HEADERS)
	$(CC) -c delta_qexp.c -o delta_qexp.o $(CFLAGS)

delta_qexp: delta_qexp.o Z_mpn.o fmpz_poly.o Zpoly.o memory-manager.o mpn_extras.o ZmodFpoly.o ZmodF.o ZmodF_mul.o
	$(CC) -o delta_qexp delta_qexp.o Z_mpn.o fmpz_poly.o Zpoly.o memory-manager.o mpn_extras.o ZmodFpoly.o ZmodF.o ZmodF_mul.o $(CFLAGS) $(LIBS)
