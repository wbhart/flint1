ifndef FLINT_TUNE
	# defaults for sage.math development machine
	FLINT_TUNE = -mtune=opteron -march=opteron

	# for the record, here's what I use on my G5 powerpc:
	# FLINT_TUNE = -m64 -mcpu=970 -mtune=970 -mpowerpc64 -falign-loops=16 -falign-functions=16 -falign-labels=16 -falign-jumps=16

	# and here's for my laptop:
	# FLINT_TUNE = 
endif

ifndef FLINT_LINK_OPTIONS
	FLINT_LINK_OPTIONS = -static
endif

# default GMP directories on sage.math development machine
ifndef FLINT_GMP_INCLUDE_DIR
	FLINT_GMP_INCLUDE_DIR = "/home/dmharvey/gmp/install/include"
endif

ifndef FLINT_GMP_LIB_DIR
	FLINT_GMP_LIB_DIR = "/home/dmharvey/gmp/install/lib"
endif

LIBS = -L$(FLINT_GMP_LIB_DIR) $(FLINT_LINK_OPTIONS) -lgmp -lpthread -lm
INCS =  -I"/usr/include" -I$(FLINT_GMP_INCLUDE_DIR)

CC = gcc -std=c99

CFLAGS = $(INCS) -funroll-loops -fexpensive-optimizations $(FLINT_TUNE) -O3

RM = rm -f

HEADERS = \
	Z.h \
	Z_mpn.h \
	Z_mpn_mul-tuning.h \
	ZmodF.h \
	ZmodF_mul-tuning.h \
	ZmodF_mul.h \
	ZmodF_poly.h \
	extras.h \
	flint.h \
	fmpz.h \
	fmpz_poly.h \
	longlong.h \
	longlong_wrapper.h \
	memory-manager.h \
	mpn_extras.h \
	mpz_poly.h \
	profiler-main.h \
	profiler.h \
	test-support.h


####### library object files

FLINTOBJ = \
	mpn_extras.o \
	Z.o \
	memory-manager.o \
	Z_mpn.o \
	ZmodF.o \
	ZmodF_mul.o \
	fmpz.o \
	fmpz_poly.o \
	mpz_poly.o \
	ZmodF_poly.o


mpn_extras.o: mpn_extras.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpn_extras.c -o mpn_extras.o

Z.o: Z.c $(HEADERS)
	$(CC) $(CFLAGS) -c Z.c -o Z.o

memory-manager.o: memory-manager.c $(HEADERS)
	$(CC) $(CFLAGS) -c memory-manager.c -o memory-manager.o

Z_mpn.o: Z_mpn.c $(HEADERS)
	$(CC) $(CFLAGS) -c Z_mpn.c -o Z_mpn.o

ZmodF.o: ZmodF.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF.c -o ZmodF.o

ZmodF_mul.o: ZmodF_mul.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul.c -o ZmodF_mul.o

fmpz.o: fmpz.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz.c -o fmpz.o

fmpz_poly.o: fmpz_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly.c -o fmpz_poly.o

mpz_poly.o: mpz_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly.c -o mpz_poly.o

ZmodF_poly.o: ZmodF_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly.c -o ZmodF_poly.o

long_extras.o: long_extras.c long_extras.h
	$(CC) $(CFLAGS) -c long_extras.c -o long_extras.o



####### test program object files

test-support.o: test-support.c $(HEADERS)
	$(CC) $(CFLAGS) -c test-support.c -o test-support.o


fmpz_poly-test.o: fmpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly-test.c -o fmpz_poly-test.o

ZmodF-test.o: ZmodF-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF-test.c -o ZmodF-test.o

ZmodF_poly-test.o: ZmodF_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly-test.c -o ZmodF_poly-test.o

mpz_poly-test.o: mpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-test.c -o mpz_poly-test.o

Z_mpn-test.o: Z_mpn-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c Z_mpn-test.c -o Z_mpn-test.o

ZmodF_mul-test.o: ZmodF_mul-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-test.c -o ZmodF_mul-test.o

long_extras-test.o: long_extras-test.c 
	$(CC) $(CFLAGS) -c long_extras-test.c -o long_extras-test.o



####### test program targets

Z_mpn-test: Z_mpn-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) Z_mpn-test.o test-support.o -o Z_mpn-test $(FLINTOBJ) $(LIBS)

fmpz_poly-test: fmpz_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz_poly-test.o test-support.o -o fmpz_poly-test $(FLINTOBJ) $(LIBS)

ZmodF-test: ZmodF-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF-test.o test-support.o -o ZmodF-test $(FLINTOBJ) $(LIBS)

ZmodF_poly-test: ZmodF_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF_poly-test.o test-support.o -o ZmodF_poly-test $(FLINTOBJ) $(LIBS)

mpz_poly-test: mpz_poly-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) mpz_poly-test.o test-support.o -o mpz_poly-test $(FLINTOBJ) $(LIBS)

ZmodF_mul-test: ZmodF_mul-test.o test-support.o $(FLINTOBJ) $(HEADERS)
	$(CC) $(CFLAGS) ZmodF_mul-test.o test-support.o -o ZmodF_mul-test $(FLINTOBJ) $(LIBS)

long_extras-test: long_extras.o long_extras-test.o test-support.o memory-manager.o
	$(CC) $(CFLAGS) long_extras.o long_extras-test.o test-support.o memory-manager.o -o long_extras-test $(LIBS)



####### profiling object files

profiler.o: profiler.c $(HEADERS)
	$(CC) $(CFLAGS) -c profiler.c -o profiler.o

profiler-main.o: profiler-main.c $(HEADERS)
	$(CC) $(CFLAGS) -c profiler-main.c -o profiler-main.o


fmpz_poly-profile-tables.o: fmpz_poly-profile.c $(HEADERS)
	python make-profile-tables.py fmpz_poly
	$(CC) $(CFLAGS) -c fmpz_poly-profile-tables.c -o fmpz_poly-profile-tables.o
	rm fmpz_poly-profile-tables.c

fmpz_poly-profile.o: fmpz_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly-profile.c -o fmpz_poly-profile.o


mpz_poly-profile-tables.o: mpz_poly-profile.c $(HEADERS)
	python make-profile-tables.py mpz_poly
	$(CC) $(CFLAGS) -c mpz_poly-profile-tables.c -o mpz_poly-profile-tables.o
	rm mpz_poly-profile-tables.c

mpz_poly-profile.o: mpz_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly-profile.c -o mpz_poly-profile.o


ZmodF_poly-profile-tables.o: ZmodF_poly-profile.c $(HEADERS)
	python make-profile-tables.py ZmodF_poly
	$(CC) $(CFLAGS) -c ZmodF_poly-profile-tables.c -o ZmodF_poly-profile-tables.o
	rm ZmodF_poly-profile-tables.c

ZmodF_poly-profile.o: ZmodF_poly-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_poly-profile.c -o ZmodF_poly-profile.o


ZmodF_mul-profile-tables.o: ZmodF_mul-profile.c $(HEADERS)
	python make-profile-tables.py ZmodF_mul
	$(CC) $(CFLAGS) -c ZmodF_mul-profile-tables.c -o ZmodF_mul-profile-tables.o
	rm ZmodF_mul-profile-tables.c

ZmodF_mul-profile.o: ZmodF_mul-profile.c $(HEADERS)
	$(CC) $(CFLAGS) -c ZmodF_mul-profile.c -o ZmodF_mul-profile.o




####### profiling program targets

PROFOBJ = $(FLINTOBJ) profiler.o profiler-main.o

fmpz_poly-profile: fmpz_poly-profile.o fmpz_poly-profile-tables.o test-support.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o fmpz_poly-profile fmpz_poly-profile.o fmpz_poly-profile-tables.o test-support.o $(PROFOBJ) $(LIBS)

mpz_poly-profile: mpz_poly-profile.o mpz_poly-profile-tables.o test-support.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o mpz_poly-profile mpz_poly-profile.o mpz_poly-profile-tables.o test-support.o $(PROFOBJ) $(LIBS)

ZmodF_mul-profile: ZmodF_mul-profile.o ZmodF_mul-profile-tables.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o ZmodF_mul-profile ZmodF_mul-profile.o ZmodF_mul-profile-tables.o $(PROFOBJ) $(LIBS)

ZmodF_poly-profile: ZmodF_poly-profile.o ZmodF_poly-profile-tables.o $(PROFOBJ)
	$(CC) $(CFLAGS) -o ZmodF_poly-profile ZmodF_poly-profile.o ZmodF_poly-profile-tables.o $(PROFOBJ) $(LIBS)


kara-profile: kara-profile.c profiler.o test-support.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o kara-profile kara-profile.c profiler.o test-support.o $(FLINTOBJ) $(LIBS)



####### example programs

delta_qexp.o: delta_qexp.c $(HEADERS)
	$(CC) $(CFLAGS) -c delta_qexp.c -o delta_qexp.o

delta_qexp: delta_qexp.o $(FLINTOBJ)
	$(CC) $(CFLAGS) -o delta_qexp delta_qexp.o $(FLINTOBJ) $(LIBS)

BLTcubes: long_extras.o BLTcubes.c
	$(CC) $(CFLAGS) -o BLTcubes BLTcubes.c long_extras.o $(LIBS)

BPTJCubes: long_extras.o memory-manager.o
	$(CC) $(CFLAGS) -o BPTJCubes BPTJCubes.c memory-manager.o long_extras.o $(LIBS)

