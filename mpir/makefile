LIBDIR=$(PREFIX)/lib
INCLUDEDIR=$(PREFIX)/include
DOCDIR=$(PREFIX)/doc

ifndef MPIR_CC
	MPIR_CC = gcc
endif

ifeq ($(MAKECMDGOALS),library)
	CC = $(MPIR_CC) -fPIC -std=c99
else
	CC = $(MPIR_CC) -std=c99
endif

ifndef MPIR_PY
	MPIR_PY = python
endif

ifndef MPIR_CPP
	MPIR_CPP = g++
endif

CPP = $(MPIR_CPP) 

LIBS = -L$(MPIR_GMP_LIB_DIR) $(MPIR_LINK_OPTIONS) -lgmp -lpthread -lm

INCS = -I$(MPIR_GMP_INCLUDE_DIR) 

CFLAGS = $(INCS) $(MPIR_TUNE) -O3

RM = rm -f

HEADERS = \
	mpir.h \
	longlong.h \
	longlong_wrapper.h \
	memory_manager.h \
	fmpz.h \
	mpn_extras.h \
	test_support.h \
	fmpz_poly.h \
	mpz_poly.h \
        fmpz_mat.h
	

####### library object files

MPIROBJ = \
	memory_manager.o \
	fmpz.o \
	fmpz_poly.o \
	mpz_poly.o \
        fmpz_mat.o

test: memory_manager-test fmpz-test fmpz_poly-test 
	./memory_manager-test
	./fmpz-test
	./fmpz_poly-test
	./fmpz_mat-test

all: test

library: $(MPIR_LIB)

libmpir.dylib: $(MPIROBJ)
	$(CC) -single_module -fPIC -dynamiclib -o libmpir.dylib $(MPIROBJ) $(LIBS)

libmpir.dll: $(MPIROBJ)
	$(CC) -fPIC -shared -o libmpir.dll $(MPIROBJ) $(LIBS)

libmpir.so: $(MPIROBJ)
	$(CC) -fPIC -shared -o libmpir.so $(MPIROBJ) $(LIBS)

memory_manager.o: memory_manager.c $(HEADERS)
	$(CC) $(CFLAGS) -c memory_manager.c -o memory_manager.o

fmpz.o: fmpz.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz.c -o fmpz.o

fmpz_poly.o: fmpz_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly.c -o fmpz_poly.o

mpz_poly.o: mpz_poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c mpz_poly.c -o mpz_poly.o

fmpz_mat.o: fmpz_mat.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_mat.c -o fmpz_mat.o

####### test program object files

test_support.o: test_support.c test_support.h $(HEADERS)
	$(CC) $(CFLAGS) -c test_support.c -o test_support.o

memory_manager-test.o: memory_manager-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c memory_manager-test.c -o memory_manager-test.o

fmpz-test.o: fmpz-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz-test.c -o fmpz-test.o

fmpz_poly-test.o: fmpz_poly-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_poly-test.c -o fmpz_poly-test.o

fmpz_mat-test.o: fmpz_mat-test.c $(HEADERS)
	$(CC) $(CFLAGS) -c fmpz_mat-test.c -o fmpz_mat-test.o

####### test program targets

memory_manager-test: memory_manager-test.o $(MPIROBJ) $(HEADERS)
	$(CC) $(CFLAGS) memory_manager-test.o -o memory_manager-test $(MPIROBJ) $(LIBS)

fmpz-test: fmpz-test.o test_support.o $(MPIROBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz-test.o -o fmpz-test $(MPIROBJ) test_support.o $(LIBS) 

fmpz_poly-test: fmpz_poly-test.o test_support.o $(MPIROBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz_poly-test.o -o fmpz_poly-test $(MPIROBJ) test_support.o $(LIBS)

fmpz_mat-test: fmpz_mat-test.o test_support.o $(MPIROBJ) $(HEADERS)
	$(CC) $(CFLAGS) fmpz_mat-test.o -o fmpz_mat-test $(MPIROBJ) test_support.o $(LIBS)

####### tuning object files 

####### tuning program targets

####### profiling object files

####### profiling program targets

####### example programs
