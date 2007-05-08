CPP = gcc -std=c99
CC = gcc -std=c99
LINKOBJ  = ssmul.o driver.o
DRIVEROBJ = Z.o Zvec.o ssfft.o mpn_extras.o profiler.o ssmul.o Z-ssmul.o driver.o
DRIVER2OBJ = mpn_extras.o radixmul.o Z.o driver2.o
DRIVER3OBJ = profiler.o ssfft.o mpn_extras.o Z-ssmul.o driver3.o
TESTOBJ = ssfft.o mpn_extras.o Z.o ssmul.o Z-ssmul.o ssmul-test.o
PROFOBJ = mpn_extras.o ssmulprof.o test.o profiler.o
FLINTOBJ = Zvec.o Z.o ssmul.o ssfft.o Z-ssmul.o mpn_extras.o flint.o 
ifndef FLINT_GMP_INCLUDE_DIR
	FLINT_GMP_INCLUDE_DIR = "/home/dmharvey/gmp/install/include"
endif

ifndef FLINT_GMP_LIB_DIR
	FLINT_GMP_LIB_DIR = "/home/dmharvey/gmp/install/lib"
endif
LIBS = -L$(FLINT_GMP_LIB_DIR) -lgmp -lpthread
INCS =  -I"/usr/include" -I$(FLINT_GMP_INCLUDE_DIR)
CXXINCS = -I"/home/dmharvey/gmp/install/include" 
BIN  = ssmul
BIN2 = radixmul
BIN3  = Zssmul
FBIN = FLINT
CFLAGS = $(INCS) -funroll-loops -fexpensive-optimizations -O3
CXXFLAGS = $(CXXINCS) -funroll-loops -fexpensive-optimizations -O3
RM = rm -f

.PHONY: all all-before all-after clean clean-custom driver test

all: all-before ssmul all-after

driver: $(DRIVEROBJ)
	$(CPP) $(DRIVEROBJ) -o $(BIN) $(LIBS)

driver2: $(DRIVER2OBJ)
	$(CPP) $(DRIVER2OBJ) -o $(BIN2) $(LIBS)

driver3: $(DRIVER3OBJ)
	$(CPP) $(DRIVER3OBJ) -o $(BIN3) $(LIBS)

test: $(TESTOBJ)
	$(CPP) $(TESTOBJ) -o $(BIN) $(LIBS)

modpmul-test: modpmul-test.o modpmul.o mpn_extras.o profiler.o
	$(CPP) modpmul-test.o modpmul.o mpn_extras.o profiler.o -o modpmul-test $(LIBS)

profile: $(PROFOBJ)
	$(CPP) $(PROFOBJ) -o $(BIN) $(LIBS)

ssmul-profile: Zvec.o Z.o  ssfft.o ssmul-profile.o profiler.o Z-ssmul.o ssmul.o mpn_extras.o
	$(CPP) -o ssmul-profile ssfft.o Zvec.o Z.o Z-ssmul.o ssmul-profile.o profiler.o ssmul.o mpn_extras.o $(LIBS)

flint: $(FLINTOBJ)
	$(CPP) $(FLINTOBJ) -o $(FBIN) $(LIBS)

clean: clean-custom
	${RM} $(FLINTOBJ) $(PROFOBJ) $(TESTOBJ) $(DRIVEROBJ) $(BIN) $(FBIN)

$(BIN): $(LINKOBJ)
	$(CPP) $(LINKOBJ) -o "ssmul" $(LIBS)

ssmulprof.o: ssmul.c
	$(CPP) -c ssmul.c -o ssmulprof.o -DPROFILE $(CXXFLAGS)

ssmul.o: ssmul.c ssmul.h
	$(CPP) -c ssmul.c -o ssmul.o $(CXXFLAGS)

Z-ssmul.o: Z-ssmul.c Z-ssmul.h
	$(CPP) -c Z-ssmul.c -o Z-ssmul.o $(CXXFLAGS)

ssmul-test.o: ssmul-test.c Zvec.o
	$(CPP) -c ssmul-test.c -o ssmul-test.o $(CXXFLAGS)

ssmul-profile.o: ssmul-profile.cpp
	$(CPP) -c ssmul-profile.cpp -o ssmul-profile.o $(CXXFLAGS)

radixmul.o: radixmul.c radixmul.h
	$(CPP) -c radixmul.c -o radixmul.o $(CXXFLAGS)

modpmul.o: modpmul.c modpmul.h
	$(CPP) -c modpmul.c -o modpmul.o $(CXXFLAGS)

modpmul-test.o: modpmul-test.c modpmul.h profiler.h
	$(CPP) -c modpmul-test.c -o modpmul-test.o $(CXXFLAGS)

mpn_extras.o: mpn_extras.c mpn_extras.h
	$(CPP) -c mpn_extras.c -o mpn_extras.o $(CXXFLAGS)

driver.o: driver.c
	$(CPP) -c driver.c -o driver.o $(CXXFLAGS)

driver2.o: driver2.c
	$(CPP) -c driver2.c -o driver2.o $(CXXFLAGS)

driver3.o: driver3.c
	$(CPP) -c driver3.c -o driver3.o $(CXXFLAGS)

flint.o: test.c
	$(CPP) -c test.c -o flint.o $(CXXFLAGS)

Zvec.o: Zvec.c
	$(CPP) -c Zvec.c -o Zvec.o $(CXXFLAGS)

Z.o: Z.c Z.h
	$(CPP) -c Z.c -o Z.o $(CXXFLAGS)

profiler.o: profiler.c profiler.h
	$(CC) -c profiler.c -o profiler.o $(CFLAGS)
	
flint-manager.o: flint-manager.c flint-manager.h
	$(CC) -c flint-manager.c -o flint-manager.o $(CFLAGS)

Zpoly.o: Zpoly.c Zpoly.h flint-manager.o 
	$(CC) -c Zpoly.c -o Zpoly.o $(CFLAGS)
	
Zpoly-test.o: Zpoly-test.c Zpoly.h flint-manager.o 
	$(CC) -c Zpoly-test.c -o Zpoly-test.o $(CFLAGS)
	
Zpoly-test: Zpoly-test.o Zpoly.o flint-manager.o 
	$(CC) Zpoly.o Zpoly-test.o flint-manager.o -o Zpoly-test $(CFLAGS) $(LIBS)
	
Zpoly_mpn.o: Zpoly_mpn.c Zpoly_mpn.h Zpoly.h flint-manager.o 
	$(CC) -c Zpoly_mpn.c -o Zpoly_mpn.o $(CFLAGS)
	
Zpoly_mpn-test.o: Zpoly_mpn-test.c Zpoly_mpn.h Zpoly.h flint-manager.o 
	$(CC) -c Zpoly_mpn-test.c -o Zpoly_mpn-test.o $(CFLAGS)
	
Zpoly_mpn-test: Zpoly_mpn-test.o mpn_extras.o Zpoly.o Zpoly_mpn.o Zpoly.h flint-manager.o 
	$(CC) Zpoly_mpn.o mpn_extras.o Zpoly.o Zpoly_mpn-test.o flint-manager.o -o Zpoly_mpn-test $(CFLAGS) $(LIBS)


ZmodF.o: ZmodF.c ZmodF.h
	$(CC) -c ZmodF.c -o ZmodF.o $(CFLAGS)
	
ZmodF-test.o: ZmodF-test.c ZmodF.h
	$(CC) -c ZmodF-test.c -o ZmodF-test.o $(CFLAGS)

ZmodF-test: ZmodF-test.o ZmodF.o
	$(CC) ZmodF-test.o ZmodF.o -o ZmodF-test $(CFLAGS) $(LIBS)

ZmodFpoly.o: ZmodFpoly.c ZmodFpoly.h ZmodF.h
	$(CC) -c ZmodFpoly.c -o ZmodFpoly.o $(CFLAGS)
	
ZmodFpoly-test.o: ZmodFpoly-test.c ZmodFpoly.h ZmodF.h
	$(CC) -c ZmodFpoly-test.c -o ZmodFpoly-test.o $(CFLAGS)

ZmodFpoly-test: Zpoly_mpn.o flint-manager.o mpn_extras.o ZmodFpoly-test.o ZmodFpoly.o ZmodF.o
	$(CC) Zpoly_mpn.o flint-manager.o mpn_extras.o ZmodFpoly-test.o ZmodFpoly.o ZmodF.o -o ZmodFpoly-test $(CFLAGS) $(LIBS)


ssfft.o: ssfft.c ssfft.h
	$(CC) -c ssfft.c -o ssfft.o $(CFLAGS)

ssfft-test.o: ssfft-test.c ssfft.h
	$(CC) -c ssfft-test.c -o ssfft-test.o $(CFLAGS)

ssfft-test: ssfft.o ssfft-test.o
	$(CC) -m64 ssfft.o ssfft-test.o -o ssfft-test $(LIBS)

ssfft-profile.o: ssfft-profile.c
	$(CC) -c ssfft-profile.c -o ssfft-profile.o $(CFLAGS)

ssfft-profile: ssfft.o ssfft-profile.o profiler.o
	$(CC) ssfft.o ssfft-profile.o profiler.o -o ssfft-profile $(LIBS)
