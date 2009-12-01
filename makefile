CC=tcc
DEF=tiny_impdef
CFLAGS=
LIBS=-L../MPIR-tcc
INCS=-I${PWD} -I../MPIR-tcc -Izn_poly/src
export

SOURCES = zn_poly/src/zn_mod.c zn_poly/src/misc.o zn_poly/src/mul_ks.o \
          zn_poly/src/pack.o zn_poly/src/mul.o zn_poly/src/mulmid.o \
          zn_poly/src/mulmid_ks.o zn_poly/src/ks_support.o \
          zn_poly/src/mpn_mulmid.o zn_poly/src/nuss.o zn_poly/src/pmf.o \
          zn_poly/src/pmfvec_fft.o zn_poly/src/tuning.o zn_poly/src/mul_fft.o \
          zn_poly/src/mul_fft_dft.o zn_poly/src/array.o zn_poly/src/invert.o \
          fmpz.c fmpz_poly.c long_extras.c memory-manager.c mpn_extras.c \
          mpz_extras.c mpz_mat.c mpz_poly.c mpz_poly-tuning.c zmod_mat.c \
          zmod_poly.c ZmodF.c ZmodF_mul.c ZmodF_mul-tuning.c ZmodF_poly.c

TESTS = mpn_extras-test.exe long_extras-test.exe ZmodF-test.exe \
	ZmodF_poly-test.exe ZmodF_mul-test.exe fmpz-test.exe mpz_poly-test.exe \
	zmod_mat-test.exe zmod_poly-test.exe fmpz_poly-test.exe

OBJS = $(patsubst %.c, %.o, $(SOURCES))

.c.o:
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

%.exe: %.c test-support.o profiler.o
	$(CC) $(CFLAGS) $(INCS) profiler.o test-support.o $< flint.def mpir.def -o $@

check: all $(TESTS)
	$(foreach name, $(TESTS), $(name);)

all: $(OBJS)
	$(CC) -shared -rdynamic *.o zn_poly/src/*.o $(LIBS) mpir.def -o flint.dll
	$(DEF) flint.dll

clean:
	rm -f *.o flint.dll flint.def zn_poly/src/*.o *.exe

.PHONY: check all clean
