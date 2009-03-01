#
# This python script is called by configure to generate the makefile
# for zn_poly.
#

# --------------------------------------------------------------------------
# various lists of modules


# These are the modules that go into the actual zn_poly library. They get
# compiled in both optimised and debug modes.
lib_modules = ["zn_mod", "misc", "mul_ks", "pack", "mul", "midmul", "tuning",
               "pmf", "mul_fft", "mul_fft_dft", "midmul_fft", "nussbaumer",
               "array", "invert"]

# These are modules containing various test routines. They get compiled
# in debug mode only.
test_modules = ["test", "mul_ks-test", "pack-test", "mul_fft-test",
                "midmul_fft-test", "nussbaumer-test", "ref_mul", "invert-test"]

# These are modules containing various profiling routines. They get compiled
# in optimised mode only.
prof_modules = ["prof_main", "profiler", "mul-profile", "negamul-profile", 
                "array-profile", "midmul-profile", "invert-profile"]

# These are modules containing profiling routines for NTL. They get compiled
# in optimised mode only, with the C++ compiler.
cpp_prof_modules = ["ntl-profile"]

# These are modules containing dummy routines replacing the NTL ones, if
# we're not compiling with NTL support.
noncpp_prof_modules = ["ntl-profile-dummy"]

# These are modules shared by the test and profiling code. They get compiled
# in both debug and optimised mode.
testprof_modules = ["support"]

# Profiling targets. Each X has a main file "X-main.c" which is linked
# against prof_main.c. They are compiled once with PROFILE_NTL defined
# and once without.
prof_progs = ["mul-profile", "midmul-profile", "negamul-profile",
              "array-profile", "invert-profile"]

# These are modules used in the tuning program; they get compiled only in
# optimised mode.
tune_modules = ["tune", "mul_ks-tune", "mul-tune", "midmul-tune",
                "nussbaumer-tune"]

# Example programs.
example_progs = ["bernoulli"]

# These are the headers that need to be copied to the install location.
install_headers = ["zn_poly.h", "wide_arith.h"]

# These are the other headers.
other_headers = ["support.h", "profiler.h", "zn_poly_internal.h"]


# --------------------------------------------------------------------------
# read command line options

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--prefix", dest="prefix", default="/usr/local")
parser.add_option("--cflags", dest="cflags", default="-O2")
parser.add_option("--ldflags", dest="ldflags", default="")
parser.add_option("--gmp-prefix", dest="gmp_prefix", default="/usr/local")
parser.add_option("--ntl-prefix", dest="ntl_prefix", default="/usr/local")
parser.add_option("--use-flint", dest="use_flint", action="store_true", default=False)
parser.add_option("--flint-prefix", dest="flint_prefix", default="/usr/local")

options, args = parser.parse_args()


gmp_include_dir = options.gmp_prefix + "/include"
gmp_lib_dir = options.gmp_prefix + "/lib"

ntl_include_dir = options.ntl_prefix + "/include"
ntl_lib_dir = options.ntl_prefix + "/lib"

if options.use_flint:
   flint_include_dir = options.flint_prefix + "/include"
   flint_lib_dir = options.flint_prefix + "/lib"

cflags = options.cflags
ldflags = options.ldflags
prefix = options.prefix

includes = "-I" + gmp_include_dir
libs = "-L" + gmp_lib_dir + " -lgmp -lm"

if options.use_flint:
   includes = includes + " -I" + flint_include_dir
   libs = libs + " -L" + flint_lib_dir + " -lflint"
   cflags = cflags + " -std=c99 -DZNP_USE_FLINT"

cpp_includes = includes + " -I" + ntl_include_dir
cpp_libs = libs + " -L" + ntl_lib_dir + " -lntl"


# --------------------------------------------------------------------------
# generate the makefile

import time

print "#"
print "# Do not edit directly -- this file was auto-generated"
print "# by makemakefile.py on " + time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())
print "#"
print

print
print "CC = gcc"
print "CFLAGS = " + cflags
print "LDFLAGS = " + ldflags
print "INCLUDES = " + includes
print "LIBS = " + libs

print
print "CPP = g++"
print "CPPFLAGS = " + cflags
print "CPP_INCLUDES = " + cpp_includes
print "CPP_LIBS = " + cpp_libs

print
print "HEADERS = " + " ".join(install_headers + other_headers)
print "LIBOBJS = " + " ".join([x + ".o" for x in lib_modules])
print "TESTOBJS = " + " ".join([x + "-DEBUG.o" for x in lib_modules + test_modules + testprof_modules])
print "PROFOBJS = " + " ".join([x + ".o" for x in lib_modules + prof_modules + noncpp_prof_modules + testprof_modules])
print "CPP_PROFOBJS = " + " ".join([x + ".o" for x in lib_modules + prof_modules + cpp_prof_modules + testprof_modules])
print "TUNEOBJS = " + " ".join([x + ".o" for x in lib_modules + tune_modules + testprof_modules + prof_modules + noncpp_prof_modules if x != "prof_main"])

print
print
print "all: libzn_poly.a"
print
print "install:"
print "\tmkdir -p %s/include/zn_poly" % prefix
print "\tmkdir -p %s/lib" % prefix
print "\tcp libzn_poly.a %s/lib" % prefix
print "\tcp zn_poly.h %s/include/zn_poly" % prefix
print "\tcp wide_arith.h %s/include/zn_poly" % prefix
print
print "clean:"
print "\trm -f *.o"
print "\trm -f libzn_poly.a"
print "\trm -f libzn_poly.dylib"
print "\trm -f libzn_poly.so"
print "\trm -f test"
print "\trm -f tune"
for x in prof_progs:
   print "\trm -f " + x
   print "\trm -f " + x + "-ntl"
for x in example_progs:
   print "\trm -f " + x


print
print
print "##### library targets"
print
print "libzn_poly.a: $(LIBOBJS)"
print "\tar -r libzn_poly.a $(LIBOBJS)"
print "\tranlib libzn_poly.a"
print
print "libzn_poly.dylib: $(LIBOBJS)"
print "\t$(CC) -single_module -fPIC -dynamiclib -o libzn_poly.dylib $(LIBOBJS) $(LIBS)"
print
print "libzn_poly.so: $(LIBOBJS)"
print "\t$(CC) -shared -o libzn_poly.so $(LIBOBJS) $(LIBS)"

print
print
print "##### test program"
print
print "test: $(TESTOBJS) $(HEADERS)"
print "\t$(CC) -g $(LDFLAGS) -o test $(TESTOBJS) $(LIBS)"

print
print
print "##### profiling programs"
print
for x in prof_progs:
   print "%s-main.o: %s-main.c $(HEADERS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(INCLUDES) -DNDEBUG -o %s-main.o -c %s-main.c" % (x, x)
   print
   print "%s: %s-main.o $(PROFOBJS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(LDFLAGS) -o %s %s-main.o $(PROFOBJS) $(LIBS)" % (x, x)
   print
   print "%s-main-ntl.o: %s-main.c $(HEADERS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(INCLUDES) -DPROFILE_NTL -DNDEBUG -o %s-main-ntl.o -c %s-main.c" % (x, x)
   print
   print "%s-ntl: %s-main-ntl.o $(CPP_PROFOBJS)" % (x, x)
   print "\t$(CPP) $(CPPFLAGS) $(LDFLAGS) -o %s-ntl %s-main-ntl.o $(CPP_PROFOBJS) $(CPP_LIBS)" % (x, x)
   print

print
print
print "##### tuning utility"
print
print "tune: $(TUNEOBJS)"
print "\t$(CC) $(CFLAGS) $(LDFLAGS) -o tune $(TUNEOBJS) $(LIBS)"


print
print
print "##### example programs"
for x in example_progs:
   print
   print "%s: %s.o $(LIBOBJS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(LDFLAGS) -o %s %s.o $(LIBOBJS) $(LIBS)" % (x, x)


print
print
print "##### object files (with debug code)"
for x in lib_modules + test_modules + testprof_modules + example_progs:
   print
   print "%s-DEBUG.o: %s.c $(HEADERS)" % (x, x)
   print "\t$(CC) -g $(CFLAGS) $(INCLUDES) -DDEBUG -o %s-DEBUG.o -c %s.c" % (x, x)

print
print
print "##### object files (no debug code)"
for x in lib_modules + prof_modules + testprof_modules + tune_modules + example_progs:
   print
   print "%s.o: %s.c $(HEADERS)" % (x, x)
   print "\t$(CC) $(CFLAGS) $(INCLUDES) -DNDEBUG -o %s.o -c %s.c" % (x, x)

print
print
print "##### object files (C++, no debug code)"
for x in cpp_prof_modules:
   print
   print "%s.o: %s.c $(HEADERS)" % (x, x)
   print "\t$(CPP) $(CPPFLAGS) $(CPP_INCLUDES) -DNDEBUG -o %s.o -c %s.c" % (x, x)


### end of file
