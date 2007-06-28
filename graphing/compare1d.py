# script for comparing FLINT 1d profile output
# requires matplotlib to be installed
#
# (C) 2007 David Harvey, GPL license yadda yadda

import sys
from optparse import OptionParser
import pylab
from math import log


parser = OptionParser(usage="python compare1d.py [options] input1 [input2 ...]")

parser.add_option("-o", "--output",
                  default="graph.png",
                  help="output filename (default graph.png)")
                  
parser.add_option("-c", "--config",
                  default="compare1d_config",
                  help="configuration file (default compare1d_config)")


(options, args) = parser.parse_args()

if len(args) == 0:
   parser.print_help()
   sys.exit()


# default configuration settings;
# see compare1d_config.py for explanation of each setting
CONFIG_dpi = 96
CONFIG_dotsize = 3
CONFIG_title = None
CONFIG_tolerance = 1.05
CONFIG_xscale = None
CONFIG_yscale = None
CONFIG_xlabel = None
CONFIG_ylabel = None
CONFIG_truncate = 2.0


# override settings from configuration file if requested
try:
   config_module = __import__(options.config)

   # import all variables starting with CONFIG_ into global namespace
   for var in dir(config_module):
      if var.startswith("CONFIG_"):
         globals()[var] = config_module.__dict__[var]

except ImportError:
   # configuration file not found
   print "Warning: could not find configuration file \"%s.py\"; using default settings" % (options.config)



if CONFIG_title is not None:
   pylab.title(CONFIG_title)


log_strings = {None: "", "log10" : " (log10)", "log2" : " (log2)"}
if CONFIG_xlabel is None:
   CONFIG_xlabel = ""
if CONFIG_ylabel is None:
   CONFIG_ylabel = ""
CONFIG_xlabel += log_strings[CONFIG_xscale]
CONFIG_ylabel += log_strings[CONFIG_yscale]

pylab.xlabel(CONFIG_xlabel)
pylab.ylabel(CONFIG_ylabel)


files = [open(arg) for arg in args]
data = []

for f in files:
   # skip up to "==========" where the data starts
   iter = f.__iter__()
   while not iter.next().startswith("====="):
      pass
   # decode each quadruple (x, min, max)
   points = {}
   for line in iter:
      fields = line.split()
      x = float(fields[0])
      y_min = float(fields[1])
      y_max = float(fields[2])

      # ignore points where the max and min are too far apart
      if y_max/y_min <= CONFIG_tolerance:
         points[x] = y_min

   data.append(points)



for i in range(len(data)):
   points = data[i]
   L = list(points.iteritems())
   L.sort()
   xvals = [x for (x, y) in L]
   yvals = [y for (x, y) in L]
   pylab.plot(xvals, yvals, label=args[i])


pylab.legend()



pylab.savefig(options.output, dpi=CONFIG_dpi)


# here's some code I was playing around with to generate a legend, just
# leaving it here for now because I don't want to lose it
#
# import matplotlib
# a = matplotlib.patches.Ellipse((0.9, 0.9), 0.2, 0.2, 360)
# a.set_facecolor((0.5, 0.5, 0.5))
# a.set_edgecolor((0.5, 0.5, 0.5))
# a.set_transform(pylab.gcf().transFigure)
# 
# b = matplotlib.patches.Ellipse((4.0, 0.5), 1.0, 1.0, 360)
# b.set_facecolor((0.5, 0.5, 0.5))
# b.set_edgecolor((0.5, 0.5, 0.5))
# 
# pylab.gca().add_patch( a )
# pylab.gca().add_patch( b )
# 
#for i in range(10):
#   print cmap(i/10.0)
