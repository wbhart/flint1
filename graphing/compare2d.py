#######################################################################
#   This file is part of FLINT.
#
#   FLINT is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   FLINT is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with FLINT; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

# script for comparing FLINT 2d profile output
# requires matplotlib to be installed
#
# (C) 2007 David Harvey, GPL license yadda yadda

import sys
import matplotlib
matplotlib.use('Agg')
from optparse import OptionParser
from matplotlib.colors import LinearSegmentedColormap
import pylab
from math import log


parser = OptionParser(usage="python compare2d.py [options] input1 input2")

parser.add_option("-o", "--output",
                  default="graph.png",
                  help="output filename (default graph.png)")
                  
parser.add_option("-c", "--config",
                  default="compare2d_config",
                  help="configuration file (default compare2d_config)")


(options, args) = parser.parse_args()

if len(args) != 2:
   parser.print_help()
   sys.exit()


# default configuration settings;
# see compare2d_config.py for explanation of each setting
CONFIG_dpi = 96
CONFIG_dotsize = 3
CONFIG_title = None
CONFIG_tolerance = 1.05
CONFIG_xscale = None
CONFIG_yscale = None
CONFIG_xlabel = None
CONFIG_ylabel = None
CONFIG_truncate = 2.0
CONFIG_min_intensity = 0.5


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


# define colour map which makes:
#   values in [0, 0.5] are blue
#   values in [0.5, 1] are red
cmap_data = {
    'red'  :  ((0.0, 0.0, 0.0),
               (0.5, 0.0, CONFIG_min_intensity),
               (1.0, 1.0, 1.0)),
    'green':  ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
    'blue' :  ((0.0, 1.0, 1.0),
               (0.5, CONFIG_min_intensity, 0.0),
               (1.0, 0.0, 0.0))
}
cmap = LinearSegmentedColormap('????', cmap_data)


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


files = [open(args[0]), open(args[1])]
data = []

for f in files:
   # skip up to "==========" where the data starts
   iter = f.__iter__()
   while not iter.next().startswith("====="):
      pass
   # decode each quadruple (x, y, min, max)
   points = {}
   for line in iter:
      fields = line.split()
      x = float(fields[0])
      y = float(fields[1])
      z_min = float(fields[2])
      z_max = float(fields[3])

      # ignore points where the max and min are too far apart
      if z_max/z_min <= CONFIG_tolerance:
         points[(x, y)] = z_min

   data.append(points)


# merge data into single list of ratios
# values in [-1.0, 0.0] means the 1st sample was faster
# values in [0.0, 1.0] means the 2nd sample was faster
ratios = []
for (x, y) in data[0]:
   z0 = data[0][x, y]
   if (x, y) in data[1]:
      z1 = data[1][x, y]

      if CONFIG_xscale == "log10":
         x = log(x)/log(10.0)
      elif CONFIG_xscale == "log2":
         x = log(x)/log(2.0)

      if CONFIG_yscale == "log10":
         y = log(y)/log(10.0)
      elif CONFIG_yscale == "log2":
         y = log(y)/log(2.0)

      ratio = log(z1 / z0) / log(CONFIG_truncate)
      if ratio > 1.0:
         ratio = 1.0
      elif ratio < -1.0:
         ratio = -1.0

      ratios.append((x, y, ratio))

ratios.sort()

pylab.scatter([x for (x, y, z) in ratios], [y for (x, y, z) in ratios],
              c = [z for (x, y, z) in ratios], vmin=-1.0, vmax=1.0,
              s = CONFIG_dotsize, faceted=False, cmap = cmap)



# write down the two filenames (this should really be the profile name)
pylab.figtext(0.05, 0.03, args[0], fontsize=10, color=(1.0, 0, 0),
              verticalalignment="bottom")
pylab.figtext(0.05, 0.03, args[1], fontsize=10, color=(0, 0, 1.0),
              verticalalignment="top")


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
