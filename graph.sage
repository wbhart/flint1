# graph generating script for FLINT development and profiling
# (C) 2006 David Harvey, GPL yadda yadda yadda

from sage.plot.plot import GraphicPrimitive

# use python's faster log function
from math import log

# reads input file into a dictionary whose keys are (deg, coeff)
def readfile(src, normalise=None):
   data = {}
   for line in src:
      (deg, coeff, val, error) = line.split()
      deg = int(deg)
      coeff = int(coeff)
      val = float(val)
      try:
         if normalise == "n":
            val = val - log(deg * coeff)/log(10)
         if normalise == "nlogx":
            val = val - (log(deg * coeff) + log(log(deg)))/log(10)
         elif normalise == "nlogy":
            val = val - (log(deg * coeff) + log(log(coeff)))/log(10)
         elif normalise == "nlogn":
            val = val - (log(deg * coeff) + log(log(deg * coeff)))/log(10)
         data[deg, coeff] = val
      except OverflowError:
         # ignore log(0) errors
         pass
   return data


#####################################################################
### optimised implementation of GraphicPrimitive_Point

class GraphicPrimitive_FastPoint(GraphicPrimitive):
   def __init__(self, x, y, pointsize, rgbcolor):
      self.x = float(x)
      self.y = float(y)
      self.pointsize = pointsize
      self.rgbcolor = rgbcolor

   def _render_on_subplot(self, subplot):
      subplot.scatter([self.x], [self.y], s=self.pointsize,
                      c=self.rgbcolor, alpha=1, faceted=False)


fast_point = GraphicPrimitive_FastPoint


#####################################################################


from optparse import OptionParser

parser = OptionParser(usage="sage graph.sage [options] input1 [input2]\n\nOne input implies single-plot mode.\nTwo inputs imply comparison plot mode.\nRun sage graph.sage -h for help.")

parser.add_option("-o", "--output",
                  default="graph.png",
                  help="output filename (default graph.png)")

parser.add_option("-v", "--preview",
                  type="int",
                  default=1,
                  help="preview mode, e.g. '-v 5' will only display 1/5th of the data points")

parser.add_option("-p", "--pixelsize",
                  type="int",
                  default=20,
                  help="pixel size (default 20)")

parser.add_option("-n", "--normalise",
                  help="normalise values by one of the functions: 'n' = x*y, 'nlogx' = x*y*log(x), 'nlogy' = x*y*log(y), 'nlogn' = x*y*log(x*y) [only used in single-plot mode]")

parser.add_option("-t", "--tolerance",
                  type="float",
                  default=1.03,
                  help="ratios less than this are considered equal (default 1.03) [only used in comparison plot mode]")

parser.add_option("-c", "--crop",
                  type="float",
                  default=2.0,
                  help="ratios larger than this are cropped (default 2.0) [only used in comparison plot mode]")

parser.add_option("-z", "--zero",
                  action="store_true",
                  dest="zero",
                  help="compare a single input to the zero data set")

(options, args) = parser.parse_args()

if len(args) == 0:
   parser.error("no input files specified")
if len(args) > 2:
   parser.error("too many input files specified")


pixelsize = int(options.pixelsize)
normalise = options.normalise
if normalise not in [None, "n", "nlogx", "nlogy", "nlogn"]:
   parser.error("normalisation must be one of n, nlogx, nlogy, nlogn")


if len(args) == 1 and not options.zero:
   # single graph case

   print
   print "reading data..."
         
   data = readfile(file(args[0]), normalise)

   # figure out the top and bottom 2% quantiles
   vals = data.values()
   vals.sort()
   cutoff_lo = vals[int(len(vals) * 0.02)]
   cutoff_hi = vals[int(len(vals) * 0.98)]

   G = Graphics()

   print "building graphics object..."
   xvals = []
   yvals = []
   count = 0
   for ((deg, coeff), val) in data.iteritems():
      count = count + 1
      if count != options.preview:
         continue
      count = 0
      
      # scale so that cutoff_hi is 90% white and cutoff_lo is 20% white
      if val < cutoff_lo:
         val = cutoff_lo
      elif val > cutoff_hi:
         val = cutoff_hi
      value = (val - cutoff_lo) / (cutoff_hi - cutoff_lo) * 0.7 + 0.2
      colour = (value, value, value)

      x = log(deg)/log(2)
      y = log(coeff)/log(2)
      xvals.append(x)
      yvals.append(y)

      p = fast_point(x, y, pixelsize, colour)
      G.append(p)

   G._Graphics__xmin = min(xvals)
   G._Graphics__xmax = max(xvals)
   G._Graphics__ymin = min(yvals)
   G._Graphics__ymax = max(yvals)

   
else:
   # comparison case (two inputs)

   print
   print "reading data..."
         
   data1 = readfile(file(args[0]))
   if options.zero:
      # make up dummy data set of all zeroes
      data2 = {}
      for (deg, coeff) in data1.keys():
         data2[deg, coeff] = 0.0
   else:
      data2 = readfile(file(args[1]))

   # collect data that's in BOTH datasets
   data = {}
   for (deg, coeff) in data1.keys():
      if (deg, coeff) in data2:
         data[deg, coeff] = (data1[deg, coeff], data2[deg, coeff])

   G = Graphics()

   print "building graphics object..."
   xvals = []
   yvals = []
   count = 0
   for (deg, coeff), (val1, val2) in data.iteritems():
      count = count + 1
      if count != options.preview:
         continue
      count = 0
      
      # scale so that:
      # crop value gets intensity 1.0
      # tolerance value gets intensity 0.3
      # anything less than tolerance is a boring grey colour
      value = abs(val2 - val1) * log(10)
      if value < log(options.tolerance):
         colour = (0.3, 0.3, 0.3)
      else:
         if value > log(options.crop):
            value = 1.0
         else:
            value = value / log(options.crop)
         value = value * 0.7 + 0.3
         if val1 < val2:
            colour = (value/3, value/3, value)
         else:
            colour = (value, value/3, value/3)

      x = log(deg)/log(2)
      y = log(coeff)/log(2)
      xvals.append(x)
      yvals.append(y)

      p = fast_point(x, y, pixelsize, colour)
      G.append(p)

   G._Graphics__xmin = min(xvals)
   G._Graphics__xmax = max(xvals)
   G._Graphics__ymin = min(yvals)
   G._Graphics__ymax = max(yvals)


print "rendering to %s...." % options.output
G.save(options.output, dpi=200, figsize=[5, 5])
if len(args) == 1:
   print "finished."
else:
   print "finished (blue = %s, red = %s)." % (args[0], args[1])
   

########### end of file
