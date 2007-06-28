#******************************************************************************
#
#  default configuration file for compare.py
#
#  Please don't edit this file; make a copy and then tell compare.py about it
#  by running e.g.
#
#      python compare.py -c my_config.py
#
#******************************************************************************


#------------------------------------------------------------------------------
# dots per inch for output file
# larger values make a bigger graphics file

CONFIG_dpi = 96


#------------------------------------------------------------------------------
# title at the top of the image

CONFIG_title = None


#------------------------------------------------------------------------------
# ignore all data points whose minimum and maximum times differ by a factor
# of more than CONFIG_tolerance

CONFIG_tolerance = 1.05


#------------------------------------------------------------------------------
# xscale and yscale determine the scaling used on each axis.
# Possible values are None, "log10", "log2"

CONFIG_xscale = None
CONFIG_yscale = None


#------------------------------------------------------------------------------
# text labels for each axis (None to leave blank)

CONFIG_xlabel = None
CONFIG_ylabel = None


################# end of config file
