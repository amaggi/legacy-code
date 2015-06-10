#!/usr/bin/env python
# encoding: utf-8

from glob import glob
from plot_dsp import read_dsp, plot_all_map

print "Plotting coverage using all dsp files for all events"

periods = [20, 50, 100, 150]

# get all the filenames
dsp_files = glob('*/*.dsp')
# create a list for all the files
dsp_list = []

# loop over the files
for fname in dsp_files:
    # read the dispersion
    dsp, dsp_dict = read_dsp(fname)
    # add the information to the list
    dsp_list.append((dsp, dsp_dict))
 
# plot the coverage map
plot_all_map(dsp_list, coverage=True)
for per in periods:
    print "Treating period %.1f"%per
    plot_all_map(dsp_list, coverage=True, period=per)
