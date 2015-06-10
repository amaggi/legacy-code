#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np
from glob import glob
from sklearn.cluster import Birch
from plot_dsp import read_dsp, plot_all_dsp, plot_all_map

# make an output directory for clusters
dirname = 'clusters'
if not os.path.isdir(dirname):
	os.mkdir(dirname)

# distance in degrees for clusters
dist = 2.0

# get all filenames
dsp_files = glob('*/*.dsp')
n_dsp = len(dsp_files)

# set up the data structures
event_array = np.empty((n_dsp,2), dtype=float)
station_list=[]
dsp_list = []

# read all the dispersion files
for i in xrange(n_dsp):
	fname = dsp_files[i]
	dsp, dsp_dict = read_dsp(fname)
	dsp_list.append((dsp, dsp_dict))

	station_list.append(dsp_dict["STA"])
	event_array[i, 0] = dsp_dict["EVLA"]
	event_array[i, 1] = dsp_dict["EVLO"]
station_array = np.array(station_list)
dsp_array = np.array(dsp_list)

# extract the unique station names
stations = np.unique(station_array)
print stations

for sta in stations:
	events = event_array[station_array==sta, :]
	dsp_shortlist = dsp_array[station_array==sta]
	print sta, events.shape, dsp_shortlist.shape

	# cluster on events so as to compare dispersion curves for nearby
	# events
	brc = Birch(branching_factor=50, n_clusters=None, threshold=dist,
		    compute_labels=True)
	brc.fit(events)
	labels = brc.predict(events)
	print np.max(labels)
	for lab in np.unique(labels):
		dsp_this_label_list = dsp_shortlist[labels==lab]
		cluster_name = os.path.join(dirname, 'cluster_%s_%03d'%(sta, lab))
		plot_all_dsp(dsp_this_label_list, legend=False,
			     fname='%s_gvel.png'%cluster_name)
		plot_all_map(dsp_this_label_list, fname='%s_map.png'%cluster_name, legend=False)
		f = open('%s_info.txt'%cluster_name, 'w')
		for (dsp, dsp_dict) in dsp_this_label_list:
			f.write('%s %s %d %03d %02d %02d %.3f %.3f\n' %
			       (dsp_dict["STA"], dsp_dict["COMP"],
			        dsp_dict["YEAR"], dsp_dict["JDAY"],
				dsp_dict["HOUR"], dsp_dict["MIN"],
				dsp_dict["EVLA"], dsp_dict["EVLO"]))
		f.close()


