#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from mpl_toolkits.basemap import Basemap
from dsp2intomodes import period2key

periods = [20, 50, 100, 150]

# set parameters for the projection
minlat = -40.0
maxlat = 60.0
minlon = -60.0
maxlon = 120.0
clat = np.average([minlat, maxlat])
clon = np.average([minlon, maxlon])
parallels = np.arange(minlat-20., maxlat+20., 20.)
meridians = np.arange(minlon-20., maxlon+20., 20.)

# for local grid : spacing between points in meters 
spacing = 5000.

m = Basemap(projection='lcc', lon_0=clon, lat_0=clat, llcrnrlon=minlon,
            llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat,
            resolution='i', area_thresh=1000.)
nx = int((m.xmax-m.xmin)/spacing)+1
ny = int((m.ymax-m.ymin)/spacing)+1

# for global grid
lats = np.linspace(-89.5, 89.5, 180)
lons = np.linspace(-179.5, 179.5, 360)

def plot_gvel(per, relative=False):

    key = period2key(per)
    gvel_filename = "iso_gvel_fin_%s.dat"%key
    figname = "iso_gvel_fin_%s.png"%key

    # read the data (this method only works for 1x1 degree grid as loadtxt only
    # works whern there are the same number of points in each row)
    print gvel_filename
    rawdata = np.loadtxt(gvel_filename, skiprows=1)
    data = rawdata.ravel().reshape(180, 360)

    # project the data from the global grid to the projected one
    data_proj = m.transform_scalar(data[::-1, :], lons, lats, nx, ny)
    if relative:
        average = np.average(data_proj)
        data_proj = (data_proj/average - 1.0) * 100.0


    # do the plot
    plt.figure()

    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(parallels, labels=[1, 0, 0, 0])
    m.drawmeridians(meridians, labels=[0, 0, 0, 1])

    im = m.imshow(data_proj, cmap='jet_r')

    if relative:
        m.colorbar(location='right', label='%')
        plt.title('Group velocity anomaly at T=%.1f s (wrt U = %.2f km/s)' %
                 (per, average))
    else:
        m.colorbar(location='right', label='U (km/s)')
        plt.title('Group velocity anomaly at T=%.1f s'% per)

    plt.savefig(figname)
    plt.close()
    
    print figname


if __name__ == '__main__':

    for per in periods:
        plot_gvel(per, relative=False)


