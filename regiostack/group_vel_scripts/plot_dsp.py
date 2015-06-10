#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d

# set some global plotting parameters
markersize = 15

minlat = -40.0
maxlat = 60.0
minlon = -60.0
maxlon = 120.0
clat = np.average([minlat, maxlat])
clon = np.average([minlon, maxlon])
parallels = np.arange(minlat-20., maxlat+20., 10.)
meridians = np.arange(minlon-20., maxlon+20., 10.)

m = Basemap(projection='lcc', lon_0=clon, lat_0=clat, llcrnrlon=minlon,
            llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat,
            resolution='i', area_thresh=1000.)

def read_dsp(fname):
    """
    Reads a dsp file and extracts the relevant information.
    """

    # set up list of names for columns in .dsp files
    names = list(["PROG", "RL", "UC", "MODE", "PER", "VEL", "ERR", "DIST",
                  "AZ", "AMP", "EVLA", "EVLO", "STLA", "STLO", "D1", "D2",
                  "D3", "D4", "STA", "COMP", "YEAR", "JDAY", "HOUR", "MIN"])
                  
    # read the file
    pd_dsp = pd.read_table(fname, sep="\s+", header=None, names=names)

    # extract the dispersion information
    dsp_names = list(["PER", "VEL", "ERR"])
    dsp = pd_dsp[dsp_names].values

    # extract the information about the station
    sta_names = (["STA", "COMP", "STLA", "STLO"])
    sta_list = pd_dsp[sta_names].values
    sta = sta_list[0, :]

    # extract the information about the event
    ev_names = (["YEAR", "JDAY", "HOUR", "MIN", "EVLA", "EVLO"])
    ev_list = pd_dsp[ev_names].values
    ev = ev_list[0, :]

    dsp_dict = {}
    for i in xrange(len(sta_names)):
        dsp_dict[sta_names[i]] = sta[i]
    for i in xrange(len(ev_names)):
        dsp_dict[ev_names[i]] = ev[i]

    return dsp, dsp_dict


def interpolate_dsp(dsp):
    # interpolates a dispersion curve
    # returns velocity and error at that period
    # raises ValueError if requested period is outside the band

    # interpret the dsp array
    per = dsp[:, 0]
    vel = dsp[:, 1] 
    err = dsp[:, 2] 

    # do the interpolation
    vfunc = interp1d(per, vel, kind='cubic', bounds_error=True,
                     assume_sorted=False)
    efunc = interp1d(per, err, kind='cubic', bounds_error=True,
                     assume_sorted=False)

    return vfunc, efunc

def plot_dsp(dsp, dsp_dict):
    plt.figure()
    plt.plot(dsp[:,0], dsp[:,1], 'bo')
    plt.errorbar(dsp[:,0], dsp[:,1], yerr=dsp[:,2], fmt="none")
    plt.xlim(0, 270)
    plt.ylim(2, 5.5)
    plt.xlabel('Period (s)')
    plt.ylabel('Group velocity (km/s)')
    plt.title("%d(%03d) %d:%d - %s.%s"%(dsp_dict["YEAR"], dsp_dict["JDAY"],
                                      dsp_dict["HOUR"], dsp_dict["MIN"],
                                      dsp_dict["STA"], dsp_dict["COMP"]))

    fname = "%d_%03d_%02d_%02d_%s_%s.png"%(dsp_dict["YEAR"], dsp_dict["JDAY"],
                                      dsp_dict["HOUR"], dsp_dict["MIN"],
                                      dsp_dict["STA"], dsp_dict["COMP"])
    print "Writing plot to file : %s"%fname
    plt.savefig(fname)
    plt.close()

def plot_all_dsp(dsp_list, legend=True, fname=None):
    plt.figure()
    for (dsp, dsp_dict) in dsp_list:
        label = "%s-%s"%(dsp_dict["STA"], dsp_dict["COMP"])
        plt.errorbar(dsp[:,0], dsp[:,1], yerr=dsp[:,2], fmt='o', label=label)
    plt.xlim(0, 270)
    plt.ylim(2, 5.5)
    plt.xlabel('Period (s)')
    plt.ylabel('Group velocity (km/s)')
    if legend:
        plt.legend()
    plt.title("%d(%03d) %02d:%02d"%(dsp_dict["YEAR"], dsp_dict["JDAY"],
                                      dsp_dict["HOUR"], dsp_dict["MIN"]))
    if fname is None:
        fname = "%d_%03d_%02d_%02d_all.png"%(dsp_dict["YEAR"], dsp_dict["JDAY"],
                                      dsp_dict["HOUR"], dsp_dict["MIN"])
 
    print "Writing plot to file : %s"%fname
    plt.savefig(fname)
    plt.close()

def plot_all_map(dsp_list, coverage=False, period=None, fname=None, legend=True):
    # make this prettier using basemap...
    plt.figure()
    #m.fillcontinents(color='coral', lake_color='aqua')
    # loop over all paths
    for (dsp, dsp_dict) in dsp_list:
        # plot a star at the event location 
        x, y = m(dsp_dict["EVLO"], dsp_dict["EVLA"])
        m.plot(x, y, 'y*', markersize=markersize)
        if coverage:
            # if coverage plot, make all the stations the same color and do not
            # add a legend
            x, y = m(dsp_dict["STLO"], dsp_dict["STLA"])
            m.plot(x, y, 'rv', markersize=markersize)
        else:
            # if not a coverage plot, color code the stations and add a legend
            label = "%s-%s"%(dsp_dict["STA"], dsp_dict["COMP"])
            x, y = m(dsp_dict["STLO"], dsp_dict["STLA"])
            m.plot(x, y, 'v', markersize=markersize,
                 label=label)
        if period is None:
            # if no period is selected, then plot a simple coverage plot
            x_ev, y_ev = m(dsp_dict["EVLO"], dsp_dict["EVLA"])
            x_st, y_st = m(dsp_dict["STLO"], dsp_dict["STLA"])
            m.plot([x_ev, x_st], [y_ev, y_st], 'k')
        else:
            # only plot the lines for this period
            vfunc, efunc = interpolate_dsp(dsp)
            try:
                v = vfunc(period)
                x_ev, y_ev = m(dsp_dict["EVLO"], dsp_dict["EVLA"])
                x_st, y_st = m(dsp_dict["STLO"], dsp_dict["STLA"])
                m.plot([x_ev, x_st], [y_ev, y_st], 'b')
                
            except ValueError:
                # if interpolation did not work, you did not have a
                # measurement, so skip
                pass
    
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(parallels, labels=[1, 1, 0, 0])
    m.drawmeridians(meridians, labels=[0, 0, 0, 1])
    # set filenames
    if coverage:
        if period is None:
            plt.title("Path coverage")
	    if fname is None:
                fname = "coverage.png"
        else:
            plt.title("Group velocity at %03.1f s"%period)
	    if fname is None:
                fname = "coverage_%03.1f.png"%period
    else:
	if legend:
       	    plt.legend(loc='lower left')
        if period is None:
            plt.title("%d(%03d) %02d:%02d"%(dsp_dict["YEAR"], dsp_dict["JDAY"],
                                            dsp_dict["HOUR"], dsp_dict["MIN"]))
	    if fname is None:
                fname = "%d_%03d_%02d_%02d_map.png"%(dsp_dict["YEAR"],
                                                 dsp_dict["JDAY"],
                                                 dsp_dict["HOUR"],
                                                 dsp_dict["MIN"])
        else:
            plt.title("%d(%03d) %02d:%02d @ %03.1fs"%(dsp_dict["YEAR"],
                                                 dsp_dict["JDAY"],
                                                 dsp_dict["HOUR"],
                                                 dsp_dict["MIN"], period))
	    if fname is None:
                fname = "%d_%03d_%02d_%02d_%03.1fs_map.png"%(dsp_dict["YEAR"],
                                                     dsp_dict["JDAY"],
                                                     dsp_dict["HOUR"],
                                                     dsp_dict["MIN"], period)


    print "Writing plot to file : %s"%fname
    plt.savefig(fname)
    plt.close()


if __name__ == '__main__':

    from glob import glob

    print "Treating all dsp files..."
    # get all the filenames
    dsp_files = glob('*.dsp')
    # create a list for all the files
    dsp_list = []

    # loop over the files
    for fname in dsp_files:
        # read the dispersion
        dsp, dsp_dict = read_dsp(fname)
        # add the information to the list
        dsp_list.append((dsp, dsp_dict))
        # plot the single dispersion
        plot_dsp(dsp, dsp_dict)

    # plot the combined dispersion
    plot_all_dsp(dsp_list)
    plot_all_map(dsp_list)

