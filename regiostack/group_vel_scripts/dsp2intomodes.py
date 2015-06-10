#!/usr/bin/env python
# encoding: utf-8

from glob import glob
from plot_dsp import read_dsp, interpolate_dsp

periods = [20, 50, 100, 150] 
step = 1.0  # smallest grid step for inversion
integration_step = 0.01 # step along path for integrating velocity
random_error = 0.05     # contribution to group velocity error from event
                        # mislocation (in km/s)
corr_distance = 200.     # correlation distance (for inversion)
sigma_m = 0.5       # a-priori model variance

def period2key(period):
    return "%3.1f"%period


if __name__ == '__main__':

    # read the dispersion files
    dsp_files = glob('*/*.dsp')
    dsp_list = []

    for fname in dsp_files:
        # read the dispersion information
        dsp, dsp_dict = read_dsp(fname)
        dsp_list.append((dsp, dsp_dict))

    # create and open one file per period (to be safe)
    intomodes_files = {}
    mean_vel = {}
    n_paths = {}
    for per in periods:
        key = period2key(per)
        fname = "intomodes_%s.dat"%key
        f = open(fname, 'w')
        intomodes_files[key] = f
        mean_vel[key] = 0.0
        n_paths[key] = 0

    # loop over the paths
    for (dsp, dsp_dict) in dsp_list:
        # get the interpolation functions for this dispersion file
        vfunc, efunc = interpolate_dsp(dsp)

        # get interpolated value
        for per in periods:
            try:
                # get the values
                v = vfunc(per)
                e = efunc(per)
                key = period2key(per)
                n_paths[key] += 1
                mean_vel[key] += v
                f = intomodes_files[key]
                f.write("%s\n"%dsp_dict["STA"])
                f.write("%.3f %.3f %.3f %.3f\n" % (dsp_dict["STLA"],
                                                   dsp_dict["STLO"],
                                                   dsp_dict["EVLA"],
                                                   dsp_dict["EVLO"]))
                f.write("%.3f\n" % v)
                f.write("%.3f\n" % e)
            except ValueError:
                # if have a ValueError, then the measurement does not exist
                pass

    # close all files
    for f in intomodes_files.values():
        f.close()

    # loop over periods again to write input files
    for per in periods:
        print per
        key = period2key(per)
        fname = "input_%s.txt"%key
        print fname
        f = open(fname, 'w')
        # get the mean value of the velocity
        mean_vel[key] = mean_vel[key] / n_paths[key]
        f.write("%f\n%f\n%d\n%f\n" % (random_error, integration_step,
                                      n_paths[key], step))
        f.write("1\n1\n1\n")
        f.write("%f\n" % per)
        f.write("intomodes_%s.dat\n" % key)
        f.write("iso_gvel_fin_%s.dat\n" % key)
        f.write("aniso_gvel_fin_%s.dat\n" % key)
        f.write("%.1f %.1f\n" % (corr_distance, corr_distance))
        f.write("0 %.2f %0.2f 0 0 0\n" % (mean_vel[key], sigma_m))
        f.close()

