#!/bin/sh

# Earthquake focal parameters
strike=30           # degrees
dip=50
rake=90
moment=1e27         # dyn/cm
source_rise_time=1  # in seconds
source_length=10    # in seconds

# Earthquake location
evlat=0.
evlon=10.
evdp=9

# Station location
stlat=0.
stlon=30.

# Length of seismogram
npts=2048
dt=1

# Output name for seismogram .vrt .rad .trn
outputname=seis

# Response type (might need to change...)
resp=ida

# This is where it all happens:
g=gram.input
echo 2 > $g
echo >> $g
echo 1990 01  1 00 00.0 $evlat $evlon $evdp >> $g
echo $strike $dip $rake 0 0 0 0 0 0 $moment $source_rise_time $source_length >> $g
echo >> $g
echo $resp >> $g
echo $stlat $stlon >> $g
echo n >> $g  # want to modify some of the defaults
echo 3 >> $g  # parameter 3 is number of points
echo $npts >> $g
echo 4 >> $g  # parameter 4 is dt
echo $dt >> $g
echo 7 >> $g  # parameter 7 is start time of seismogram (default 40)
echo -40 >> $g
echo 0 >> $g  # we are happy with the parameters now
echo e >> $g
echo rayleigh >> $g
echo love >> $g
echo 0 >> $g  # we use all the modes in the earth-created files
echo $outputname >> $g
echo y >> $g # y=no straingrams, n=straingrams
echo n >> $g # program will segfault because of bad file handling!



gram < $g

