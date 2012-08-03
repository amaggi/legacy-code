#!/bin/bash

#
#       $Id: data_prep.sh,v 1.2 2005/02/02 01:54:42 alessia Exp $
#
#	This script is a quick and dirty way of doing all that the program
#	data_prep does, with very little error checking, and only using the
#	station file if absolutelly necessary.
#
#
#	Arguments: 1) directory name 2) length of record desired in seconds
#	Optional arguments 3) non zero if want to use cmt lat, lon info
#

dname=$1
time=$2

#	Read event information from .txt file.
#	Format: first line (ignored)
#	second line: 
#	year (jdy) mm dd hh mm ss.ss  lat lon dep mb 

tail -1 $dname.txt > temp
year=`awk '{print $1}' < temp `
jday=`awk '{print substr($2,2,3)}' < temp `
month=`awk '{print $3}' < temp`
day=`awk '{print $4}' < temp`
hour=`awk '{print $5}' < temp`
min=`awk '{print $6}' < temp`
sec=`awk '{print $7}' < temp | awk -F. '{print $1}' `
hsec=`awk '{print $7}' < temp | awk -F. '{print $2}' `
evlat=`awk '{print $8}' < temp`
evlon=`awk '{print $9}' < temp`
evdep=`awk '{print $10}' < temp`
mb=`awk '{print $11}' < temp`

rm -f temp


#	Add this information to all SAC headers
#	by writing a sac macro to take care of things

echo "echo on" > sac_mac.m
echo r "$"1  >> sac_mac.m
echo rmean >> sac_mac.m
echo rtrend >> sac_mac.m
echo ch o >> sac_mac.m
echo ch T0 >> sac_mac.m
echo ch T1 >> sac_mac.m
echo ch T2 >> sac_mac.m
echo ch T3 >> sac_mac.m
echo ch T4 >> sac_mac.m
echo ch T5 >> sac_mac.m
echo ch T6 >> sac_mac.m
echo ch T7 >> sac_mac.m
echo ch T8 >> sac_mac.m
echo ch T9 >> sac_mac.m
echo w over >> sac_mac.m
echo r "$"1  >> sac_mac.m
echo ch evla $evlat evlo $evlon evdp $evdep >> sac_mac.m
echo "evaluate to msec" $hsec "* 100" >> sac_mac.m
echo ch o gmt $year $jday $hour $min $sec "%msec"  >> sac_mac.m
echo w over >> sac_mac.m
echo r "$"1 >> sac_mac.m
echo setbb otime "&1,o" >> sac_mac.m
echo "evaluate to tshift -1 * %otime" >> sac_mac.m
echo ch allt %tshift >> sac_mac.m
echo ch user7 $mb >> sac_mac.m
echo ch lovrok true >> sac_mac.m
echo ch iztype io >> sac_mac.m
echo w over >> sac_mac.m
echo cut o 0 $time >> sac_mac.m
# use cuterr usebe if you are sure that at least some of the data
# lies within the cut window, and you don't like filling with zeros
# if some seismograms may contain data wholly without of your data
# window then you must use cuterr fillz
echo cuterr usebe >> sac_mac.m
#echo cuterr fillz >> sac_mac.m
echo r >> sac_mac.m
echo cut off >> sac_mac.m
echo rmean >> sac_mac.m
echo w over >> sac_mac.m
echo quit >> sac_mac.m

