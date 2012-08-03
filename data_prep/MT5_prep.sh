#!/software/bin/bash

#	MT5_prep.sh
#
#	This script takes directories output from SEED_prep.sh and 
#	preprocesses the SAC files so they can be reduced to MT5
#	format using sac2mt5.  
#	This script assumes you want to deconvolve the station
#	responses from the seismograms and turn them into wwlpbn
#	traces.

#	The script assumes the following directory stucture:
#
#	yyyyjjjhhmmss
#		yyyyjjjhhmmss.txt
#		yyyyjjjhhmmss_NET_STA.CHA (for a number of stations)
#		RESP
#			RESP.NET.STA..BHZ (where the extra dot is for loc code)
#		PZ
#			SAC_PZs_NET_STA_CHA
#
#	The script needs the environment variable DATA_PREP to be set 
#	to the same value as it is set for SEED_prep.

#	Usage:
#		MT5_prep.sh yyyyjjjhhmmss

dirname=$1

# Convert all files to wwlpbn
cd $dirname

echo echo on > sac.mac
for i in *_*_*.??? ; do

	sta=`echo $i | awk -F"_" '{ print $3}' | awk -F"." '{print $1}'`
	cha=`echo $i | awk -F"_" '{ print $3}' | awk -F"." '{print $2}'`
	net=`echo $i | awk -F"_" '{ print $2}'` 
	pz=`ls PZ/SAC_PZs* | grep $net | grep $sta | grep $cha `

	echo r $i >> sac.mac
	echo dec 5 >> sac.mac
	echo dec 4 >> sac.mac
	echo w append _wwlpbn >> sac.mac
	echo r ${i}_wwlpbn >> sac.mac
	echo trans from polezero sub $pz to wwlpbn >> sac.mac
	echo w over >> sac.mac

done
echo quit >> sac.mac

#/applications/local/bin/sac2000 sac.mac
/applications/local/bin/sac2000 sac.mac

# Keep a copy of the original broad band files for time picking
mkdir ${dirname}_bb
mv *_*_*.??? ${dirname}_bb

# Pick the P and S times out of the seismograms and write
# sac macros to cut around them

echo echo on > cut_p.mac
for i in *_*_*.??Z_wwlpbn ; do
	index=`/home/maggi/bin/b-612/pick_index $i P`
	echo P is at index $index in file $i
        echo $index | awk '{if($1 != "Cannot" ) printf "cut t%d -30 90\nr %s\ncut off\nch kcmpnm LHZ\nw append _cut\n", $1 , file}' file=$i >> cut_p.mac
#        echo $index | awk '{print "cut t" $1, "-30 90"}' >> cut_p.mac
#	echo r $i >> cut_p.mac
#	echo cut off >> cut_p.mac
#	echo ch kcmpnm LHZ >> cut_p.mac
#	echo w append _cut >> cut_p.mac
done
echo quit >> cut_p.mac

echo echo on > cut_s.mac
for i in *_*_*.??E_wwlpbn ; do
	index=`/home/maggi/bin/b-612/pick_index $i S`
	echo S is at index $index in file $i
        echo $index | awk '{if ($1 != "Cannot")  printf "cut t%d -30 90\nr %s\ncut off\nch kcmpnm LHE\nw append _cut\n", $1 , file}' file=$i >> cut_s.mac
#        echo $index | awk '{print "cut t" $1, "-30 90"}' >> cut_s.mac
#	echo r $i >> cut_s.mac
#	echo cut off >> cut_s.mac
#	echo ch kcmpnm LHE >> cut_s.mac
#	echo w append _cut >> cut_s.mac
done
for i in *_*_*.??N_wwlpbn ; do
	index=`/home/maggi/bin/b-612/pick_index $i S`
	echo S is at index $index in file $i
        echo $index | awk '{if ($1 != "Cannot") printf "cut t%d -30 90\nr %s\ncut off\nch kcmpnm LHN\nw append _cut\n", $1 , file}' file=$i >> cut_s.mac
#        echo $index | awk '{print "cut t" $1, "-30 90"}' >> cut_s.mac
#	echo r $i >> cut_s.mac
#	echo cut off >> cut_s.mac
#	echo ch kcmpnm LHN >> cut_s.mac
#	echo w append _cut >> cut_s.mac
done
echo quit >> cut_s.mac

# run the sac macros
/applications/local/bin/sac2000 cut_p.mac
/applications/local/bin/sac2000 cut_s.mac

# clean up
mkdir ${dirname}_wwlpbn
mv *_*_*.???_wwlpbn ${dirname}_wwlpbn

mkdir ${dirname}_mt5
mv *_wwlpbn_cut ${dirname}_mt5

rm [a-z]*

#gzip -r ${dirname}_bb
#gzip -r ${dirname}_wwlpbn
