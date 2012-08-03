#!/bin/sh

#	SEED_prep.sh
#
#	This script fills in the header information and cuts the
#	time series for IRIS/PASSCAL broadband data.  The first 
#	argument is the directory name where the data is and the 
#	second argument is the length of the time series in seconds.
#	The third optional argument can be used to do further analysis
#	
#	It needs the following environment variable to be set:
#
#	DATA_PREP ->	the directory containing the IASPEI related files 
#			necessary for finding travel times (e.g. iasp91.hed
#			iasp91.tbl)



DATA_PREP=${HOME}/code/share

#echo ${DATA_PREP}

dirname=$1
time=$2
cmt=1

# if cmt != 0 then 

#	Uncompress all files in the directory
#	Add lines here if the compression algorithm was anything
#	other than gzip.

gunzip -r $dirname

#	Change into the directory
#	The directory names are in the following convention:
#	YYYYJJJHHMMSS

cd $dirname

#	This next group of statements change from the IRIS file
#	naming convention to the Cambridge naming convention.
#	YYYYJJJHHMMSS_NET_STA.CHA


#	The following is written assuming you have used version 4.17+ of
#	rdseed, which also contains location codes.
#	Gets the network name, the station name and the channel from the
#	rdseed outputted sac filename, and renames the sac file 
#	according to the Cambridge convention

for sacfile in *.SAC ; do
	network=`echo $sacfile | awk -F. '{print $7}'`
	station=`echo $sacfile | awk -F. '{print $8}'`
	loc=`echo $sacfile | awk -F. '{print $9}'`
	channel=`echo $sacfile | awk -F. '{print $10}'`


#	If there is a gap in the seismogram in the SEED file, and
#	this gap is greater than the tollerence set by the 
#	XXXXXXXXX environment variable, then rdseed will make separate
#	files for the bits either side of the gap.  
#	Here we fine out how many files we have for each 
#	network-station-channel grouping.
#	If there are many files per seismogram 
#	then there probably was a timing problem.  You could either 
#	decide to get rid of the seismogram entirely or increase 
# 	the rdseed tollerence set by XXXXXXXX.
#	This script will move multiple files to a separate directory
#	so you can decide what to do with them later.

	nfiles=`ls | grep ${network}.${station}.${loc} | grep $channel | wc | awk '{print $1}'`
#	echo $network $station $channel $nfiles
	if  [ $nfiles -gt 1 ] ; then 
		mkdir Multiple_Seismograms
		mv *${network}.${station}.${loc}.${channel}.SAC Multiple_Seismograms ;
	fi

#	Rename the sac files
	mv $sacfile ${dirname}_${network}_${loc}_${station}.${channel}
done


#	This next group of statements use the source parameters
#	from the pde file to fill in the header of each SAC file.
#	Each SAC file is cut from the origin for the requested
#	number of seconds.
#
#	This is a silly way of doing things, as a do loop in the 
#	sac macro would be much faster, but the darned thing does not
#	work on all systems (yet!), so we have to go with the old fashioned way.

echo $dirname
#
data_prep.sh $dirname $time 
for file in *.?[L,H,V]? ; do
	sac2000 sac_mac.m $file
done
#
#	This next group of statements add the IASPIE P-wave and 
#	S-wave arrival times to t0 and t1 in the SAC header.
#
#	First set up the relevent links
ln -s $DATA_PREP/iasp91.hed iasp91.hed
ln -s $DATA_PREP/iasp91.tbl iasp91.tbl
#
#	Now for each seismogram write the relevent get_iasp.m files
#	and run them
for file in *.?[L,V,H]? ; do
	echo echo on > get_iasp.m
	echo rh $file >> get_iasp.m
	echo setbb dp "&1,evdp">> get_iasp.m
	echo setbb del "&1,gcarc">> get_iasp.m
	echo "$"run ttimes >> get_iasp.m
	echo P >> get_iasp.m
	echo S >> get_iasp.m
	echo >> get_iasp.m
	echo %dp >> get_iasp.m
	echo %del >> get_iasp.m
	echo -10 >> get_iasp.m
	echo -10 >> get_iasp.m
	echo "$"endrun>> get_iasp.m
	echo quit >> get_iasp.m

	sac2000 get_iasp.m $file

	# Now all the relevent travetime information is in the dat file
	# Findit reads the information from dat and rewrites it to
	# dat1, with -12345 as the non-defined value

	findit 

	# Now we have to read the information
	# Apologies for the stupid way this is being done:
	# it is FAR too late in the day for intelligence!

	id1=`cat dat1 | awk '{ if ( FNR == 1 ) print $1}'` 
	id2=`cat dat1 | awk '{ if ( FNR == 2 ) print $1}'`
	id3=`cat dat1 | awk '{ if ( FNR == 3 ) print $1}'`
	id4=`cat dat1 | awk '{ if ( FNR == 4 ) print $1}'`
	id5=`cat dat1 | awk '{ if ( FNR == 5 ) print $1}'`
	id6=`cat dat1 | awk '{ if ( FNR == 6 ) print $1}'`
	id7=`cat dat1 | awk '{ if ( FNR == 7 ) print $1}'`
	id8=`cat dat1 | awk '{ if ( FNR == 8 ) print $1}'`
	tim1=`cat dat1 | awk '{ if ( FNR == 9 ) print $1}'`
	tim2=`cat dat1 | awk '{ if ( FNR == 10 ) print $1}'`
	tim3=`cat dat1 | awk '{ if ( FNR == 11 ) print $1}'`
	tim4=`cat dat1 | awk '{ if ( FNR == 12 ) print $1}'`
	tim5=`cat dat1 | awk '{ if ( FNR == 13 ) print $1}'`
	tim6=`cat dat1 | awk '{ if ( FNR == 14 ) print $1}'`
	tim7=`cat dat1 | awk '{ if ( FNR == 15 ) print $1}'`
	tim8=`cat dat1 | awk '{ if ( FNR == 16 ) print $1}'`
	pslw=`cat dat1 | awk '{ if ( FNR == 17 ) print $1}'`
	sslw=`cat dat1 | awk '{ if ( FNR == 18 ) print $1}'`
	ppslw=`cat dat1 | awk '{ if ( FNR == 19 ) print $1}'`

	echo echo on > get_iasp.m
	echo rh $file >> get_iasp.m
	echo ch t1 $tim1 kt1 $id1 >> get_iasp.m
	echo ch t2 $tim2 kt2 $id2 >> get_iasp.m
	echo ch t3 $tim3 kt3 $id3 >> get_iasp.m
	echo ch t4 $tim4 kt4 $id4 >> get_iasp.m
	echo ch t5 $tim5 kt5 $id5 >> get_iasp.m
	echo ch t6 $tim6 kt6 $id6 >> get_iasp.m
	echo ch t7 $tim7 kt7 $id7 >> get_iasp.m
	echo ch t8 $tim8 kt8 $id8 >> get_iasp.m
	echo evaluate to NEWP $pslow / 111.125 >> get_iasp.m
	echo ch user0 %NEWP kuser0 'pslow' >> get_iasp.m
	echo wh over >> get_iasp.m
	echo quit >> get_iasp.m

	sac2000 get_iasp.m $file

	rm -f dat1 dat 
done
# The end!
