#!/bin/sh


DATA_PREP=${CODE}/src/data_prep


dirname=$1


cd $dirname

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
