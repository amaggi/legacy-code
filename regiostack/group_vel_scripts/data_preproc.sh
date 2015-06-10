#!/bin/bash

#####
#	script for processing IRIS Digital Data
#	modified by Alessia for Yacouba
#####

# read the source information from event_data.txt
source event_data.txt

MIN_FREQ=0.004	# 250s

#####

#####
# extract data and response from all the .seed files present
for f in *.seed
do
echo $f
rdseed -R -d -o 1 -f $f
done
#####

#####
# START LOOP OVER SAC FILES
for i in *.SAC
do
	# clean up any residual response files
	rm -f AMP* PHASE*
	# get information about the file
	NET=`saclhdr -KNETWK $i`
	LOC=`saclhdr -KHOLE $i`
	if [  ${LOC} ]
	then
		if [ ${LOC} = "-12345" ]
		then
			LOC=""
		fi
	fi
	KSTNM=`saclhdr -KSTNM $i`
	KCMPNM=`saclhdr -KCMPNM $i`
	DELTA=`saclhdr -DELTA $i`
	DOY=`saclhdr -NZJDAY $i`
	
	# create response file using evelresp (to be used to transfer to velocity)
	echo RESP.${NET}.${KSTNM}.${LOC}.${KCMPNM} resp
	cp RESP.${NET}.${KSTNM}.${LOC}.${KCMPNM} resp
	evalresp ${KSTNM} ${KCMPNM} ${YEAR} ${DOY} ${MIN_FREQ} 50.0 1000 -u 'vel' -f resp
	mv AMP* afile
	mv PHASE* pfile

	# set upper limits for deconvolution according to sampling interval
	F4=`echo $DELTA | awk '{print 1.0/($1 * 4.0)}'`
	F3=`echo $DELTA | awk '{print 1.0/($1 * 8.0)}'`
	echo $KSTNM $LOC $KCMPNM $F3 $F4

#####
# run sac to do the processing
# step 1 add the event information to the sac header
# step 2 remove the trend
# step 3 deconvolve instrument response
# step 4 low pass filter at the appropriate frequency and decimate to 0.25 s
# step 5 write final sac file with extension .sac
#####
gsac << EOF
r $i
ch EVLA $LAT EVLO $LON EVDP $DEP
ch OGMT $YEAR $DOY $HR $MN $SEC $MSEC
ch lovrok true
ch lcalda true
wh
rtr
transfer from eval subtype  afile pfile TO NONE FREQLIMITS 0.002 ${MIN_FREQ} ${F3} ${F4}
lp c 2 np 3 p 2
interpolate delta 0.25
w ${KSTNM}_${LOC}_${KCMPNM}.sac
quit
EOF

done
# END LOOP OVER SAC FILES
#####

# clean up
rm -f *.SAC RESP.* afile pfile rdseed.err_log resp
