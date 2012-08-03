#!/bin/bash

#
#       semd_prep.sh
#
#       Transforms SEM output seismogram files into sac files, after convolution
#	with the source function defined by the CMTSOLUTION_with_times file.  Also
#	adjusts the sac header to contain event and station information.
#
#	Expects the following tree structure
#	semd_prep.sh  ---  pac_stations.dat --- model_1 --- model_2 --- model_3 ...
#                                                 |
#                                                 |
#              SPECFEM3D_GLOBE --- event_1 --- event_2 --- event_3 ...
#                |                    |
#  xconvolve_source_timefunction      |
#                            CMTSOLUTION_with_times  ---- semd 
#                                                          |
#                                                          |
#                                                  sta.nw.cmp.semd ... 
#
####################################################################################


event=$1
length_s=$2

here=`pwd`
stations=STATIONS
sfile=stations.dat
dir=${event}
qmxd=${dir}/qmxd
cmt=${dir}/CMTSOLUTION_with_times
cmt_base=${dir}/CMTSOLUTION
conv=${here}/xconvolve_source_timefunction

data=${HOME}/pacific/data
pzdir=${data}/${event}/PZ

do_conv=1	# set =0 if do not want to convolve the source from the
                # CMTSOLUTION file
stype=".true."	# set = .true. for triangle or .false. for gaussian

do_resp=1	# set =0 if do not want to convolve the station response

####################################################################################
#     READ THE CMTSOLUTION FILE FOR EVENT PARAMETERS
####################################################################################

# just make sure solution file exists
if [ ! -e $cmt ] ; then
  echo "Cannot find $cmt.  Will stop now."
  exit
fi

# get information from CMTSOLUTION_with_times file
tshift=`grep "time shift" $cmt | awk '{print $NF}'`
hdur=`grep "half duration" $cmt | awk '{print $NF}'`

cmtlat=`grep "latitude" $cmt | awk '{print $NF}'`
cmtlon=`grep "longitude" $cmt | awk '{print $NF}'`
cmtdep=`grep "depth" $cmt | awk '{print $NF}'`

line=`head -1 $cmt`

ep_year=`echo $line | awk '{print $2}'`
ep_month=`echo $line | awk '{print $3}'`
ep_day=`echo $line | awk '{print $4}'`

ep_jday=`jday $ep_year $ep_month $ep_day | grep "jday" | awk '{print $NF}'`

ep_hour=`echo $line | awk '{print $5}'`
ep_min=`echo $line | awk '{print $6}'`
ep_sec=`echo $line | awk '{print $7}' | awk -F"." '{print $1}'`
ep_hsec=`echo $line | awk '{print $7}' | awk -F"." '{print $2}'`
ep_msec=`echo $ep_hsec | awk '{print $1*10}'`

####################################################################################
#     CALCULATE THE RELEVANT PREM SYNTHETICS
####################################################################################

echo "Finding stations for which there is data"

if [ -e $sfile ] ; then ( rm $sfile ) ; fi

for file in ${data}/${event}/*.BHZ ; do
  sta=`echo $file | awk -F"_" '{print $NF}' | awk -F"." '{print $1}'`  
  net=`echo $file | awk -F"_" '{print $2}'`  
  grep "$sta " $stations | grep " $net " >> $sfile
done
sort $sfile | uniq > t1
mv t1 $sfile

echo "Calculating normal mode synthetics"
calc_modes_qmxd.pl -m $cmt_base -n $length_s -s 1/s -O $qmxd -S $sfile -c LHZ -C
mv RECORDHEADERS $dir/RECORDHEADERS_LHZ
calc_modes_qmxd.pl -m $cmt_base -n $length_s -s 1/s -O $qmxd -S $sfile -c LHR -C
mv RECORDHEADERS $dir/RECORDHEADERS_LHR
calc_modes_qmxd.pl -m $cmt_base -n $length_s -s 1/s -O $qmxd -S $sfile -c LHT -C
mv RECORDHEADERS $dir/RECORDHEADERS_LHT
mv RECORDHEADERS $sfile 

####################################################################################
#     DO SOURCE CONVOLUTION AND TRANSORM TO SAC FORMAT
####################################################################################

if [ ! -e $conv ] ; then
  ln -s ../sem_synth/$conv .
fi

if [ $do_conv -eq 1 ] ; then
  echo "Doing source convolution"
  for file in $qmxd/*.qmxd ; do
    nlines=`wc -l $file | awk '{print $1}'`
    echo $nlines > input_convolve_code.txt
    echo $hdur >> input_convolve_code.txt
    echo $stype >> input_convolve_code.txt
    echo "Convolving $file with hdur = $hdur s"
    $conv < $file > $file.convolved
    ascii2sac.csh $file.convolved
  done
  rm input_convolve_code.txt
else
ascii2sac.csh $qmxd/*.qmxd
fi

mv ${qmxd}/*.sac $dir

####################################################################################
#     SETUP SAC HEADERS, APPLY CENTROID TIME SHIFT AND RENAME FILES
####################################################################################

for sac in $dir/*.sac ; do

  # use filename to extract station information from stations file
  # filname convention :
  # sss.nn.ccc.semd.sac or sss.nn.ccc.semd.convolved.sac

  echo $sac

  file=`echo $sac | awk -F"/" '{print $NF}'`
  sta=`echo $file | awk -F"." '{print $1}'`
  net=`echo $file | awk -F"." '{print $2}'`
  comp=`echo $file | awk -F"." '{print $3}'`

  staline=`grep "$sta " $stations | grep " $net "`
  stlat=`echo $staline | awk '{print $3}'`
  stlon=`echo $staline | awk '{print $4}'`
  stel=`echo $staline | awk '{print $5}'`


  # Write sac macro to change headers
  echo "echo on" > sac.mac
  echo "r $sac" >> sac.mac
  # put location information in first
  echo "ch stla $stlat stlo $stlon stel $stel" >> sac.mac
  echo "ch evla $cmtlat evlo $cmtlon evdp $cmtdep" >> sac.mac
  echo "ch kstnm '$sta'" >> sac.mac
  echo "ch kevnm '$event'" >> sac.mac
  echo "ch kcmpnm '$comp'" >> sac.mac
  echo "w over" >> sac.mac
  echo "r " >> sac.mac
  # first apply time shift so b, e and o refer to the (as yet unset) PDE origin time
  echo "ch allt $tshift" >> sac.mac
  echo "ch o 0.0" >> sac.mac
  # now set the reference time to the epicenter (PDE) origin time
  echo "ch nzyear $ep_year nzjday $ep_jday nzhour $ep_hour nzmin $ep_min nzsec $ep_sec nzmsec $ep_msec" >> sac.mac
  echo "w over" >> sac.mac
  echo "r " >> sac.mac
  echo "ch lovrok true" >> sac.mac
  echo "ch iztype io" >> sac.mac
  echo "w over" >> sac.mac
  echo "quit " >> sac.mac

  sac sac.mac

  mv $sac ${dir}/$event.$file
  sac=${dir}/${event}.${file}

####################################################################################
#     IF REQUIRED, CONVOLVE WITH INSTRUMENT RESPONSE
####################################################################################
 
  if [ $do_resp -eq 1 ] ; then

    dcomp=`echo $comp | sed s/L/B/ `
    pz_base=${pzdir}/SAC_PZs_${net}_${sta}_${dcomp}

    echo "echo on" > resp.mac

    # There may be more than one location code per station - network combination
    # Do convolution for all available PZ files
    ls ${pz_base}* > pz_list
    npz=`wc -l pz_list | awk '{print $1}'`

    n=1
    while [ $n -le $npz ] ; do
      pz=`tail +$n pz_list | head -1`
      loc=`echo $pz | awk -F"_" '{print $NF}'`
      newsac=`echo ${sac}_${loc} | sed s/semd/sem/`
      echo "Creating $newsac"

      echo "r $sac" >> resp.mac 
      echo "rtrend" >> resp.mac
      echo "rmean" >> resp.mac
      echo "taper" >> resp.mac
      echo "transfer from none to polezero s $pz" >> resp.mac      
      echo "w $newsac" >> resp.mac

      n=`echo $n | awk '{print $1+1}'`
    done

    echo "quit" >> resp.mac

    sac resp.mac

  fi

done

rm sac.mac resp.mac

chmod 444 $qmxd/*.qmxd
exit

