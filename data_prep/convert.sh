#!/software/bin/bash

# This script converts files which are the 
# output of SEED_prep.sh to displacement
# ready for input into Newton5. It also 
# creates ascii archives of the converted
# files.
dir=$1

cd $dir

echo echo on > sac.mac
for i in *_*_*.??Z ; do

	sta=`echo $i | awk -F"_" '{ print $3}' | awk -F"." '{print $1}'`
	cha=`echo $i | awk -F"_" '{ print $3}' | awk -F"." '{print $2}'`
	net=`echo $i | awk -F"_" '{ print $2}'` 
	pz=`ls PZ/SAC_PZs* | grep $net | grep $sta | grep $cha `

	echo r $i >> sac.mac
	echo trans from polezero sub $pz to none freq 0.001 0.005 1.0 2.0 >> sac.mac
	echo w over >> sac.mac
	echo convert from sac $i to alpha over >> sac.mac


done
echo quit >> sac.mac

/applications/local/bin/sac2000 sac.mac

rm -rf PZ RESP
rm -rf [a-z]*

cd ..

tar cvzf $dir.tgz $dir
rm -rf $dir
