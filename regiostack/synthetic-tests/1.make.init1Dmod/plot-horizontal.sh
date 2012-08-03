#!/bin/sh

infile=$1
out=`echo $infile | sed s/.xyz/.eps/`
grdfile=`echo $infile | sed s/.xyz/.grd/`

gmtset BASEMAP_TYPE plain MEASURE_UNIT inch PAPER_MEDIA a4+

region="121/301/-89/89"
proj="X5d"
ticks="15g/15g"

cpt=vs.cpt

xyz2grd $infile -: -G$grdfile -R$region -I2/2 -V

psbasemap -R$region -J$proj -B$ticks -K -V -P -Y2 > $out
grdimage $grdfile -R$region -J$proj -C$cpt -O -K >> $out
pscoast -R$region -J$proj -Dl -W2 -A5000/0 -O -K >> $out
psscale -D2.5/-0.5/5/0.25h -C$cpt -B1 -O -K >> $out
psbasemap -R$region -J$proj -B$ticks -O -V   >> $out

gv $out

