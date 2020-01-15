#!/bin/bash
source utils.sh
sharedir=${AMK}/share
name=$1
inputfile=amk.dat
cwd=$PWD
###reading input
read_input
###
sed 's/thermo/thermo('$temperature','$temperature')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template >  $tsdirll/IRC/minf_${name}.mop
sed 's/thermo/thermo('$temperature','$temperature')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template >  $tsdirll/IRC/minr_${name}.mop
get_geom_mopac.sh $tsdirll/IRC/${name}_ircf.out | awk '{if(NF==4) print $0}' >> $tsdirll/IRC/minf_${name}.mop
get_geom_mopac.sh $tsdirll/IRC/${name}_ircr.out | awk '{if(NF==4) print $0}' >> $tsdirll/IRC/minr_${name}.mop

