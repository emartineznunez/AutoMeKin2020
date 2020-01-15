#!/bin/bash
#
source utils.sh

cwd=$PWD
inputfile=amk.dat
exe=$(basename $0)

###reading input file
read_input
###
##Create table for ts rxyz files
sqlite3 ${tsdirll}/TSs/ts.db "drop table if exists ts; create table ts (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT,sigma INTEGER);"
##

#remove tmp files
tmp_files=($tsdirll/IRC/*.arc $tsdirll/IRC/*.mop $tsdirll/TSs/*.mop)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT
# first of all, see if all ircs are correct


# analysis
for name in $(awk '{print $3}' $tslistll)
do
  echo "Analyzing $name"

  if [ -f ${tsdirll}/ts_mopac_failed ]; then
    skip=$(awk 'BEGIN{skip=0};/'$name'/{skip=1};END{print skip}' ${tsdirll}/ts_mopac_failed)
    if [ $skip == 1 ]; then
       echo "TS $name has ben previously discarded-->(skip because mopac cannot optimize it)"
       continue
    fi
  fi

#First the tss
  geom="$(get_geom_thermo_mopac.sh $tsdirll/TSs/${name}_thermo.out)"
  if [ "$program_opt" != "mopac" ]; then
     freq="$(get_freq_mopac.sh $tsdirll/${name}_mop.out)"
  else
     freq="$(get_freq_mopac.sh $tsdirll/${name}.out)"
  fi
  e=$(awk 'BEGIN{e=0};/HEAT OF FORMATION =/{e=$5};END{print e}' $tsdirll/TSs/${name}_thermo.out )
  zpe=$(awk 'BEGIN{zpe=0};/          ZERO POINT ENERGY/{zpe=$4};END{print zpe}' $tsdirll/TSs/${name}_thermo.out )
  sigma=$(awk '/SYMMETRY NUMBER/{print $NF;exit}' $tsdirll/TSs/${name}_thermo.out)
  g_corr=$(awk 'BEGIN{zpe=0;h=0;s=0;t='$temperature'}
     /          ZERO POINT ENERGY/{zpe=$4}
     /CALCULATED THERMODYNAMIC PROPERTIES/{ok=1}
  {if(ok==1 && $1 == '$temperature') {
  getline;getline;getline;getline;
  h=$3/1000;s=$5/1000;print zpe+h-t*s;exit}
  }' $tsdirll/TSs/${name}_thermo.out )
##insert into ts.db
  sqlite3 ${tsdirll}/TSs/ts.db "insert into ts (natom,name,energy,zpe,g,geom,freq,sigma) values ($natom,'$name',$e,$zpe,$g_corr,'$geom','$freq',$sigma);"

#Now the minima
  namef=minf_${name}
  geomf="$(get_geom_thermo_mopac.sh $tsdirll/IRC/${namef}.out)"
  freqf="$(get_freq_mopac.sh $tsdirll/IRC/${namef}.out)"
  ef=$(awk 'BEGIN{e=0};/HEAT OF FORMATION =/{e=$5};END{print e}' $tsdirll/IRC/${namef}.out )
  zpef=$(awk 'BEGIN{zpe=0};/          ZERO POINT ENERGY/{zpe=$4};END{print zpe}' $tsdirll/IRC/${namef}.out )
  sigmaf=$(awk '/SYMMETRY NUMBER/{print $NF;exit}' $tsdirll/IRC/${namef}.out)
  g_corrf=$(awk 'BEGIN{zpe=0;h=0;s=0;t='$temperature'}
     /          ZERO POINT ENERGY/{zpe=$4}
     /CALCULATED THERMODYNAMIC PROPERTIES/{ok=1}
  {if(ok==1 && $1 == '$temperature') {
  getline;getline;getline;getline;
  h=$3/1000;s=$5/1000;print zpe+h-t*s;exit}
  }' $tsdirll/IRC/${namef}.out )
#The minima might have failed in the opt process. In that case empty the freq column 
  if [ -z "$freqf" ]; then
     echo "Problems with this minimum: $namef-->thermo calc failed"
     zpef=0
     g_corrf=0
     ircn0=$(echo $namef | sed 's@min@@;s@_@ @') 
     ircnf=$(echo $ircn0 | awk '{print $2"_irc"$1".xyz"}')
     ef=$(awk '/HEAT OF FORMATION/{heat=$(NF-4)};END{print heat}' $tsdirll/IRC/$ircnf)
     geomf="$(awk '/HEAT OF FORMATION/{natom=0};{if(NF==4) {++natom;line[natom]=$0} };END{i=1;while(i<=natom){print line[i];i++}}' $tsdirll/IRC/$ircnf)"
     freqf=""
     sigmaf=1
  fi
##insert into min.db
  sqlite3 ${tsdirll}/MINs/min.db "insert or ignore into min (natom,name,energy,zpe,g,geom,freq,sigma) values ($natom,'$namef',$ef,$zpef,$g_corrf,'$geomf','$freqf',$sigmaf);"

  namer=minr_${name}
  geomr="$(get_geom_thermo_mopac.sh $tsdirll/IRC/${namer}.out)"
  freqr="$(get_freq_mopac.sh $tsdirll/IRC/${namer}.out)"
  er=$(awk 'BEGIN{e=0};/HEAT OF FORMATION =/{e=$5};END{print e}' $tsdirll/IRC/${namer}.out )
  zper=$(awk 'BEGIN{zpe=0};/          ZERO POINT ENERGY/{zpe=$4};END{print zpe}' $tsdirll/IRC/${namer}.out )
  sigmar=$(awk '/SYMMETRY NUMBER/{print $NF;exit}' $tsdirll/IRC/${namer}.out)
  g_corrr=$(awk 'BEGIN{zpe=0;h=0;s=0;t='$temperature'}
     /          ZERO POINT ENERGY/{zpe=$4}
     /CALCULATED THERMODYNAMIC PROPERTIES/{ok=1}
  {if(ok==1 && $1 == '$temperature') {
  getline;getline;getline;getline;
  h=$3/1000;s=$5/1000;print zpe+h-t*s;exit}
  }' $tsdirll/IRC/${namer}.out )
#The minima might have failed in the opt process. In that case empty the freq column 
  if [ -z "$freqr" ]; then
     echo "Problems with this minimum: $namer-->thermo calc failed"
     er=0
     zper=0
     g_corrr=0
     ircn0=$(echo $namer | sed 's@min@@;s@_@ @')
     ircnr=$(echo $ircn0 | awk '{print $2"_irc"$1".xyz"}')
     er=$(awk '/HEAT OF FORMATION/{heat=$(NF-4)};END{print heat}' $tsdirll/IRC/$ircnr)
     geomr="$(awk '/HEAT OF FORMATION/{natom=0};{if(NF==4) {++natom;line[natom]=$0} };END{i=1;while(i<=natom){print line[i];i++}}' $tsdirll/IRC/$ircnr)"
     freqr=""
     sigmar=1
  fi
##insert into min.db
  sqlite3 ${tsdirll}/MINs/min.db "insert or ignore into min (natom,name,energy,zpe,g,geom,freq,sigma) values ($natom,'$namer',$er,$zper,$g_corrr,'$geomr','$freqr',$sigmar);"
done
