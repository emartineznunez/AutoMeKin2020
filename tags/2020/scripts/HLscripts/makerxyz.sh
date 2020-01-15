#!/bin/bash

source utils.sh
#remove tmp files
tmp_files=(ee* freq* gcorr* geom* zpe* tmp*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
cwd=$PWD
inputfile=amk.dat

##
read_input
##

working=$tsdirhl/PRODs/CALC/working
sqlite3 ${tsdirhl}/PRODs/CALC/prodfhl.db "drop table if exists prodfhl; create table prodfhl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT,formula TEXT );"
ext=0
for line in $(awk 'NR>1{print NR}' ${tsdirhl}/PRODs/PRlist_frag)
do
    ((ext=ext+1))
    name="PR"$(awk 'NR=='$line'{print $2"_"$3}' ${tsdirhl}/PRODs/PRlist_frag)
    echo "Doing calcs for prod # $ext"
    namedb=$(basename $name .rxyz)
    rm -rf ee$ext freq$ext zpe$ext gcorr$ext geom$ext geom0$ext
    nfrag=$(awk 'NR=='$line'{for(i=4;i<=NF;i++) if($i!="+")++nf};END{print nf}' ${tsdirhl}/PRODs/PRlist_frag)
    for j in $(seq $nfrag)
    do
       frag=$(awk 'NR=='$line'{for(i=4;i<=NF;i++) if($i!="+") {++nf;fr[nf]=$i}};END{print fr['$j']}' ${tsdirhl}/PRODs/PRlist_frag)
       file0=${tsdirhl}/PRODs/CALC/${frag}.log
#EMN. If it is repeated, look for the original file
       if [ ! -f $file0 ]; then
          file0a=$(basename $file0 .log)
          file0b=$(cat $working/fraglist | awk '{if($2=="'$file0a'") {print $3;exit}}' ) 
          file=${tsdirhl}/PRODs/CALC/${file0b}.log
       else
          file=$file0
       fi
#EMN
       get_energy_g09_$HLcalc.sh $file $noHLcalc >> ee$ext 
       get_geom_g09.sh $file > geom0$ext
       awk '{print $1,$2+5*'$j',$3,$4}' geom0$ext >> geom$ext
       get_ZPE_g09.sh $file >> zpe$ext 
       get_freq_g09.sh $file >> freq$ext
#EMN
##create temp file tmp_rxyz from sqlite tables
       cat geom0$ext | wc -l  >tmp_rxyz 
       sigma=$(awk 'BEGIN{IGNORECASE=1};/SYMMETRY NUMBER/{print $NF;exit}' $file | sed 's@\.@@' )
       echo E= $(get_energy_g09_$HLcalc.sh $file $noHLcalc) zpe= $(get_ZPE_g09.sh $file) g_corr= 0 sigma= $sigma >>tmp_rxyz
       cat geom0$ext >> tmp_rxyz
       get_freq_g09.sh $file | awk '{print sqrt($1*$1)}'  >> tmp_rxyz
###Calculate g using saulo's thermochem.py 
#      get_G_g09.sh $file >> gcorr$ext
       mult="$(awk '/Multiplicity/{print $NF}' $file)"
       thermochem.py tmp_rxyz $temperature hl $mult | awk '/Thermal correction to Gib/{print $9}' >> gcorr$ext
    done
    e=$(awk '{e+=$1};END{printf "%14.9f\n",e}' ee$ext)
    gcorr=$(awk '{e+=$1};END{printf "%14.9f\n",e}' gcorr$ext)
    zpe=$(awk '{e+=$1};END{printf "%8.2f\n",e}' zpe$ext)
    natom=$(awk '{if(NF==4) ++natom};END{print natom}' geom$ext)
    geom=$(cat geom$ext)
    freq=$(cat freq$ext)
    formula="$(sqlite3 $tsdirhl/PRODs/prodhl.db "select formula from prodhl where name='$namedb'")"
###insert into prodfhl
    sqlite3 ${tsdirhl}/PRODs/CALC/prodfhl.db "insert into prodfhl (natom,name,energy,zpe,g,geom,freq,formula) values ($natom,'$namedb',$e,$zpe,$gcorr,'$geom','$freq','$formula');"
done

echo "End of calcs"
