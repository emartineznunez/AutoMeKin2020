#!/bin/bash
#
source utils.sh

cwd=$PWD
sharedir=${AMK}/share
if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
echo "Input file" $inputfile
exe=$(basename $0)

##Reading HL stuff
read_input
##
##max values 0.001 and 1 
avgerr=$(echo $avgerr | awk '{avg=$1;if(avg>0.001) avg=0.001;print avg}' )
bigerr=$(echo $bigerr | awk '{big=$1;if(big>1) big=1;print big}' )

#On exit remove tmp files
tmp_files=($tsdirhl/IRC/*.chk black* minfailed_list labels mingeom sprint.* deg* tmp* ConnMat ScalMat) 
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT
if [ ! -d "$tsdirhl/PRODs" ]; then
   echo "$tsdirhl/PRODs does not exist. It will be created"
   mkdir $tsdirhl/PRODs
else
   echo "$tsdirhl/PRODs already exists. Remove PR files"
   rm -f $tsdirhl/PRODs/PR*
fi
if [ ! -d "$tsdirhl/MINs/norep" ]; then
   echo "$tsdirhl/MINs/norep does not exist. It will be created"
   mkdir $tsdirhl/MINs/norep
else
   echo "$tsdirhl/MINs/norep already exists."
   rm -r $tsdirhl/MINs/norep
   mkdir $tsdirhl/MINs/norep
fi

#remove some stuff at the beginning (in case of repeating the same script)
rm -f $tsdirhl/MINs/*min*ts*
#
echo "Nomenclature" > $tsdirhl/MINs/names_of_minima
echo "Screening" > $tsdirhl/MINs/minlist_screened

### and prodhl minnrhl table
sqlite3 ${tsdirhl}/PRODs/prodhl.db "drop table if exists prodhl; create table prodhl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, formula TEXT);"
sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "drop table if exists minnrhl; create table minnrhl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, sigma INTEGER);"

echo "List of bad MINIMA" > minfailed_list
echo "PR list" > $tsdirhl/PRODs/PRlist
n=0
nmin=0
npro=0
nrm=0
cp $tsdirhl/min0.log $tsdirhl/IRC
for i in $(ls $tsdirhl/IRC/min*.log)
do 
  ((n=n+1))
  name=$(basename $i .log)
# See if it did not crashed and grab geometires 

  ok=$(awk 'BEGIN{ok=0};/Frequencies/{++nfreq;if($3>0 && $4>0 && nfreq==1) ok=1};END{print ok}' $i)
#Force min0 to be ok
  if [ $ok -eq 0 ] && [ $name = "min0" ]; then ok=1 ; fi

  if [ $ok -eq 1 ]; then
     echo "$name optimized correctly"
# insert all minima except min0
     if [ $name != "min0" ]; then
        zpe=$(get_ZPE_g09.sh $i)
        g=$(get_G_g09.sh $i)
        geom="$(get_geom_g09.sh $i)"
        freq="$(get_freq_g09.sh $i)"
        energy=$(get_energy_g09_$HLcalc.sh $i $noHLcalc)
###EMN
        sqlite3 ${tsdirhl}/MINs/minhl.db "insert or ignore into minhl (natom,name,energy,zpe,g,geom,freq) values ($natom,'$name',$energy,$zpe,$g,'$geom','$freq');"
     fi
#Now we screen the list to rename duplicates
     echo  $natom > mingeom
     echo '' >> mingeom
     get_geom_g09.sh $i >> mingeom
     echo "1" $natom > sprint.dat
     createMat.py mingeom 3 $nA
     cat ConnMat >> sprint.dat
     sprint2.exe <sprint.dat >sprint.out

     paste <(awk 'NF==4{print $1}' mingeom) <(deg.sh) >deg.out
     deg_form.sh > deg_form.out
###
##use absolute energy instead of relative one
     echo $energy >> $tsdirhl/MINs/${name}_data
##
     format.sh $name $tsdirhl/MINs ${nfrag_th} 
     ndis=$(awk '{ndis=$1};END{print ndis}' $tsdirhl/MINs/${name}_data )
### mv MINs where there is 2 or more fragments already formed
     if  [[ ("$ndis" -gt "1") ]] && [[ $name != "min0" ]] 
     then 
        ((npro=npro+1)) 
        echo "Products=" $ndis $name 
        namepr=PR${npro}_${name}
###remove this later on
        sqlite3 "" "attach '${tsdirhl}/MINs/minhl.db' as minhl; attach '${tsdirhl}/PRODs/prodhl.db' as prodhl; insert into prodhl (natom,name,energy,zpe,g,geom,freq) select natom,'$namepr',energy,zpe,g,geom,freq from minhl where name='$name';delete from minhl where name='$name'"
        echo "PROD" $npro $name.rxyz >> $tsdirhl/PRODs/PRlist
     else
        ((nmin=nmin+1))
        echo "min" $nmin "-->" $name.rxyz >> $tsdirhl/MINs/names_of_minima
        echo "min"$nmin "data" >> $tsdirhl/MINs/minlist_screened
        cat $tsdirhl/MINs/${name}_data >> $tsdirhl/MINs/minlist_screened
     fi
###
  else
     echo "$name check later on-->this is a product or problem in opt"
     echo $name.rxyz >> minfailed_list
     continue
  fi
  rm -rf $tsdirhl/MINs/${name}_data

done 
echo "Total number of minima" $nmin

#reduce output
reduce.sh $tsdirhl/MINs min
awk '{if($NF==1) print $0}' $tsdirhl/MINs/minlist_screened.red >  $tsdirhl/MINs/minlist_screened.redconn
awk '{if($NF> 1) print $0}' $tsdirhl/MINs/minlist_screened.red >> $tsdirhl/MINs/minlist_disconnected
diffGT.sh $tsdirhl/MINs/minlist_screened.redconn $tsdirhl/MINs min $avgerr $bigerr
#remove repeated structures in this dir and also in MINs
cp $tsdirhl/MINs/names_of_minima $tsdirhl/MINs/names_of_minima_norep
if [ -f "black_list.out" ]; then 
   for i in $(awk '{print $0}' black_list.out)
   do
     echo "Structure min"$i "repeated"
     ((nrm=nrm+1))
     nomnr="$(awk '{if($2 != '$i') print $0}' $tsdirhl/MINs/names_of_minima_norep)"
     echo "$nomnr" > $tsdirhl/MINs/names_of_minima_norep 
   done
else
   echo "No repetitions"
fi
###
nn=0
for name in $(awk '{if(NR>1) print $4}' ${tsdirhl}/MINs/names_of_minima_norep)
do
  ((nn=nn+1))
  namenrxyz=$(basename $name .rxyz)
  number=$(awk 'NR=='$nn'+1,NR=='$nn'+1{print $2}'  ${tsdirhl}/MINs/names_of_minima_norep)
  namenr=$(basename min${number}_${name} .rxyz)
##insert data into minnrhl.db from minhl.db
#  cp ${tsdirhl}/MINs/${name} ${tsdirhl}/MINs/norep/min${number}_${name}
##
  sqlite3 "" "attach '${tsdirhl}/MINs/minhl.db' as minhl; attach '${tsdirhl}/MINs/norep/minnrhl.db' as minnrhl; insert into minnrhl (natom,name,energy,zpe,g,geom,freq,sigma) select natom,'$namenr',energy,zpe,g,geom,freq,sigma from minhl where name='$namenrxyz';"
done


###
((nfin=nmin-nrm))
echo $nmin "have been optimized at the HL, of which" $nrm "removed because of repetitions"
echo "Current number of MINs optimized at the HL=" $nfin
#


##################
##################

echo "Now running MINFAILED.sh"

##################
##################

cp $tsdirhl/PRODs/PRlist $tsdirhl/PRODs/PRlist.old
npro=$(awk '{npro=$2};END{print npro}' $tsdirhl/PRODs/PRlist )
#file to screen the failed minima and/or products

nfail=$(wc -l minfailed_list | awk '{print $1-1}')
if [ $nfail -eq 0 ]; then
   echo "You are lucky. All minima have been optimized correctly. Exit here"
   exit
else
   echo "number of minima that failed and/or are products" $nfail
fi

for i in $(awk 'NR>1{print $0}' minfailed_list)
do 
  name=$(basename $i .rxyz)
  geom="$(get_geom_g09.sh ${tsdirhl}/IRC/${name}.log)"
#Now we screen the list to rename duplicates
  echo $natom > mingeom
  echo '' >> mingeom
  get_geom_g09.sh ${tsdirhl}/IRC/${name}.log >>mingeom
##If the calc. was not done skip this minimum
  anlf=$(wc -l mingeom | awk '{print $1}')
  nlmg=$(($natom+2))
  if [ $anlf -lt $nlmg ]; then
     echo "Double check this opt: $name"
     continue
  fi
##
  echo "1" $natom > sprint.dat
  createMat.py mingeom 3 $nA
  cat ConnMat >> sprint.dat
  sprint2.exe <sprint.dat >sprint.out
  paste <(awk 'NF==4{print $1}' mingeom) <(deg.sh) >deg.out
  deg_form.sh > deg_form.out
##
  echo "This is a just to see if there is more than one fragment" > $tsdirhl/MINs/${name}_data
#
  format.sh $name $tsdirhl/MINs ${nfrag_th} 
  ndis=$(awk '{ndis=$1};END{print ndis}' $tsdirhl/MINs/${name}_data )
#  echo $name $ndis
### mv MINs where there is 2 or more fragments already formed
  if  [[ ("$ndis" -gt "1") ]]
  then 
     ((npro=npro+1)) 
##remove this later on
     namepr="PR"$npro"_"$name
##EMNinsert into prodhl.db
     sqlite3 ${tsdirhl}/PRODs/prodhl.db "insert into prodhl (natom,name,energy,zpe,g,geom,freq) values ($natom,'$namepr',0,0,0,'$geom',0);"
     lp=$(awk 'BEGIN{lp=1};/'$name'/{lp=0};END{print lp}' $tsdirhl/PRODs/PRlist )
     if  [[ ("$lp" -eq "1") ]]
     then
         echo "PROD" $npro $name.rxyz >> $tsdirhl/PRODs/PRlist
         echo "Move structure $name to PROD # $npro"
     else
         echo "Structure $name already moved to PROD" 
     fi
  else
     echo "Double check this opt: $name"
  fi
done 
