#!/bin/bash
source utils.sh
#remove tmp files
tmp_files=(tmp*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
inputfile=amk.dat
cwd=$PWD

###reading input files
read_input
###

# rate=0 print DG
# rate=1 print DE
if [ $rate -eq 0 ] ; then
   echo "Sort Gibbs Free Energy differences"
elif [ $rate -eq 1 ]; then
   echo "Sort Energy(+ZPE) differences"
else
   echo "Specify a type of rate: canonical or microcanonical"
   exit
fi

if [ $rate -eq 1 ]; then
   dum_min0=$(sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "select energy,zpe from minnrhl where name like '%min0%'" | sed 's@|@ @g')
elif [ $rate -eq 0 ]; then
   dum_min0=$(sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "select energy,g from minnrhl where name like '%min0%'" | sed 's@|@ @g')
fi

# First of all, we sort the TSs
if [ ! -d "$tsdirhl/TSs/SORTED" ]; then
   echo "$tsdirhl/TSs/SORTED does not exist. It will be created"
   mkdir $tsdirhl/TSs/SORTED
else
   rm -r $tsdirhl/TSs/SORTED
   mkdir $tsdirhl/TSs/SORTED
   echo "$tsdirhl/TSs/SORTED already exists but it was created again"
fi
if [ ! -d "$tsdirhl/MINs/SORTED" ]; then
   echo "$tsdirhl/MINs/SORTED does not exist. It will be created"
   mkdir $tsdirhl/MINs/SORTED
else
   rm -r $tsdirhl/MINs/SORTED
   mkdir $tsdirhl/MINs/SORTED
   echo "$tsdirhl/MINs/SORTED already exists but it was created again"
fi
###create sorted tables
sqlite3 ${tsdirhl}/TSs/SORTED/tsshl.db "drop table if exists tsshl; create table tsshl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,lname TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT);"
sqlite3 ${tsdirhl}/MINs/SORTED/minshl.db "drop table if exists minshl; create table minshl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,lname TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT);"

# We start with the TSs
if [ -f tmp_oren ]; then rm tmp_oren ; fi
if [ -f tmp_nonoren ]; then rm tmp_nonoren ; fi


echo "non-ordered energies"
for name in $(sqlite3 ${tsdirhl}/TSs/tshl.db "select name from tshl")
#for i in $(ls $tsdirhl/TSs/ts*.rxyz)
do
  echo $name
  namer=${name}.rxyz
#  file=$i
#  name="$(basename $i)" 
  if [ $rate -eq 1 ]; then
     dum_en=$(sqlite3 ${tsdirhl}/TSs/tshl.db "select energy,zpe from tshl where name='$name'" | sed 's@|@ @g')
  elif [ $rate -eq 0 ]; then
     dum_en=$(sqlite3 ${tsdirhl}/TSs/tshl.db "select energy,g from tshl where name='$name'" | sed 's@|@ @g')
  fi
  echo "$dum_min0"   >tmp_en
  echo "$dum_en"    >>tmp_en
  if [ $rate -eq 1 ]; then
     awk '{e[NR]=$1;ezpe[NR]=$2};END{printf "%10s %10.6f\n","'$namer'",(e[2]-e[1])*627.51+ezpe[2]-ezpe[1]}' tmp_en >> tmp_nonoren
  elif [ $rate -eq 0 ]; then
     awk '{e[NR]=$1;ezpe[NR]=$2};END{printf "%10s %10.6f\n","'$namer'",(e[2]-e[1]+ezpe[2]-ezpe[1])*627.51}' tmp_en >> tmp_nonoren
  fi
done
file=tmp_nonoren
echo "ordered energies"
awk '{a[++d]=$2}
END{
n = asort(a,b)
for (i=1; i<=n; i++){
 print b[i]
 }
}' $file > tmp_oren
echo "final energies"
paste tmp_oren tmp_nonoren | awk '{oe[NR]=$1;flag[NR]=$2;e[NR]=$3} 
END{
i=1
while(i<=NR){
  j=1
  while(j<=NR){
    if(oe[i]==e[j]) print "TS",i,flag[j],e[j]
    ++j
    }
  ++i
  }
}'  > $tsdirhl/TSs/SORTED/TSlist_sorted
for i in $(awk '{print $2}' $tsdirhl/TSs/SORTED/TSlist_sorted)
do
    name=$(awk 'NR=='$i',NR=='$i'{print $3}' $tsdirhl/TSs/SORTED/TSlist_sorted)
    namenr=$(basename $name .rxyz)
    echo $i $name
##keep this for the moment
    name0=TS$i
    names=TS${i}_${name}
#    cp $tsdirhl/TSs/${name} $tsdirhl/TSs/SORTED/${names}
##insert data into tsshl
    sqlite3 "" "attach '${tsdirhl}/TSs/tshl.db' as tshl; attach '${tsdirhl}/TSs/SORTED/tsshl.db' as tsshl; insert into tsshl (natom,name,lname,energy,zpe,g,geom,freq) select natom,'$name0','$names',energy,zpe,g,geom,freq from tshl where name='$namenr';"
done

# We now go on with the MINs
if [ -f tmp_oren ]; then rm tmp_oren ; fi
if [ -f tmp_nonoren ]; then rm tmp_nonoren ; fi

echo "non-ordered energies"
for name in $(sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "select name from minnrhl")
#for i in $(ls $tsdirhl/MINs/norep/min*.rxyz)
do
  echo $name
  namer=${name}.rxyz
  if [ $rate -eq 1 ]; then
     dum_en=$(sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "select energy,zpe from minnrhl where name='$name'" | sed 's@|@ @g')
  elif [ $rate -eq 0 ]; then
     dum_en=$(sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "select energy,g from minnrhl where name='$name'" | sed 's@|@ @g')
  fi
  echo "$dum_min0"   >tmp_en
  echo "$dum_en"    >>tmp_en
  if [ $rate -eq 1 ]; then
     awk '{e[NR]=$1;ezpe[NR]=$2};END{printf "%10s %10.6f\n","'$namer'",(e[2]-e[1])*627.51+ezpe[2]-ezpe[1]}' tmp_en >> tmp_nonoren
  elif [ $rate -eq 0 ]; then
     awk '{e[NR]=$1;ezpe[NR]=$2};END{printf "%10s %10.6f\n","'$namer'",(e[2]-e[1]+ezpe[2]-ezpe[1])*627.51}' tmp_en >> tmp_nonoren
  fi
done
file=tmp_nonoren
echo "ordered energies"
awk '{a[++d]=$2}
END{
n = asort(a,b)
for (i=1; i<=n; i++){
 print b[i]
 }
}' $file > tmp_oren
echo "final energies"
paste tmp_oren tmp_nonoren | awk '{oe[NR]=$1;flag[NR]=$2;e[NR]=$3} 
END{
i=1
while(i<=NR){
  j=1
  while(j<=NR){
    if(oe[i]==e[j]) print "MIN",i,flag[j],e[j]
    ++j
    }
  ++i
  }
}' > $tsdirhl/MINs/SORTED/MINlist_sorted

sed 's/_min/ min/g' $tsdirhl/MINs/SORTED/MINlist_sorted >tmp_mls
for i in $(awk '{print $2}' $tsdirhl/MINs/SORTED/MINlist_sorted)
do
    name=$(awk 'NR=='$i',NR=='$i'{print $3}' $tsdirhl/MINs/SORTED/MINlist_sorted)
    namenr=$(basename $name .rxyz)
    name2=$(awk 'NR=='$i',NR=='$i'{print $4}' tmp_mls)
    echo $i $name2
##keep this for the moment
    name0=MIN$i
    names=MIN${i}_${name2}
##insert data into mins
    sqlite3 "" "attach '${tsdirhl}/MINs/norep/minnrhl.db' as minnrhl; attach '${tsdirhl}/MINs/SORTED/minshl.db' as minshl;
    insert into minshl (natom,name,lname,energy,zpe,g,geom,freq) select natom,'$name0','$names',energy,zpe,g,geom,freq from minnrhl where name='$namenr';"
done

sed 's/_min/ min/g' $tsdirhl/MINs/SORTED/MINlist_sorted | awk '{print $1,$2,$4,$5}' > $tsdirhl/MINs/SORTED/MINlist_sorted.log
