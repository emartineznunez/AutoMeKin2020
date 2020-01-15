#!/bin/bash
source utils.sh
#On exit remove tmp files
tmp_files=(minn* mingeom*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
sharedir=${AMK}/share
inputfile=$1
cwd=$PWD

elements=${sharedir}/elements
###reding input
read_input
###

en=$(awk 'BEGIN{en=100;if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile )

cp ${frA}.xyz ${frA}.xyz_backup
rm -f minn minn.log
dirhl=$PWD"/tsdirHL_"$frA
if [ ! -d $dirhl ]; then
   echo $dirhl "does not exist and, therefore, the structure" ${frA}.xyz "will be employed for that fragment"
   rm ${frA}.xyz_backup
   exit
fi
#Getting other possibilities from tsdirHL
kmcfile=$PWD"/tsdirHL_"$frA"/KMC/RXNet_long.cg_groupedprods"
minfile=$PWD"/tsdirHL_"$frA"/MINs/SORTED/MINlist_sorted"
mindir=$PWD"/tsdirHL_"$frA"/MINs/SORTED"
confile=$PWD"/tsdirHL_"$frA"/working/conf_isomer.out"
minn=$(awk '/min0/{print $2}' $minfile)
minok=$(awk 'BEGIN{min='$minn'}
{for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
j=1
while(j<=iso[NR]){
   if('$minn'==m[NR,j]) min=m[NR,1]
   j++
   }
}
END{print min}' $confile )
factor=1
get_minn.sh $kmcfile $minok $en $factor

cat minn $confile | awk 'BEGIN{n=0} 
{ 
   ++n
   if(NF==1) lem[n]=$1
   if(NF>1) {++m;linea[m]=$0;linea1[m]=$1}
   }
END{
  j=1
  while(j<=n){
     p=1
     k=1
     while(k<=m){
        if(linea1[k]==lem[j]) {print linea[k];p=0}
        k++
        }
     if(p==1) print lem[j]
     j++
     }
}' >minn.log
if [ -d test_assoc ]; then
   rm -f test_assoc/*
else
   mkdir test_assoc
fi
n=0
awk '{for(i=1;i<=NF;++i) print $i,NR}' minn.log >minn.out
for i in $(awk '{for(i=1;i<=NF;++i) if(NR<=3) print $i}' minn.log)
do
   names="MIN"$i
   ((n=n+1))
   if [ $n -eq 1 ]; then
      sqlite3 $mindir/minshl.db "select geom from minshl where name='$names'" > mingeom0
      met_label=$(awk '{if( NR == FNR) {l[NR]=$1;tne=NR}}
      {if(FNR > 1 && NR >FNR ) {
         IGNORECASE = 1
         i=1
         while(i<=tne){
            if( $1 == l[i] && i>=21 && i<=30) print FNR
            if( $1 == l[i] && i>=39 && i<=48) print FNR
            if( $1 == l[i] && i>=72 && i<=80) print FNR
            i++
            }
        }
      }' $elements mingeom0)
   fi
   sed 's/method/'"$method"' charge='$charge' bonds/g' ${sharedir}/freq_template1 > test_assoc/${i}.mop
   cat mingeom0 >>test_assoc/${i}.mop
   sed 's/method/'"$method"' charge='$charge'/g' ${sharedir}/freq_template2 >> test_assoc/${i}.mop
   mopac test_assoc/${i}.mop 2>/dev/null
   val=$(awk 'BEGIN{huge=10^10}
   /BOND ORDERS AND VALENCIES/{getline
   getline
   getline
   getline
   i=1
   while(i<=huge){
     getline
     if($2=='$met_label') val=$NF
     if($2=="Cite") {print val;exit}
     i++
     }
   }' test_assoc/$i.out)
   iso=$( awk 'NR=='$n',NR=='$n'{print $2}' minn.out )
   echo $i $val $iso 
   echo $i $val $iso >> test_assoc/min_selected
done
smin=$(awk 'BEGIN{max=-10^10}
{p=10^3*2^(-$2)
if(p>max && $3<=3) {max=p;n=$1}
}
END{
print n
}' test_assoc/min_selected)
echo "For "$frA "we select minimum=" $smin "of dir "$dirhl
echo "That minimum will now be called= "${frA}.xyz
get_geom_mopac.sh test_assoc/$smin".out" > ${frA}.xyz  


