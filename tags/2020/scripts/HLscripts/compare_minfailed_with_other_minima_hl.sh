#!/bin/bash
#this script serves to identify the failed minima (those that remain even after using MINFAILED.sh)
# it compares the structure of the failed minimum with those of the MIN already optimized
sharedir=${AMK}/share
source utils.sh
#remove tmp files
tmp_files=(tmp* deg* ConnMat labels mingeom ScalMat sprint.* listtss)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
inputfile=amk.dat
cwd=$PWD

#reading input file
read_input
###
awk '/failed/{print $3}' $tsdirhl/KMC/RXNet_long.cg > listtss
nfail=$(wc -l listtss | awk '{print $1}')
if [ $nfail -eq 0 ]; then
   echo "You are lucky. No problems with the structures. Go head"
   exit
else
   echo "For" $nfail "paths the minima are problematic" 
   echo "You should edit KMC/RXNet_long.cg according to the suggestions below"
fi

elements=${sharedir}/elements

echo "This file is just for the values of the computed minima" >tmp_values_ref

minpath=$tsdirhl/MINs/SORTED

#Do the calcs for the minima only once

totmin=`awk 'END{print $2}' $minpath/MINlist_sorted`
for file in $(sqlite3 ${minpath}/minshl.db "select lname from minshl")
do
   i=$(sqlite3 ${minpath}/minshl.db "select id from minshl where lname='$file'")
   echo $i "out of" $totmin "file=" $file
   
   sqlite3 ${minpath}/minshl.db "select natom,geom from minshl where lname='$file'"  | sed 's@|@\n\n@g' >mingeom
   createMat.py mingeom 1 $nA
##
##
   echo "1 $natom" | cat - ConnMat | sprint.exe |  awk '/Natom/{natom=$2}
   /Adjace/{i=1
   while(i<=natom){
     getline
     print $0
     i++
     }
   }' > tmp_adja


   awk '{if( NR == FNR) {l[NR]=$1;n[NR]=NR/10+1;tne=NR}}
   {if(FNR > 1 && NR >FNR ) {
      IGNORECASE = 1
      i=1
      while(i<=tne){
         if( $1 == l[i]) print n[i]
         i++
         }
     }
   }' $elements labels > tmp_atn

   cat tmp_atn tmp_adja >tmp_wrk

   awk '{if(NF==1) n[NR]=$1}
   {
   if(NF>1) {++i; for(j=1;j<=NF;j++) a[i,j]=$j}
   }
   END{
   ii=1
   while(ii<=i) {
     sumi=0
     j=1
     while(j<=i) {
        k=1
        while(k<=i) {
           sumi+=a[ii,j]*a[j,k]*n[ii]*n[j]*n[k]
           ++k
           }
        ++j
        }
     sum+=(sumi)^n[ii]
     ++ii
     }
   print '$i',sum}'  tmp_wrk >> tmp_values_ref


done


sed 's/.rxyz//g' minfailed_list >tmp01


for i in $(awk '{print $0} ' listtss)
do 
  ((n=n+1))
  name=$(awk '/'$i'/{print $0} ' minfailed_list)
  nmbr=$(awk '/'$i'/{print NR} ' minfailed_list)
  if [ "$nmbr" == "" ]; then
     echo "the structure that failed is not in the list of failed opts. Check the IRC"
     exit
  fi 
  name2=$(awk 'NR=='$nmbr',NR=='$nmbr'{print $0} ' tmp01)
  echo $natom > minfailed_structure
  echo '' >> minfailed_structure
  get_geom_g09.sh $tsdirhl/IRC/${name2}.log  >> minfailed_structure

# do the calc for the minfailed_structure

   createMat.py minfailed_structure 1 $nA
   echo "1 $natom" | cat - ConnMat | sprint.exe |  awk '/Natom/{natom=$2}
   /Adjace/{i=1
   while(i<=natom){
     getline
     print $0
     i++
     }
   }'  >tmp_adja

   awk '{if( NR == FNR) {l[NR]=$1;n[NR]=NR/10+1;tne=NR}}
   {if(FNR > 1 && NR >FNR ) {
      IGNORECASE = 1
      i=1
      while(i<=tne){
         if( $1 == l[i]) print n[i]
         i++
         }
     }
   }' $elements labels > tmp_atn


   cat tmp_atn tmp_adja >tmp_wrk
   awk '{if(NF==1) n[NR]=$1}
   {
   if(NF>1) {++i; for(j=1;j<=NF;j++) a[i,j]=$j}
   }
   END{
   ii=1
   while(ii<=i) {
     sumi=0
     j=1
     while(j<=i) {
        k=1
        while(k<=i) {
           sumi+=a[ii,j]*a[j,k]*n[ii]*n[j]*n[k]
           ++k
           }
        ++j
        }
     sum+=(sumi)^n[ii]
     ++ii
     }
   print sum}'  tmp_wrk > tmp_value_structure

cat tmp_value_structure tmp_values_ref | awk 'BEGIN{print "This structure is similar to MIN"} 
{if(NF==1) comp=$1
if(NF==2 && comp==$2) printf "%s ", $1
}
END{
print "\n"
}'

done

