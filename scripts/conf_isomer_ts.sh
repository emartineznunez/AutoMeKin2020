#!/bin/bash
sharedir=${AMK}/share

source utils.sh
#remove tmp files
tmp_files=(ConnMat deg* mingeom ScalMat sprint.out wrk cits c_i_ts ${working}/*_data)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
cwd=$PWD

#reading input
read_input
###

if [ $1 -eq 0 ]; then
   tsdir=${tsdirll}
   database=tss
else
   tsdir=${tsdirhl}
   database=tsshl
fi
#redefine working if necessary
working=${tsdir}/working
echo "Screening" > $working/TSlist_screened

elements=${sharedir}/elements

# Loop over the TS structures
for name in $(sqlite3 ${tsdir}/TSs/SORTED/${database}.db "select lname from ${database}")
do 
   echo $name
   sqlite3 ${tsdir}/TSs/SORTED/${database}.db "select natom,geom from ${database} where lname='$name'" | sed 's@|@\n\n@g'  >mingeom

   createMat.py mingeom 1 $nA

   echo "1 $natom" | cat - ConnMat | sprint.exe >sprint.out

   awk '{if( NR == FNR) {l[NR]=$1;n[NR]=NR/10+1;tne=NR}}
   {if(FNR > 1 && NR >FNR ) {
      IGNORECASE = 1
      i=1
      while(i<=tne){
         if( $1 == l[i]) print n[i]
         i++
         }
     }
   }' $elements mingeom > wrk 

   awk '/Natom/{natom=$2}
   /Adjace/{i=1
   while(i<=natom){
     getline
     print $0
     i++
     }
   }' sprint.out >>wrk

   paste <(awk 'NF==4{print $1}' mingeom) <(deg.sh) >deg.out

   deg_form.sh > deg_form.out
   format.sh $name $working $thdiss
   echo $name "data" >>  $working/TSlist_screened
   cat $working/$name"_data" >> $working/TSlist_screened

   awk '{if(NF==1) n[NR]=$1}
   {if(NF>1) {++i; for(j=1;j<=NF;j++) a[i,j]=$j;a[i,i]=n[i]} }
   END{
   print i
   for(ii=1;ii<=i;ii++) {
      for(j=1;j<=i;j++) 
        printf("%2.1f ",a[ii,j])
        print "" 
      }
   }' wrk | diag.exe >>$working/TSlist_screened
done
#Loooking for conf. isomers
echo "Looking for conf. isomers"
reduce2.sh $working TS 

echo "Printing c_i_ts"
awk '{min[NR]=$1;
for(i=2;i<=NF;i++) {n[NR,i]=$i}
if(NR>1){
j=1
while(j<NR){
err=0
for(i=2;i<=NF;i++) {
    dif=n[NR,i]-n[j,i]
    err+=dif*dif
    }
  if(err==0) {
     ++ijk
     io[ijk]=min[NR]
     fl=1
     ii=1
     while(ii<=ijk){
       if(min[j]==io[ii]) fl=0
       ii++
       }
     if(fl==1) print min[j],min[NR],err
  }
  j++
  }
}
}' $working/TSlist_screened.red | awk 'BEGIN{ORS=" "} 
{
++n[int($1)];iso[int($1),n[int($1)]+1]=int($2)
}
END{
h=1
while(h<=1000){
  if(n[h]>0) {
    iso[h,1]=h
    for(i=1;i<=n[h]+1;i++) print iso[h,i] 
    print "\n"
    }
  h++
  }
}' | awk 'BEGIN{ORS=" "} 
{
if(NF==0) exit
for(i=1;i<=NF;i++) a[i]=$i
for(i=NF+1;i<=1000;i++) a[i]=-9999
n=asort(a)
for (i=1; i<=n; i++) if(a[i] != -9999) print a[i]
print "\n"
}' > c_i_ts
###If c_i_ts is empty exit here after creating an empty conf_isomer_ts.out file
###
file=c_i_ts
if [ -f $file ]; then
   nof=$(awk 'BEGIN{nf=0};NR==1{nf=NF};END{print nf}' $file)
else
   nof=0
fi
###

if [ $nof -ge 1 ]; then
   echo "Detecting conformational isomers"
else
   echo -n > $working/conf_isomer_ts.out 
   echo No conformational isomers of the ts
   exit 0 
fi
#Check that conf_isomer.out is not empty
###
cire=$working/conf_isomer.out
if [ -f $cire ]; then
   nof=$(awk 'BEGIN{nf=0};NR==1{nf=NF};END{print nf}' $cire)
else
   nof=0
fi
###
##Now split c_i_ts if $cire exists
if [ $nof -ge 1 ]; then
   echo "Now split c_i_ts"
else
   echo -n > $working/conf_isomer_ts.out 
   echo No conformational isomers of min nor ts 
   exit 0 
fi


rxnf=${tsdir}/KMC/RXNet
rm -f cits
for i in $(awk '{print NR}' c_i_ts)
do
  ni=$(awk 'NR=='$i',NR=='$i'{print NF}' c_i_ts)
  declare -a its=( $(for i in $(seq 1 $ni); do echo 0; done) )
  declare -a min01=( $(for i in $(seq 1 $ni); do echo 0; done) )
  declare -a min02=( $(for i in $(seq 1 $ni); do echo 0; done) )
  n=0
  ns=1
  for j in $(seq 1 $ni)
  do
     its[$j]=$(awk 'NR=='$i',NR=='$i'{print $'$j'}' c_i_ts)
     min1=$(awk 'BEGIN{min=0};{if($2=='${its[$j]}' && $7=="MIN") min=$8};END{print min}' $rxnf) 
     min2=$(awk 'BEGIN{min=0};{if($2=='${its[$j]}' && $10=="MIN") min=$11};END{print min}' $rxnf) 
     min01[$j]=$(awk 'BEGIN{min0=0};{for(i=1;i<=NF;i++) {if($i=='$min1') min0=$1 }};END{print min0}' $cire) 
     min02[$j]=$(awk 'BEGIN{min0=0};{for(i=1;i<=NF;i++) {if($i=='$min2') min0=$1 }};END{print min0}' $cire) 
     if [ ${min01[$j]} -eq 0 ]; then  
        unset its[$j] 
        unset min01[$j] 
        unset min02[$j] 
     fi
  done
  lia=${#its[@]}
  if [ $lia -gt 1 ]; then
     echo "${its[@]:1:$ni}"   >> cits
     echo "${min01[@]:1:$ni}" >> cits
     echo "${min02[@]:1:$ni}" >> cits
     echo "analyze" >> cits
  fi
done
#Check that cits exists
###
file=cits
if [ -f $file ]; then
   nof=$(awk 'BEGIN{nf=0};NR==1{nf=NF};END{print nf}' $file)
else
   nof=0
fi
###
if [ $nof -ge 1 ]; then
   echo "Printing conformational isomers"
else
   echo -n > $working/conf_isomer_ts.out 
   echo No conformational isomers of the ts
   exit 0 
fi
#
#initialize na
na=1
i=0
while [ $na -gt 0 ]; do
  ((i=i+1))
  echo "Iter # $i" 
###############
  output="$(awk 'BEGIN{na=0;tot=0}
  {if($1=="analyze") {
    na=0
    ok0=0
    if(tot>0) {
      for(i=1;i<=la;i++) {
         ok=1
         if(arr[2,1]!=arr[2,i] || arr[3,1]!=arr[3,i]) ok=0
         if(ok==1) printf "%s ",arr[1,i]
         if(ok==0 ) {++ok0;n1[ok0]=arr[1,i];n2[ok0]=arr[2,i];n3[ok0]=arr[3,i]}
         }
         print "","++"
         if(ok==0) {
         for(j=1;j<=ok0;j++) printf "%s ",n1[j]
         print ""
         for(j=1;j<=ok0;j++) printf "%s ",n2[j]
         print ""
         for(j=1;j<=ok0;j++) printf "%s ",n3[j]
         print ""
         print "analyze"}
      }
    }
  }
  {if($NF=="++") print $0}
  {if(NF>0 && $1!="analyze" && $NF!="++") {++tot;++na;la=NF
  for(i=1;i<=NF;i++) arr[na,i]=$i
   }
  }' cits)"
  echo "$output" > cits 
###############
  na=$(grep analyze cits | wc -l)
done
awk '{for(i=1;i<=(NF-1);i++) {if(NF>2) printf "%s ",$i};{if(NF>2) print ""}}' cits > $working/conf_isomer_ts.out
echo "End of the calc to determine the TS conf. isomers"

