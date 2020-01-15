#!/bin/bash
source utils.sh

if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
exe="simplifyRXN.sh"
cwd=$PWD

#reading input
read_input
###

# ptgr is the percent of the total number of processes to be considered a relevant path
if [ $1 -eq 0 ]; then
   tsdir=${tsdirll}
   database=tss
else
   tsdir=${tsdirhl}
   database=tsshl
fi
lastmin=$(awk '{lm=$2};END{print lm}' $tsdir/MINs/SORTED/MINlist_sorted )
#
post=$( awk '{if ($1=="Temperature") print "T";if($1=="Energy") print "E" }' $inputfile )
value=$(awk '{if($1=="Energy") print $2};{if($1=="Temperature") print $2}' $inputfile)
if [ ! -f $tsdir/KMC/kmc$post$value.out ] || [ ! -f $tsdir/KMC/processinformation ]; then
   echo "You should run rxn_network2.sh first"
   exit
fi
suma=$(awk 'BEGIN{nn=1e20};/counts per process/{nn=NR};{if(NR>nn) suma+=$2};END{print suma}' $tsdir/KMC/kmc$post$value.out )
if [ $suma -eq 0 ]; then
   echo "Not even one process. Please, increase the simulation time"
   exit	
fi
awk 'BEGIN{nn=1e20};/counts per process/{nn=NR};{if(NR>nn)print $0}' $tsdir/KMC/kmc$post$value.out >tmp
cat $tsdir/KMC/processinformation >>tmp

echo "Threshold to get rid of paths (%)" $ptgr
awk '{if(NF==2) {max=$1
  npc[$1]=$2
  }
}
{if(NF==6) ts[$2]=$4}
END{i=1
while(i<=max){
  tot+=npc[i]
  i++
  }
i=1
while(i<=max){
  p=npc[i]/tot*100
  if(p>='$ptgr') print ts[i]
  i++
  }
}' tmp >tmp_grothzts

cat tmp_grothzts $tsdir/KMC/RXNet_long.cg_groupedprods > tmp2
awk '{if(NF==1){++igr;gr[igr]=$1}}
{if($1=="KMC") print $0}
{if($1=="number") print $0}
{
p=0
i=1
while(i<=igr){
   if($2==gr[i]) {p=1;break}
   i++
   }
if($1=="TS" && p==1) print $0}' tmp2  > $tsdir/KMC/RXNet.relevant

