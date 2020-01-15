#!/bin/bash
#default sbatch FT2 resources
#SBATCH --time=04:00:00
#SBATCH -n 4
#SBATCH --output=TS-%j.log
#_remove_this_in_ft_SBATCH --partition=cola-corta,thinnodes
#SBATCH --ntasks-per-node=2
#SBATCH -c 12
#
# SBATCH -p shared --qos=shared
# SBATCH --ntasks-per-node=2
# SBATCH -c 10


sharedir=${AMK}/share
source utils.sh

if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing"
   exit
fi

molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
charge=`awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' $inputfile `
mult=`awk 'BEGIN{mu=1};{if($1=="mult") mu=$2};END{print mu}' $inputfile `
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirHL_'$molecule'"}' $inputfile`

if [ ! -d "$tsdirhl/MINs/BO" ]; then
   echo "$tsdirhl/MINs/BO does not exist. It will be created"
   mkdir $tsdirhl/MINs/BO
else
   echo "$tsdirhl/MINs/BO already exists"
fi
##Reading HL stuff
read_input
##
sed 's/opt=(ts,noeigentest,calcall,noraman)/pop=nboread/;s@level1@'$level1'@;s/charge/'$charge'/;s/mult/'$mult'/;s@iop@'$iop'@' ${sharedir}/hl_input_template  > hl_input_nbo
sed -i '/temperature=tkmc/d' hl_input_nbo

kmcfilehl=$tsdirhl/KMC/RXNet_long.cg_groupedprods
minfilehl=$tsdirhl/MINs/SORTED/MINlist_sorted
confilehl=$tsdirhl/working/conf_isomer.out
rate=` awk '{if ($1=="Temperature") print "0";if($1=="Energy") print "1" }' $inputfile `
energy=` awk 'BEGIN{e=0};{if($1=="Energy") e=$2};END{print e}'  $inputfile `
en=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile `
factor=1.5

if [ -f $minfilehl ] && [ -f $kmcfilehl ]; then
   echo "Select the relevant minima to locate variational TSs"
   minn=`awk '/min0/{print $2}' $minfilehl`
   if [ ! -f $confilehl ]; then
      echo "File $confilehl does not exist and we cannot proceed"
      exit
   fi
   minok=`awk 'BEGIN{min='$minn'}
   {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
   j=1
   while(j<=iso[NR]){
      if('$minn'==m[NR,j]) min=m[NR,1]
      j++
      }
   }
   END{print min}' $confilehl `
fi

get_minn.sh $kmcfilehl $minok $en $factor

set `awk '{print $0}' minn`
for i
do 
if [ -f $tsdirhl/MINs/BO/$i"_bo.log" ]; then
   echo "NBO for" $i "already done"
else
  echo $i
  echo "%chk="$i"_bo" > $tsdirhl/MINs/BO/$i"_bo.dat"
  cat hl_input_nbo >> $tsdirhl/MINs/BO/$i"_bo.dat"
  awk '{if(NR>=3 && NF==4) print $0}' $tsdirhl/MINs/SORTED/MIN$i"_"*.rxyz >> $tsdirhl/MINs/BO/$i"_bo.dat"
  echo "" >> $tsdirhl/MINs/BO/$i"_bo.dat"
  echo "$"nbo bndidx "$"end  >> $tsdirhl/MINs/BO/$i"_bo.dat"
  echo "" >> $tsdirhl/MINs/BO/$i"_bo.dat"
  cp ${AMK}/share/sbatch_four_node_10h $tsdirhl/MINs/BO/job$i.sh
  echo "cd $PWD/$tsdirhl" >> $tsdirhl/MINs/BO/job$i.sh
  echo "g09<"$i"_bo.dat>"$i"_bo.log" >>  $tsdirhl/MINs/BO/job$i.sh
  sbatch $tsdirhl/MINs/BO/job$i.sh
fi
done
