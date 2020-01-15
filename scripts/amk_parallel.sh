#!/bin/bash
# default sbatch FT2
#SBATCH --output=amk_parallel-%j.log
#SBATCH --time=04:00:00
# partition selection

#_remove_this_in_ft_SBATCH -p shared --qos=shared
#SBATCH -c 1 --mem-per-cpu=2048
#SBATCH -n 8

# SBATCH --partition=cola-corta,thinnodes
# SBATCH -c 1
# SBATCH -n 24


#exe=$(basename $0)
# under batchs systems the scripts are copied to a generic script (in slurm slurm_script)
source utils.sh

#current working dir
cwd=$PWD
sharedir=${AMK}/share
exe="amk_parallel.sh"
# Printing the references of the method
print_ref

#Checking the input of this script and also the presence of some files
re='^[0-9]+$'
if [ $# -eq 2 ]; then
   if [[ ! $2 =~ $re ]]; then usage "Second argument must be a number";fi
   noj1=$(( $(sort -nr <(find . -maxdepth 1 -type d -print | grep 'batch' | sed 's@batch@@;s@./@@') | head -n1) +1 ))
   nojf=$(( $noj1 + $2 -1))
elif [ $# -eq 3 ]; then
   if [[ ! $3 =~ $re ]]; then usage "Third argument must be a number";fi
   noj1=$2
   nojf=$3
else
   usage "At least 2 arguments required"
fi
inputfile=$1
if [ "$inputfile" == "amk.dat" ]; then
   echo ""
   echo "READING amk.dat"
   echo ""
else
   echo ""
   echo "READING $inputfile"
   echo ""
   ln -sf $inputfile amk.dat
fi

#input file exists?
if [ ! -f $inputfile ]
then
   usage "The first argument of this script must be the inputfile"
fi

###reading input file
read_input
###
###checks. amk_parallel only with md>0
if [ $md -lt 0 ]; then
   echo "  =================================================================="
   echo "   $exe can only be employed for samplings involving MD:          "
   echo "                   MD, MD-micro and BXDE                          "
   echo "  =================================================================="
   exit 1
fi
###keywords check: molecule and rate
keywords_check
###
### ${molecule}.xyz and/or ${frA}.xyz ${frB}.xyz file check
### and create thdist and (eventually) thdist_vdw reference files
xyzfiles_check
##
sampling_calcs
##
amkscript=0
print_method_screening
##
echo ""
echo "CALCULATIONS START HERE"
#####for samplings 31 get the association complexes
if [ $sampling -eq 31 ]; then
   if [ ${xyz_exists} -eq 0 ]; then
      echo "Selecting a ${frA}-${frB} structure"
      exec_assoc
###screening repeated structures out
      echo "Screening the structures"
      screening_assoc.sh $inputfile
      rm -rf black_list* rotate.* tmp_*
   else
      echo "${frA}-${frB} structure detected in the working directory"
   fi
fi
###
cp ${molecule}.xyz ${molecule}_ref.xyz
####
kmcfilell=${tsdirll}/KMC/RXNet_long.cg_groupedprods
minfilell=${tsdirll}/MINs/SORTED/MINlist_sorted
if [ -f $kmcfilell ] && [ -f $minfilell ] && [ $mdc -ge 1 ]; then
   awk '{if($1!="temp") print $0}' $inputfile > tmp && mv tmp $inputfile  
fi
####
###Create tsdirll folder
if [ ! -d "$tsdirll" ]; then mkdir $tsdirll ; fi
###
sqlite3 ${tsdirll}/track.db "create table if not exists track (id INTEGER PRIMARY KEY,nts  INTEGER, noj1 INTEGER, nojf INTEGER, ntraj INTEGER, emin REAL, emax REAL, permin INTEGER, permax INTEGER);"

if [ $md -eq 0 ]; then
   et="temperatures"
   uet="K"
   flag="temp"
elif [ $md -eq 1 ]; then
   et="energies"
   uet="kcal/mol"
   flag="etraj"
else
   et="temperatures"
   uet="K"
   flag="temp"
fi
erange="$(awk '{if($1=="'$flag'") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile | sed 's/-/ /')"
nf="$(echo "$erange" | awk '{print NF}')"
echo ""
if [ $nf -eq 1 ]; then
   emin0="$(echo $erange | awk '{printf "%8.2f",$1}')"
   emax0="$(echo $erange | awk '{printf "%8.2f",$1}')"
elif [ $nf -eq 2 ]; then
   emin0="$(echo $erange | awk '{printf "%8.2f",$1}')"
   emax0="$(echo $erange | awk '{printf "%8.2f",$2}')"
elif [ $nf -eq 0 ]; then
   s=$(echo "3*$natom-6" | bc )
   emin0=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
   emax0=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
   if [ $md -eq 2 ]; then
      lstnm=$(awk 'BEGIN{lstnm=0};{if($1=="atoms" && NF==3) {lstnm=$3}};END{print lstnm}' $inputfile )
      thmass=$(awk 'BEGIN{thmass=0};{if($1=="thmass") {thmass=$2}};END{print thmass}' $inputfile )
      nlms=$(awk 'BEGIN{atoms=0};{if($1=="atoms" && $2!="all") atoms=$2;if($1=="atoms" && $2=="all") atoms=0};END{print atoms}' $inputfile )

      awk '{print $0}
      END{
      print "100"
      print '$nlms'
      if('$nlms'>0) print '$lstnm'
      print '$thmass'
      }' ${molecule}.xyz | termo.exe > /dev/null

      natefin=$(awk '/Number of atoms to be excited/{print $NF}' fort.66)
      rm -rf fort.66
      emin0=$(echo "335.51*$emin0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
      emax0=$(echo "335.51*$emax0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
   fi
else
   echo Check the value of $flag in $inputfile
   exit 1
fi
#set value of factor
id=$(sqlite3 ${tsdirll}/track.db "select max(id) from track" | awk '{print $1+1-1}')
if [ $id -gt 0 ]; then
   permin=$(sqlite3 ${tsdirll}/track.db "select permin from track where id='$id'")
   permax=$(sqlite3 ${tsdirll}/track.db "select permax from track where id='$id'")
   if [ $permin -gt 60 ]; then
      factormin=0.9
   elif [ $permin -ge 0 ] && [ $permin -le 60 ];then
      factormin=$(echo "10/9" | bc -l)
   else
      factormin=1
   fi
   if [ $permax -gt 60 ]; then
      factormax=$(echo "10/9" | bc -l)
   elif [ $permax -ge 0 ] && [ $permax -le 60 ];then
      factormax=0.9
   else
      factormax=1
   fi
else
   factormin=1
   factormax=1
fi
emin=$(echo "$emin0*$factormin" | bc -l | awk '{printf "%8.2f",$1}')
emax=$(echo "$emax0*$factormax" | bc -l | awk '{printf "%8.2f",$1}')
#faf is employed to avoid emin>emax situations
faf=$(echo "$emax-$emin" | bc -l | awk 'BEGIN{faf=1};{if($1<0) faf=0};END{print faf}')
if [ $faf -eq 0 ]; then
   s=$(echo "3*$natom-6" | bc )
   emin_sug=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
   emax_sug=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
   if [ $md -eq 2 ]; then
      emin_sug=$(echo "335.51*$emin_sug/$natom" | bc -l | awk '{printf "%8.2f",$1}')
      emax_sug=$(echo "335.51*$emax_sug/$natom" | bc -l | awk '{printf "%8.2f",$1}')
   fi
   echo You may consider changing the values of $flag in $inputfile
   echo Suggested range of ${et}: ${emin_sug}-${emax_sug}
   emin=$emin0
   emax=$emax0
fi
if [ $md -eq 0 ]; then
   sqlite3 ${tsdirll}/track.db "insert into track (noj1,nojf,emin,emax,permin,permax) values ($noj1,$nojf,$emin,$emax,-1,-1);"
else
   echo Range of ${et}: ${emin}-${emax} ${uet}
   sqlite3 ${tsdirll}/track.db "insert into track (noj1,nojf,emin,emax,permin,permax) values ($noj1,$nojf,$emin,$emax,-1,-1);"
   tmpinp="$(awk '{if($1=="sampling")
      {print $0
      print "'$flag' erange"}
   else if($1!="'$flag'") print $0}' $inputfile)"
   if [[ "$emin" == "$emax" ]]; then
      echo "$tmpinp" | sed 's/erange/'"$emin"'/' >$inputfile
   else
      echo "$tmpinp" | sed 's/erange/'"$emin"'-'"$emax"'/' >$inputfile
   fi
fi

echo ""
echo "A progress bar will pop up very shortly"
#The last job is tors.sh
nojf=$((nojf + 1))
#Then we submit the nojf-noj1+1 trajectory jobs
#ft2 slurm
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  if (( $nojf-$noj1+1 < $SLURM_NTASKS )); then 
    echo "WARNING: Number of trajectory jobs ($nojf-$noj1+1) lower than allocated tasks ($SLURM_NTASKS)."
  fi
fi
doparallel "runGP.sh {1} $molecule $cwd $nojf" "$(seq $noj1 $nojf)"

