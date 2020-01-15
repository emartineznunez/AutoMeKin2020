#!/bin/bash

# default sbatch FT2
#SBATCH --output=irc-%j.log
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

exe="irc.sh"
cwd=$PWD
sharedir=${AMK}/share
source utils.sh
#check the arguments of the script
if [ $# -gt 0 ]; then
   ci=$1
else
   ci="proceed"
fi

if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
if [ $ci != "screening" ] && [ $ci != "proceed" ]; then
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "                  Wrong argument                        "
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "To check what screening has done execute this script as:"
   echo "$exe screening"
   echo ""
   echo "To proceed with the irc execute this script as:"
   echo "$exe"
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   exit
fi
###Do screnning before anything else
screening.sh  $inputfile
if [ $ci == "screening" ]; then 
   echo ""
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "    Please check redundant and fragmented structures indicated in screening.log     "
   echo " If they are not what you expected you might change MAPEmax, BAPEmax and/or eigLmax "
   echo "Then, you can carry on with the IRC calculations, run this script without argument  "
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo ""
   exit
fi
###read input file
read_input
###
if [ ! -f ${tsdirll}/ts_mopac_failed ] && [ "$program_opt" != "mopac" ]; then
   echo "TSs not optimized with mopac" >  ${tsdirll}/ts_mopac_failed
fi
##
#remove tmp files
tmp_files=(tmp* bbfs.* *.arc *.mop coordir $tsdirll/*_mop.mop $tsdirll/*_mop.arc)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

if [ ! -d "$tsdirll/MINs" ]; then
   mkdir $tsdirll/MINs
fi
##create table for min
sqlite3 ${tsdirll}/MINs/min.db "create table if not exists min (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, sigma INTEGER,unique(name));"

# Optimize minref and calculate freqs
echo "$min_template"                       > ${molecule}_freq.mop
awk 'NF==4{print $0}' ${molecule}_ref.xyz >> ${molecule}_freq.mop
echo "$freq_template"                     >> ${molecule}_freq.mop
mopac ${molecule}_freq.mop 2>/dev/null

# First we copy min0 in MIN directory 
echo "Moving min0 to its final location"
if [ -f $tsdirll/MINs/min0.out ]; then
   echo "Calcs completed for min0"
else
   geom="$(get_geom_mopac.sh ${molecule}_freq.out | awk '{if(NF==4) print $0}')"
   sed 's/thermo/thermo('$temperature','$temperature')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template >  $tsdirll/MINs/min0.mop
   echo "$geom"  >> $tsdirll/MINs/min0.mop
   mopac $tsdirll/MINs/min0.mop 2>/dev/null
   e0=$(awk '/HEAT OF FORMATION =/{e=$5};END{print e}' $tsdirll/MINs/min0.out )
   zpe0=$(awk '/          ZERO POINT ENERGY/{zpe=$4};END{print zpe}' $tsdirll/MINs/min0.out )
   g_corr0=$(awk 'BEGIN{t='$temperature'}
      /          ZERO POINT ENERGY/{zpe=$4}
      /CALCULATED THERMODYNAMIC PROPERTIES/{ok=1}
   {if(ok==1 && $1 == '$temperature') {
   getline;getline;getline;getline;
   h=$3/1000;s=$5/1000;print zpe+h-t*s;exit}
   }' $tsdirll/MINs/min0.out )
   name=min0_0
   freq="$(get_freq_mopac.sh $tsdirll/MINs/min0.out)"
   sigma=$(awk '/SYMMETRY NUMBER/{print $NF;exit}' $tsdirll/MINs/min0.out)
   sqlite3 ${tsdirll}/MINs/min.db "insert into min (natom,name,energy,zpe,g,geom,freq,sigma) values ($natom,'$name',$e0,$zpe0,$g_corr0,'$geom','$freq',$sigma);"
fi

# Now we do things specific of IRC 
if [ ! -d "$tsdirll/IRC" ]; then mkdir $tsdirll/IRC ; fi
if [ ! -d "$tsdirll/TSs" ]; then mkdir $tsdirll/TSs ; fi
m=0
sqlite3 ${tsdirll}/inputs.db "drop table if exists mopac; create table mopac (id INTEGER PRIMARY KEY,name TEXT, unique(name));"
for name in $(awk '{print $3}' $tslistll)
do
  if [ -f $tsdirll/TSs/${name}_thermo.out ] && [ -f $tsdirll/IRC/${name}_ircf.out ] && [ -f $tsdirll/IRC/${name}_ircr.out ]; then
     calc1=$(awk 'BEGIN{calc=1};/MOPAC DONE/{calc=0};END{print calc}' $tsdirll/TSs/${name}_thermo.out)
     calc2=$(awk 'BEGIN{calc=1};/MOPAC DONE/{calc=0};END{print calc}' $tsdirll/IRC/${name}_ircf.out)
     calc3=$(awk 'BEGIN{calc=1};/MOPAC DONE/{calc=0};END{print calc}' $tsdirll/IRC/${name}_ircr.out)
     if [ $calc1 -eq 0 ] && [ $calc2 -eq 0 ] && [ $calc3 -eq 0 ]; then
        calc=0
     else
        calc=1
     fi
  else
     calc=1
  fi

  if [ $calc -eq 0 ]; then
    echo "Calcs completed for" $name
  else
     if [ "$program_opt" != "mopac" ]; then
        skip=$(awk 'BEGIN{skip=0};/'$name'/{skip=1};END{print skip}' ${tsdirll}/ts_mopac_failed)
        if [ $skip == 1 ]; then
           echo "TS $name has ben previously discarded-->(skip because mopac cannot optimize it)"
           continue
        else
           echo "$ts_template"                   > ${tsdirll}/${name}_mop.mop
           get_geom_g09.sh $tsdirll/${name}.out >> ${tsdirll}/${name}_mop.mop
           echo "$freq_template"                >> ${tsdirll}/${name}_mop.mop
           mopac ${tsdirll}/${name}_mop.mop 2>/dev/null
           fe="$(mopac_freq_ts.sh ${tsdirll}/${name}_mop.out 1)"
           fx="$(echo "$fe" | awk '{printf "%10.0f",$1}')"
           if [[ ("$fx" -gt "0") ]]; then
              get_geom_mopac.sh ${tsdirll}/${name}_mop.out | awk '{if(NF==4) print $0}' > tmp_geom
           else
              echo ${name}  >> ${tsdirll}/ts_mopac_failed
              rm ${tsdirll}/${name}_mop.out
              continue
           fi
        fi
     else
        get_geom_mopac.sh $tsdirll/${name}.out | awk '{if(NF==4) print $0}' > tmp_geom
     fi
     ((m=m+1))
     sed 's/thermo/thermo('$temperature','$temperature')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template > $tsdirll/TSs/${name}_thermo.mop
     cat tmp_geom >> $tsdirll/TSs/${name}_thermo.mop 
     sed 's/method/'"$method"' charge='$charge' irc= 1/g' $sharedir/freq_template1 > $tsdirll/IRC/${name}_ircf.mop
     sed 's/method/'"$method"' charge='$charge' irc=-1/g' $sharedir/freq_template1 > $tsdirll/IRC/${name}_ircr.mop
     cat tmp_geom >> $tsdirll/IRC/${name}_ircf.mop
     cat tmp_geom >> $tsdirll/IRC/${name}_ircr.mop
     sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2 | sed 's/force/cycles=5000 recalc=1/g'  >> $tsdirll/IRC/${name}_ircf.mop
     sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2 | sed 's/force/cycles=5000 recalc=1/g'  >> $tsdirll/IRC/${name}_ircr.mop
     echo -e "insert or ignore into mopac values (NULL,'$name');\n.quit" | sqlite3 ${tsdirll}/inputs.db
  fi
done
echo Performing a total of $m irc calculations
#Perform m parallel calculations
if [ $m -gt 0 ]; then
#ft2 slurm
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  if (( $m < $SLURM_NTASKS )); then 
    echo "WARNING: Number of irc calculations ($m) lower than allocated tasks ($SLURM_NTASKS)."
  fi
fi
   doparallel "runirc.sh {1} $tsdirll" "$(seq $m)"
fi

