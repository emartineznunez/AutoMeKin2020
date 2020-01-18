#!/bin/bash

#Function for usage of amk_parallel
function usage {
   echo $*
   echo "Execute this script as in this example:"
   if [ "$exe" == "slurm_script" ]; then exe="sbatch amk_parallel.sh";fi
   echo ""
   echo " $exe FA.dat 100"
   echo ""
   echo " where FA.dat is the inputfile and 100 is the total number of tasks"
   echo ""
   exit 1
}

#Function for usage of amk
function usages {
   echo $*
   echo "Execute this script as in this example:"
   echo "  $exe amk.dat "
   echo "where amk.dat is the inputfile "
   exit 1
}

#Function for usage of llcalcs
function usagell {
   echo $*
   echo "Execute this script as in this example:"
   echo "  $exe amk.dat nbatches niter"
   echo "where amk.dat is the inputfile"
   echo "nbatches is the number of batches"
   echo "and niter the number of interactions"
   exit 1
}

#Function to printout the references
function print_ref {
   build="$(awk '{print $1}' ${AMK}/share/amk_build)"
   echo "****************************************"
   echo "                                        "
   echo "             AutoMeKin2020              "
   echo "                                        "
   echo "          version number ${build}       "
   echo "                                        "
   echo "****************************************"
}

#Function to submit jobs using slurm
function slurm {
#lets the user specify memory
#$SLURM_MEM_PER_CPU defined when option --mem-per-cpu= is used
#$SLURM_MEM_PER_NODE defined when option --mem= is used
   MEMPERCORE=$(sinfo -e -n $SLURM_NODELIST  -N -o "%m %c" -h | awk '{if(NR==1){min=$1/$2}else{new=$1/$2;if(new<min)min=new}}END{print min}')
   corespertask=${SLURM_CPUS_PER_TASK=1}
   if [ ! -z $SLURM_MEM_PER_NODE ]
   then
     #--ntasks-per-node= compulsory
     if [ ! -z $SLURM_NTASKS_PER_NODE ]
     then
       MEMPERCORE=$(( $SLURM_MEM_PER_NODE/ ($SLURM_NTASKS_PER_NODE * $corespertask) ))
     else
       echo "Please specify --ntasks-per-node= at the sbatch invocation"
       echo "or use --mem-per-cpu= option"
       exit 1
     fi
   fi
   if [ ! -z $SLURM_MEM_PER_CPU ]
   then
     MEMPERCORE=$SLURM_MEM_PER_CPU
   fi

#  SRUN="srun --exclusive -N1 -n1 --mem-per-cpu=$MEMPERCORE"
   SRUN="srun -N1 -n1 --mem=$(( $MEMPERCORE*$corespertask )) -c $corespertask --cpu_bind=none"
   runningtasks=$SLURM_NTASKS
}


#Function to submit jobs in parallel
function doparallel {
if [ ! -d amk_parallel-logs ];then mkdir amk_parallel-logs;fi

#slurm job?
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]
then
  slurm
else
# use as many concurrent tasks as number of cores-1 (if it is not defined or ill-defined)
  if [ -z $runningtasks ] || [ $runningtasks -gt $(nproc --ignore=1) ]; then
     runningtasks=$(nproc --ignore=1)
  fi
fi

# --delay 0.2 prevents overloading the controlling node
# -j is the number of tasks parallel runs
# --joblog makes parallel create a log of tasks that it has already run
# --resume-failed makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off

#progress bar only in interactive mode
if [ -z $inter ] && [ -z "$SRUN" ]; then
   parallel="parallel --bar --delay 0.2 -j $runningtasks --joblog amk_parallel-logs/${exe}-${iter}-tasks.log"
else
   parallel="parallel --delay 0.2 -j $runningtasks --joblog amk_parallel-logs/${exe}-${iter}-task.log"
fi
# this runs the parallel command we want
# in this case, we are running a script named runGP.sh
# parallel uses ::: to separate options. Here {0..99} is a shell expansion
# so parallel will run the command passing the numbers 0 through 99
# via argument {1}
COMMAND=$1
TASKS=$2
#$parallel "runGP.sh {1} $molecule" ::: $(seq $noj1 $nojf)
if [ -z "$SRUN" ]; then
   if [ -z $inter ]; then
      if [ -z $iter ]; then
         $parallel $COMMAND ::: $TASKS 2> >(zenity --progress --auto-close --no-cancel --title="parallel progress bar $exe" --width=500 --height=100 2> /dev/null) &
      else
         $parallel $COMMAND ::: $TASKS 2> >(zenity --progress --auto-close --no-cancel --title="parallel progress bar $exe iter=$iter" --width=500 --height=100 2> /dev/null)  &
      fi
   else
      nohup $parallel $COMMAND ::: $TASKS  >/dev/null 2>&1 & 
      echo $! > .parallel.pid
      wait
   fi
   echo $! > .parallel.pid
else
   $parallel $SRUN $COMMAND ::: $TASKS
fi
}

#function to printout the line in which the error occurred
function err_report {
    echo "Error on line $1 of $exe"
    rm -rf ${tmp_files[@]} 
    exit 1
}

function err_report2 {
    if [[ $exe == "amk.sh" ]]; then
       if [[ $1 != $2 ]]; then
          echo "Error on line $1 of $exe"
          rm -rf ${tmp_files[@]}
          exit 1
       fi
    else
          echo "Error on line $1 of $exe"
          rm -rf ${tmp_files[@]}
          exit 1
    fi
}

#function to cleanup on exit
function cleanup {
    rm -rf ${tmp_files[@]} 
    echo ""
    echo "Cleaning up tmp files and exiting $exe"
}

function cleanup2 {
    rm -rf ${tmp_files[@]} 
}


function read_input {
   mdc=0
   sharedir=${AMK}/share
   nb=$(basename $cwd | awk '{if(length($1)>8) print "wrkdir" ;else print $1}')
   srandseed=$(echo $nb | awk '/batch/{print $0}' | sed 's@batch@@' | awk '{print $1+100}')
   molecule=$(awk '{if($1=="molecule") print $2}'  $inputfile)
   if [ -f ${molecule}.xyz ]; then
      natom=$(awk 'NR==1{print $1}' ${molecule}.xyz)
   fi
   thd=$(awk 'BEGIN{f=0};{if($1=="NOcreatethdist") f=1};END{print f}' $inputfile)
#   fragments=$(awk 'BEGIN{f=0};{if($1=="fragments" && $2=="yes") f=1};{if($1=="fragments" && $2=="no") f=0};END{print f}' $inputfile)
   charge=$(awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' $inputfile)
   multiple_minima=$(awk 'BEGIN{mm=1};{if($1=="multiple_minima" && $2=="yes") mm=1};{if($1=="multiple_minima" && $2=="no") mm=0};END{print mm}' $inputfile)
   mult=$(awk 'BEGIN{mult=1};{if($1=="mult") mult=$2};END{print mult}' $inputfile)
# sampling:			md:	
# 0  ---> BXDE			0
# 1  ---> MD-micro		1
# 2  ---> MD			2
# 30 ---> association	       -1
# 31 ---> vdW			0
# 4  ---> external	       -1	
   sampling=$(awk 'BEGIN{sa=2};{if($1=="sampling" && $2=="BXDE") sa=0; if($1=="sampling" && $2=="MD-micro") sa=1;if($1=="sampling" && $2=="MD") sa=2;if($1=="sampling" && $2=="association") sa=30;if($1=="sampling" && $2=="vdW") sa=31; if($1=="sampling" && $2=="external") sa=4};END{print sa}'  $inputfile)
   md=$(awk 'BEGIN{md=2};{if($1=="sampling" && $2=="BXDE") md=0; if($1=="sampling" && $2=="MD-micro") md=1;if($1=="sampling" && $2=="MD") md=2;if($1=="sampling" && $2=="association") md=-1;if($1=="sampling" && $2=="vdW") md=0; if($1=="sampling" && $2=="external") md=-1};END{print md}'  $inputfile)
   rate=$(awk 'BEGIN{rate=-1};{if ($1=="Temperature") rate=0;if($1=="Energy") rate=1};END{print rate}' $inputfile)
   temperature=$(awk 'BEGIN{t=298};{if($1=="Temperature") t=$2};END{print t}' $inputfile)
   energy=$(awk 'BEGIN{e=0};{if($1=="Energy") e=$2};END{print e}'  $inputfile) 
   method=$(awk 'BEGIN{llcalc="pm7 threads=1"};{if($1=="LowLevel") {$1="";llcalc=$0" threads=1"}};END{print llcalc}' $inputfile)
   tsdirhl=$(awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirHL_'$molecule'"}' $inputfile)
   wrkmode=$(awk 'BEGIN{mode=1};{if($1=="post_proc" && $2=="bbfs" && NF==4) mode=$4};END{if(mode!=1) mode=0;print mode}' $inputfile)
##templates for mopac calcs
   min_template="$(cat $sharedir/freq_template1 | sed 's/method/'"$method"' charge='$charge'/g')"
   if [ $wrkmode -eq 0 ]; then 
      ts_template="$(cat $sharedir/ts_templateslow | sed 's/method/'"$method"' charge='$charge'/g')" 
      nppp=3
   elif [ $wrkmode -eq 1 ]; then 
      ts_template="$(cat $sharedir/ts_templatefast | sed 's/method/'"$method"' charge='$charge'/g')" 
      nppp=1
   fi
   freq_template="$(cat $sharedir/freq_template2 | sed 's/method/'"$method"' charge='$charge'/g')"
   bo_template="$(sed 's/method/'"$method"' charge='$charge' BONDS INT/g' $sharedir/freq_template1)"
##
   itrajn=$(awk 'BEGIN{tr=1};{if($1=="ntraj") tr=$2};END{print tr}'  $inputfile)
   nfs=$(awk 'BEGIN{time=500};{if($1=="fs") time=$2};END{print time}'  $inputfile)
   use_nfs=$(awk 'BEGIN{u=0};{if($1=="fs") u=1};END{print u}'  $inputfile)
   imag=$(awk 'BEGIN{imag=0};{if($1=="imagmin") imag=$2};END{print imag}'  $inputfile )
   tsdirll=$(awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$molecule'"}' $inputfile)
   postp_alg=$(awk 'BEGIN{p=1};{if($1=="post_proc" && $2=="bbfs") p=1};{if($1=="post_proc" && $2 =="bots") p=2};{if($1=="post_proc" && $2=="no") p=0};END{if('$sampling'==30) p=0;if('$sampling' ==1 || '$sampling' ==2) p=1; print p}' $inputfile)
   irange=$(awk 'BEGIN{irange=20};{if($1=="post_proc" && $2=="bbfs" && NF>=3) irange=$3};END{print irange}'  $inputfile)
   irangeo2=$(echo "scale=0; $irange/2" | bc )
   cutoff=$(awk 'BEGIN{co=200};{if($1=="post_proc" && $2=="bots" && NF>=3) co=$3};END{print co}' $inputfile)
   stdf=$(awk 'BEGIN{stdf=2.5};{if($1=="post_proc" && $2=="bots" && NF==4) stdf=$4};END{print stdf}' $inputfile)
   tslistll=${tsdirll}/tslist
   working=${tsdirll}/working
   workinghl=${tsdirhl}/working
   bu_ts=${tsdirll}/backup
   avgerr=$(awk 'BEGIN{a=0};{if($1=="MAPEmax") a=$2};END{print a}' $inputfile)
   bigerr=$(awk 'BEGIN{b=0};{if($1=="BAPEmax") b=$2};END{print b}' $inputfile)
   thdiss=$(awk 'BEGIN{t=0};{if($1=="eigLmax") t=$2};END{print t}' $inputfile)
   nmol=$(awk 'BEGIN{nmol=1000};{if($1=="nmol") nmol=$2};END{print nmol}'  $inputfile)
   imin=$(awk 'BEGIN{imin="min0"};{if($1=="imin") imin=$2};END{print imin}'  $inputfile)
   step=$(awk 'BEGIN{step=10};{if ($1=="Stepsize") step=$2};END{printf "%8.0f",step}'  $inputfile)
   emaxts=$(awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print 1.5*en}' $inputfile)
   frA=$(awk '{if($1=="fragmentA") {print $2;exit}}' $inputfile)
   frB=$(awk '{if($1=="fragmentB") {print $2;exit}}' $inputfile)
   if [ $sampling -ge 30 ];then
      if [ -f ${frA}.xyz ]; then
         nA=$(awk 'NR==1{print $1}' ${frA}.xyz)
      fi
      if [ -f ${frB}.xyz ]; then
         nB=$(awk 'NR==1{print $1}' ${frB}.xyz)
      fi
   else
      nA=$natom
      nB=0
   fi
   nassoc=$(awk 'BEGIN{n=100};{if($1=="Nassoc") n=$2};END{print n}' $inputfile)
   ptgr=$(awk 'BEGIN{impa=0.1};{if($1=="ImpPaths") impa=$2};END{print impa}' $inputfile)
   program_opt=$(awk 'BEGIN{popt="mopac"};{if($1=="LowLevel_TSopt") {popt=$2}};END{print tolower(popt)}' $inputfile)
   method_opt=$(awk 'BEGIN{llcalc="pm7"};{if($1=="LowLevel_TSopt") {llcalc=$3}};END{print tolower(llcalc)}' $inputfile)
   LLcalc=$(echo "$method_opt" | sed 's@/@ @g;s@u@@g' | awk 'BEGIN{IGNORECASE=1};{if($1=="hf") m="HF";else if($1=="mp2") m="MP2"; else if($1=="ccsd(t)") m="CCSDT";else m="DFT"};END{print m}' )
   atom1rot=$(awk 'BEGIN{ff=-1};{if($1=="rotate") ff=$2;if(ff=="com") ff=-1};END{print ff}' $inputfile)
   atom2rot=$(awk 'BEGIN{ff=-1};{if($1=="rotate") ff=$3;if(ff=="com") ff=-1};END{print ff}' $inputfile)
   dist=$( awk 'BEGIN{d=4.0};{if($1=="rotate") d=$4};END{print d}' $inputfile)
   distm=$(awk 'BEGIN{d=1.5};{if($1=="rotate") d=$5};END{print d}' $inputfile)
   factorflipv=$( awk '{if($1=="factorflipv") factor=$2};END{print factor}'  $inputfile )
   nbondsfrozen=$( awk 'BEGIN{nbf=0};{if($1=="nbondsfrozen") nbf=$2};END{print nbf}'  $inputfile ) 
   if [ $nbondsfrozen -gt 0 ]; then
      rm -rf fort.67 
      for i in $(seq 1 $nbondsfrozen)
      do
         bf="$(awk '{if($1=="nbondsfrozen") {i=1;while(i<='$i'){getline;if(i=='$i') print $1,$2;++i}}}' $inputfile)"
         echo "$bf" >> fort.67
      done
   fi
   nbondsbreak=$( awk 'BEGIN{nbb=0};{if($1=="nbondsbreak") nbb=$2};END{print nbb}'  $inputfile )
   if [ $nbondsbreak -gt 0 ]; then
      mdc=1
      awk '{if($1=="nbondsbreak") {nbb=$2;for(i=1;i<=nbb;i++) {getline;print $0}}}' $inputfile > fort.68
   fi
   nbondsform=$( awk 'BEGIN{nbfo=0};{if($1=="nbondsform") nbfo=$2};END{print nbfo}'  $inputfile )
   if [ $nbondsform -gt 0 ]; then
      mdc=1
      awk '{if($1=="nbondsform") {nbb=$2;for(i=1;i<=nbb;i++) {getline;print $0}}}' $inputfile > fort.69
   fi 
#HL stuff
   HLstring0="$(awk '{if($1=="HighLevel") print $2}' $inputfile)"
   HLstring="$(echo "$HLstring0" | sed 's@//@ @')"
   reduce=$(awk 'BEGIN{red=-1};{if($1=="HL_rxn_network") {if($2=="complete") red=0;if($2=="reduced" && NF==3) red=$3}};END{print red}' $inputfile)
   noHLcalc=$(echo $HLstring | awk 'BEGIN{nc=0};{nc=NF};END{print nc}')
   IRCpoints=$(awk 'BEGIN{np=100};{if($1=="IRCpoints") np=$2};END{print np}' $inputfile)
   iop=$(awk '{if($1=="iop") print $2}' $inputfile)
   level1=$(echo $HLstring | awk '{print $NF}')
   HLcalc1=$(echo "$level1" | sed 's@/@ @g;s@u@@g' | awk 'BEGIN{IGNORECASE=1};{if($1=="hf") m="HF";else if($1=="mp2") m="MP2"; else if($1=="ccsd(t)") m="CCSDT";else m="DFT"};END{print m}' )
#   echo High level calculations: "$HLstring0"
   if [ $noHLcalc -eq 1 ]; then
      HLcalc=$HLcalc1
   elif [ $noHLcalc -eq 2 ]; then
     level2=$(echo $HLstring | awk '{print $1}')
     HLcalc2=$(echo "$level2" | sed 's@/@ @g;s@u@@g' | awk 'BEGIN{IGNORECASE=1};{if($1=="hf") m="HF";else if($1=="mp2") m="MP2"; else if($1=="ccsd(t)") m="CCSDT";else m="DFT"};END{print m}' )
     HLcalc=$HLcalc2
   fi
###some few constants
   nfrag_th=0.005
}

##Function to run the association complexes 
function exec_assoc {
   assocdir=${cwd}/assoc_${frA}_${frB}
   if [ ! -d "$assocdir" ]; then mkdir $assocdir ; fi
###First, we need to select the best possible frA
#   select_AandB.sh $inputfile
###
   n="$(echo $nA $nB | awk '{print $1+$2}')"
   echo $n $nA $dist $distm > rotate.dat
   echo $atom1rot $atom2rot >> rotate.dat
   awk '{if(NF==4) print $0}' ${frA}.xyz >>rotate.dat
   awk '{if(NF==4) print $0}' ${frB}.xyz >>rotate.dat
   rm -rf ${assocdir}/structures
   for i in $(seq 1 $nassoc)
   do
      sed 's/method/'"$method"' charge='$charge' bonds/g' $sharedir/freq_template1 > ${assocdir}/assoc${i}.mop
      rotate.exe <rotate.dat>>${assocdir}/assoc${i}.mop
      echo $n >> ${assocdir}/structures
      echo '' >> ${assocdir}/structures
      rotate.exe <rotate.dat>>${assocdir}/structures
      sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2 >>${assocdir}/assoc${i}.mop
   done
   inter=0
   echo "Running $nassoc optimizations"
   doparallel "runAS.sh {1} $assocdir" "$(seq 1 $nassoc)"
}

function keywords_check {
###molecule keyword must be always present
   if [ -z $molecule ]; then
      echo " Your molecule keyword is missing in the input file"
      exit 1
   fi
##Kinetics section (except for association) and frA frB checks (association and vdw)
   if [ $sampling -ne 30 ]; then
      if [ $rate -eq -1 ]; then echo "Please provide a value for keywords Energy or Temperature in the Kinetics section"; exit; fi
      if [ $rate -eq 1 ] && [ $energy -eq 0  ] ; then
            echo "Please provide a value for the Energy"
            exit 1
      fi
   fi
   if [ $sampling -ge 30 ]; then
      if [ -z $frA ]; then
         echo keyword fragmentA is mandatory with association sampling
         exit 1
      fi
      if [ -z $frB ]; then
         echo keyword fragmentB is mandatory with association sampling
         exit 1 
      fi
   fi
###warning
}

function hl_print {
echo ""
echo "GENERAL    "
echo "Name of the system    =" $molecule
echo "Charge                =" $charge
echo "Multiplicity          =" $mult 
echo "High-level method     =" "$HLstring0"
echo ""
}


function xyzfiles_check  {
echo "GENERAL    "
echo "Name of the system    =" $molecule
echo "Charge                =" $charge
if [ $sampling -lt 30 ]; then
##check that xyz file is present
   if [ ! -f ${molecule}.xyz ]; then
      echo "${molecule}.xyz file does not exist"
      xyz_exists=0
      exit 1
   else
      xyz_exists=1
##remove second line if it exists
      awk 'NR==1{natom=$1;print natom"\n";getline
           for(i=1;i<=natom;i++) {getline; print $1,$2,$3,$4} }' ${molecule}.xyz > tmp && mv tmp ${molecule}.xyz 
##create reference distances : cov
      createthdist.py 0 $molecule $natom ${AMK}
   fi
else
   if [ ! -f ${frA}.xyz ]; then
      echo $frA".xyz does not exist"
      exit 1
   fi
   if [ ! -f ${frB}.xyz ]; then
      echo $frB".xyz does not exist"
      exit 1
   fi
   if [ -f ${molecule}.xyz ]; then
      xyz_exists=1
      xyzfile=${molecule}
      natom=$(awk 'NR==1{print $1}' ${molecule}.xyz)
   else
      xyz_exists=0
      cat ${frA}.xyz ${frB}.xyz | awk '{if(NF==4) print $0}' > tmp_ABe 
      natom=$(wc -l tmp_ABe | awk '{print $1}' )
      echo $natom > tmp_AB.xyz
      echo "" >> tmp_AB.xyz
      cat tmp_ABe >> tmp_AB.xyz
      xyzfile=tmp_AB
   fi
##check_vdw_atoms
   ok=$(check_vdw_atoms.py $xyzfile ${AMK})
   if [ $ok -eq 0 ];then
      echo Some of the atoms in your structure cannot be treated with this sampling
      exit
   fi
##create reference distances : cov and vdw
   createthdist.py 0 $xyzfile $natom ${AMK}
   mv thdist thdist_full
   createthdist.py 0 $xyzfile $nA ${AMK}
   mv thdist_full thdist
fi
echo "Number of atoms       =" $natom
met=$(echo $method | sed 's/threads=1//')
if [ $md -ge 0 ]; then
   echo "Low-level MD simul.   = mopac" $met
fi
if [ $sampling -ne 30 ]; then
   if [ "$program_opt" = "mopac" ]; then
      echo "Low-level TS optim.   =" $program_opt $met
   else
      echo "Low-level TS optim.   =" $program_opt "$method_opt"
   fi
else
   echo "Low-level A-B optim.  = mopac" $met
fi
}

function generate_dynamics_template {
deltat=1
dnc=2
ncycles=$(echo "scale=0; $dnc*$nfs" | bc)
##template for MD simulations
dynamics_template="$(cat $sharedir/dynamics_template)"
##This is stuff for biased dynamics*********************************
if [ ! -z $factorflipv ]; then
   echo "MD constraint: a factor of $factorflipv is employed to flip velocities"
   tmp0="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug vflip='$factorflipv'/')"
   dynamics_template="$tmp0"
fi
if [ $nbondsfrozen -gt 0 ]; then
   echo "MD constraint: $nbondsfrozen bonds are constrained using AXD"
   tmp0="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug nbondsfrozen='$nbondsfrozen'/')"
   dynamics_template="$tmp0"
fi
if [ $nbondsbreak -gt 0 ] && [ $mdc -eq 1 ]; then
   echo "MD constraint: a force is applied to $nbondsbreak bonds to promote their breakage"
   tmp1="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug nbondsbreak_ase='$nbondsbreak'/')"
   dynamics_template="$tmp1"
fi
if [ $nbondsform -gt 0 ] && [ $mdc -eq 1 ]; then
   echo "MD constraint: a force is applied to $nbondsform pairs of atoms to promote bond formation"
   tmp2="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug nbondsform_ase='$nbondsform'/')"
   dynamics_template="$tmp2"
fi
echo ""
dytem0="$(echo "$dynamics_template" | sed 's/method/'"$method"' charge='$charge'/g')"
dytem1="$(echo "$dytem0" | sed 's/ncycles/'$ncycles'/;s/deltat/'$deltat'/')"
}

function sampling_calcs {
if [ $sampling -eq 1 ]; then
###Energy can be a single value or a range of values
   erange="$(awk '{if($1=="etraj") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile | sed 's/-/ /')"
   nf="$(echo "$erange" | awk '{print NF}')"
   if [ $nf -eq 1 ]; then
      excite="$(echo $erange | awk '{print $1}')"
   elif [ $nf -eq 2 ]; then
      data="$(echo "$erange" | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')" 
      excite=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=$(echo "scale=0; $irange/2" | bc )
   elif [ $nf -eq 0 ]; then
      s=$(echo "3*$natom-6" | bc )
      emin0=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
      emax0=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
      data="$(echo $emin0 $emax0 | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')" 
      excite=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=$(echo "scale=0; $irange/2" | bc )
   else
      echo "Check the value of etraj"
      exit
   fi
###
   if [ -z "$excite" ]; then
      echo "Please provide an energy for the trajectories using keyword etraj"
      exit
   fi 
   lstnm=$(awk 'BEGIN{lstnm=0};{if($1=="modes" && NF==3) {lstnm=$3}};END{print lstnm}' $inputfile )
   nlms=$(awk 'BEGIN{modes=0};{if($1=="modes" && $2!="all") modes=$2;if($1=="modes" && $2=="all") modes=0};END{print modes}' $inputfile )
   seed=$(awk 'BEGIN{seed=0};/seed/{seed=$2};END{print seed}' $inputfile )
elif [ $sampling -eq 2 ]; then
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
   trange="$(awk '{if($1=="temp") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile | sed 's/-/ /')"
   nf="$(echo "$trange" | awk '{print NF}')"
   if [ $nf -eq 1 ]; then
      excite="$(echo $trange | awk '{print $1}')"
   elif [ $nf -eq 2 ]; then
      data="$(echo "$trange" | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')"
      excite=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=$(echo "scale=0; $irange/2" | bc )
   elif [ $nf -eq 0 ]; then
      s=$(echo "3*$natom-6" | bc )
      emin0=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
      emax0=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
      tmin0=$(echo "335.51*$emin0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
      tmax0=$(echo "335.51*$emax0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
      data="$(echo $tmin0 $tmax0 | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')"
      excite=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=$(echo "scale=0; $irange/2" | bc )
   else
      echo "Check the value of temp"
      exit
   fi
   if [ -z "$excite" ]; then
      echo "Please provide a temperature for the trajectories using keyword temp"
      exit
   fi 
fi
}

function print_method_screening {
echo ""
echo "METHOD     "
if [ $sampling -eq 0 ]; then
   echo "BXDE sampling" 
   if [ $postp_alg -eq 1 ];then
      printf "BBFS algorithm details:\nTime window (fs)      = $irange \nAttempts/single path  = $nppp \n" 
   elif [ $postp_alg -eq 2 ];then
      printf "BOTS algorithm details:\nCutoff freq (cm-1)    = $cutoff \nNumber of std devs.   = $stdf \n" 
   fi
elif [ $sampling -eq 1 ]; then
   echo "MD-micro sampling" 
   if [ $amkscript -eq 1 ]; then
      printf "BBFS algorithm details:\nTime window (fs)      = $irange \nAttempts/single path  = $nppp \n" 
   else
      printf "BBFS algorithm details:\nAttempts/single path  = $nppp \n" 
   fi
   eprint="$(awk '{if($1=="etraj") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile)"
   if [ -z $eprint ];then
      eprint="Automatically selected"
   fi
   if [ $nlms -gt 0 ]; then
      echo "# of modes excited    =" $nlms 
      echo "Modes excited:" $lstnm 
   else
      echo "All normal modes are excited "
   fi
   if [ $amkscript -eq 1 ]; then
      echo "Energy (kcal/mol)     =" $excite 
   else
      echo "Energy (kcal/mol)     =" $eprint
   fi
elif [ $sampling -eq 2 ]; then
   echo "MD sampling" 
   if [ $amkscript -eq 1 ]; then
      printf "BBFS algorithm details:\nTime window (fs)      = $irange \nAttempts/single path  = $nppp \n" 
   else
      printf "BBFS algorithm details:\nAttempts/single path  = $nppp \n" 
   fi
   tprint="$(awk '{if($1=="temp") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile)"
   if [ -z $tprint ];then
      tprint="Automatically selected"
   fi
   if [ $nlms -gt 0 ]; then
      echo "# of atoms excited    =" $nlms 
      echo "Atoms excited:" $lstnm 
      echo "Atoms with masses greater than" $thmass "will receive kinetic energy"
   else
      echo "All atoms are excited "
   fi
   if [ $amkscript -eq 1 ]; then
      echo "Temperature (K)       =" $excite 
   else
      echo "Temperature (K)       =" $tprint
   fi
elif [ $sampling -eq 30 ]; then
   echo "Association sampling"
   echo ""
elif [ $sampling -eq 31 ]; then
   echo "vdW sampling"
   if [ $postp_alg -eq 1 ];then
      printf "BBFS algorithm details:\nTime window (fs)      = $irange \nAttempts/single path  = $nppp \n" 
   elif [ $postp_alg -eq 2 ];then
      printf "BOTS algorithm details:\nCutoff freq (cm-1)    = $cutoff \nNumber of std devs.   = $stdf \n" 
   fi
elif [ $sampling -eq 4 ]; then
   echo "external sampling"
   if [ $postp_alg -eq 1 ];then
      printf "BBFS algorithm details:\nTime window (fs)      = $irange \nAttempts/single path  = $nppp \n" 
   elif [ $postp_alg -eq 2 ];then
      printf "BOTS algorithm details:\nCutoff freq (cm-1)    = $cutoff \nNumber of std devs.   = $stdf \n" 
   fi
   echo ""
else
   echo "No sampling provided. Please check your inputfile"
   echo ""
   exit
fi
if [ $md -ge 0 ]; then
   echo "Number of trajs       =" $itrajn
## BXDE with its own totaltime
   if [ $md -ne 0 ]; then
      echo "Total time (fs)       =" $nfs
   else
      if [ $use_nfs -eq 0 ]; then
         echo "Total time (fs)       = 5000" 
      else
         echo "Total time (fs)       =" $nfs
      fi
   fi
fi
echo ""
echo "SCREENING   "
if [ $postp_alg -ge 1 ]; then
   echo "Min imag freq (cm-1)  =" $imag
   echo "Max energy (kcal/mol) =" $emaxts
   if [ ! -d "$tsdirll" ]; then
      mkdir $tsdirll 2>tmp_err
      if [ -s tmp_err ]; then
         echo "check the path of tsdirll folder"
         exit
      else
         rm -rf tmp_err
      fi
   fi
fi
echo "Max value of MAPE     =" $avgerr
echo "Max value of BAPE     =" $bigerr
echo ""
}

function frag_check {
   createMat.py ${molecule}.xyz 1 $nA
###
   ndis=$(echo "1 $natom" | cat - ConnMat | sprint.exe | awk '/Results for the Laplacian/{getline; nconn=0
   for(i=6;i<=NF;i++) if($i<='${nfrag_th}') {++nconn} ; print nconn}' )
   if [ $ndis -gt 1 ]; then
      echo ""
      echo "************************************************************************"
      echo "Warning:                                                                "
      echo "Your input structure (file ${molecule}.xyz) contains $ndis fragments    "
      echo "The Kinetics results (is available) are meaningless                     "
      echo "************************************************************************"
   fi
}
