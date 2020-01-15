#!/bin/bash
# default sbatch FT2
#SBATCH --output=llcalcs-%j.log
#SBATCH --time=04:00:00
# partition selection

#_remove_this_in_ft_SBATCH -p shared --qos=shared
#SBATCH -c 1 --mem-per-cpu=2048
#SBATCH -n 32

# SBATCH --partition=cola-corta,thinnodes
# SBATCH -c 1
# SBATCH -n 48

# first  arg is inputfile
# second arg is nbatches (200 is a good number)
# third  arg is niter  
#if no arguments are provided, then a gui pops up
exe="llcalcs.sh"
cwd="$PWD"
iter=0
sharedir=${AMK}/share
source utils.sh
# Printing the references of the method
print_ref
#
if [ $# -eq 0 ]; then
   FILE="$(zenity --file-selection --filename="$PWD/*.dat" --file-filter="*.dat" --title="Select the input file" 2> /dev/null)"
   inputfile="$(basename $FILE)"
   echo "Selected input file: $inputfile"
   answer="$(zenity --forms --title="llcalcs.sh GUI" --text="Add input data" \
      --add-entry="Number of tasks" \
      --add-entry="Number of iterations" \
      --add-entry="Max number of running tasks" 2>/dev/null | awk 'BEGIN{FS="|"};{print $1,$2,$3}' )"
   nbatch=$(echo "$answer" | awk '{print $1}')
   niter=$(echo "$answer" | awk '{print $2}')
   runningtasks=$(echo "$answer" | awk '{print $3}')
elif [ $# -eq 3 ]; then
   inputfile=$1
   nbatch=$2
   niter=$3
   if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
      runningtasks=$SLURM_NTASKS
   else
     echo "With three arguments it must be run under the SLURM batch system:"
     echo "sbatch $exe inputfile ntasks niter"
     exit 1
   fi
elif [ $# -eq 4 ]; then
   inputfile=$1
   nbatch=$2
   niter=$3
   runningtasks=$4
else
   echo You must provide zero or four arguments:
   echo "nohup $exe >llcalcs.log 2>&1 &"
   echo or
   echo "nohup $exe inputfile ntasks niter runningtasks >llcalcs.log 2>&1 &" 
   exit 1
fi
export runningtasks

###Are we in the right folder?
if [ ! -f $inputfile ];then
   echo "$inputfile is not in this folder"
   exit 1
fi
if [ -z $nbatch ] || [ -z $niter ]; then
   echo "Number of batches and/or number of iterations have not been set"
   exit 1
fi
#EMN. If nbatch=0 do not run dynamics.
if [ $nbatch -eq 0 ]; then niter=1 ; fi
#EMN
read_input
###
echo ""
echo "Number of iterations  = $niter"
echo "Tasks per iteration   = $nbatch"
echo ""
###checks and writing stuff
xyzfiles_check 
###
sampling_calcs
###
amkscript=0
print_method_screening
###
echo ""
echo "CALCULATIONS START HERE"
echo ""
iter=1
#set interactive mode to 0
inter=0
export inter
echo $$ > .script.pid
#
while [ $iter -le $niter ]; do
   export iter
   echo "*****************"
   echo "   Iter: ${iter}/${niter}"
   echo "*****************"
   echo "$iter/$niter" > iter.txt
   if [ $nbatch -gt 0 ]; then
      amk_parallel.sh $inputfile $nbatch >/dev/null
   fi
##check that tslistll file exists
   if [ ! -f $tslistll ]; then
      echo ""
      echo "ERROR:"
      echo "File ${tslistll} does not exist"
      echo "Check your batch directories"
      exit
   fi
   echo "    IRC calcs    "
   irc.sh > /dev/null
   echo "Optimizing minima"
   min.sh  > /dev/null
   echo " rxnetwork calcs "
   if [ $sampling -eq 31 ]; then
      rxn_network.sh allstates >/dev/null
   else
      rxn_network.sh >/dev/null
   fi
   echo "    KMC calcs    "
   kmc.sh > /dev/null
   echo "  End of iter $iter"
   echo ""  
   ((iter=iter+1))
done
echo "Making final folder: FINAL_LL_${molecule}"
echo ""
echo "TS structures/iteration:"
echo ""
track_view.sh
final.sh > /dev/null
echo ""
echo "END OF THE CALCULATIONS"
echo ""
