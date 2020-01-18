#!/bin/bash
#default sbatch FT2 resources
#SBATCH --time=12:00:00
#SBATCH -n 4
#SBATCH --output=hlcalcs-%j.log
#_remove_this_in_ft_SBATCH --partition=cola-corta,thinnodes
#SBATCH --ntasks-per-node=2
#SBATCH -c 12
#
# SBATCH -p shared --qos=shared
# SBATCH --ntasks-per-node=2
# SBATCH -c 10


exe="hlcalcs.sh"
cwd="$PWD"
sharedir=${AMK}/share
source utils.sh
###
print_ref
##
#if no arguments are provided, then a gui pops up 
if [ $# -eq 0 ]; then
   FILE="$(zenity --file-selection --filename="$PWD/*.dat" --file-filter="*.dat" --title="Select the input file" 2> /dev/null)"
   inputfile="$(basename $FILE)"
   echo "Selected input file: $inputfile"
   runningtasks="$(zenity --forms --title="hlcalcs.sh GUI" --text="Add input data" \
      --add-entry="Max number of running tasks" 2>/dev/null  )"
elif [ $# -eq 1 ]; then
   inputfile=$1
   if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
      runningtasks=$SLURM_NTASKS
   else
     echo "With one argument it must be run under the SLURM batch system:"
     echo "sbatch $exe inputfile"
     exit 1
   fi
elif [ $# -eq 2 ]; then
   inputfile=$1
   runningtasks=$2
else
   echo You must provide zero or two arguments:
   echo "nohup $exe >hlcalcs.log 2>&1 &"
   echo or
   echo "nohup $exe inputfile runningtasks >hlcalcs.log 2>&1 &"
   exit 1
fi
export runningtasks

###Are we in the right folder?
if [ ! -f $inputfile ];then
   echo "$inputfile is not in this folder"
   exit 1
fi
###
read_input
###
##Printing some stuff
hl_print
##start of calcs
echo ""
echo "CALCULATIONS START HERE"
echo ""
#set interactive mode to 0
inter=0
export inter
echo $$ > .script.pid
#
system="$(basename $inputfile .dat)"
echo "   Optimizing TSs    "
TS.sh $inputfile > /dev/null
echo "      IRC calcs      " 
IRC.sh  >/dev/null 
echo "  Optimizing minima  "
MIN.sh  >/dev/null 
echo "   rxnetwork calcs   "
if [ $sampling -eq 31 ]; then
   RXN_NETWORK.sh allstates >/dev/null
else
   RXN_NETWORK.sh >/dev/null
fi
echo "      KMC calcs      "
KMC.sh >/dev/null
echo "  Optimizing frags   "
PRODs.sh >/dev/null
echo ""
echo "Making final folder: FINAL_HL_${molecule}"
FINAL.sh >/dev/null 
echo ""
echo "END OF THE CALCULATIONS"
echo ""

