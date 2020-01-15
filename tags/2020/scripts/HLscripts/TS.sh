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

#exe=$(basename $0)
exe="TS.sh"
cwd=$PWD
source utils.sh
#On exit remove tmp files
tmp_files=(tmp*)
#trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

#current working dir

# Printing the references of the method
print_ref
#Make sure the script is run with one argument
if [ $# -eq 0 ]; then 
   echo "One argument is required" 
   exit 1
else
   inputfile=$1
fi
#Make sure g09 is submitted to slurm in ft2
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  t=$(srun -N 1 -n 1 g09<${sharedir}/g09_test.dat | awk 'BEGIN{t=0};/Normal ter/{t=1};END{print t}')
else
  t=$(g09<${sharedir}/g09_test.dat | awk 'BEGIN{t=0};/Normal ter/{t=1};END{print t}')
fi

if [ $t -eq 0 ]; then
   echo "Please check that gaussian09 is installed in your computer and it can be invoked as g09"
   exit 1
fi

#Make sure the inputfile has not been deleted 
if [ ! -f $inputfile ]; then
   echo "The file $inputfile does not exist"
   exit
fi
if [ "$inputfile" == "amk.dat" ]; then
   echo "The name of the input file is amk.dat"
   echo ""
else
   echo "The name of the input file is $inputfile"
   echo "The file will be copied with the name amk.dat"
   echo ""
   ln -sf $inputfile amk.dat
fi
###Do screnning before anything else (just in case)
screening.sh  $inputfile
#####
##Reading High Level stuff
read_input
##checking that the HL stuff is read
if [ $noHLcalc -eq 0 ]; then echo Please, provide HighLevel keyword ; exit ; fi
min0=${molecule}
###Make $tsdirhl folder
if [ ! -d "$tsdirhl" ]; then
   echo "$tsdirhl does not exist. It will be created"
   mkdir $tsdirhl 2>tmp_err
   if [ -s tmp_err ]; then
      echo "check the path of tsdirll folder"
      exit
   fi
else
   echo "$tsdirhl already exists."
fi
###Make $tsdirhl/MINs folder
if [ ! -d "$tsdirhl/MINs" ]; then
   echo "$tsdirhl/MINs does not exist. It will be created"
   mkdir $tsdirhl/MINs
else
   echo "$tsdirhl/MINs already exists."
fi
##if using gaussian por low-level, path is complete
if [ "$program_opt" != "mopac" ]; then
   reduce=0
fi
if [ -z "$reduce" ]; then
   echo Keyword HL_rxn_network has not been specified
   exit 1
else
   if [ $reduce -lt 0 ]; then
      echo "Running the HL calculations for a subset of the TSs (bimolecular channels excluded)"
      reduce_RXNet.py $molecule 
   elif [ $reduce -gt 0 ]; then
      echo "Running the HL calculations for a subset of the TSs (bimolecular channels excluded)"
      reduce_RXNet.py $molecule $reduce
   else
      echo "Running the HL calculations for the whole reaction network"
   fi
fi
echo "Molecule name" $min0
echo "tsdirll is " $tsdirll

m=0
file=${tsdirll}/tslist
sqlite3 ${tsdirhl}/inputs.db "drop table if exists gaussian; create table gaussian (id INTEGER PRIMARY KEY,name TEXT, input TEXT, unique(name));"
for name in $(awk '{print $3}' $file)
do
  if [ -f $tsdirhl/${name}.log ]; then
     calc=$(awk 'BEGIN{calc=1;nt=0};/Normal termi/{++nt};/Error termi/{calc=0};END{if(nt=='$noHLcalc') calc=0;print calc}' $tsdirhl/${name}.log)
  else
     calc=1
  fi
  
  if [ $calc -eq 0 ]; then
    echo $tsdirhl/$name "already optimized"
  else
    ((m=m+1))
    echo $name "not optimized"
#construct g09 input file
    chk="$(echo %chk=$name)"
    cal="$(sed 's/tkmc/'$temperature'/;s@iop@'$iop'@;s@level1@'$level1'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/hl_input_template)" 
    if [ "$program_opt" != "mopac" ]; then
       geo="$(get_geom_g09.sh $tsdirll/$name.out | awk '{if(NF==4) print $0};END{print ""}')"
    else
       geo="$(get_geom_mopac.sh $tsdirll/$name.out | awk '{if(NF==4) print $0};END{print ""}')"
    fi

    ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo")"
    if [ $noHLcalc -eq 2 ]; then
       spc="$(sed 's/chk=/chk='$name'/;s@level2@'$level2'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/sp_template)"
       ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo"'\n\n'"$spc")"
    fi
    echo -e "insert or ignore into gaussian values (NULL,'$name','$ig09');\n.quit" | sqlite3 ${tsdirhl}/inputs.db
#    names[$m]=$name
  fi
done 
echo "$m TS opt calculations"

#now the initial minimum
if [ -f $tsdirhl/min0.log ]; then
   calc=$(awk 'BEGIN{calc=1;nt=0};/Normal termi/{++nt};/Error termi/{calc=0};END{if(nt=='$noHLcalc') calc=0;print calc}' $tsdirhl/min0.log)
else
   calc=1
fi

if [ $calc -eq 0 ]; then
   echo "min0 already optimized"
else
   ((m=m+1))
   echo "min0 not optimized"
#construct g09 input file
   chk="$(echo %chk=min0)"
   cal="$(sed 's/ts,noeigentest,//;s@iop@'$iop'@;s/charge/'$charge'/;s/mult/'$mult'/;s/tkmc/'$temperature'/;s@level1@'$level1'@' $sharedir/hl_input_template)"
   geo="$(awk '{if(NF==4) print $0};END{print ""}' ${min0}_ref.xyz)"
   ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo")"
   if [ $noHLcalc -eq 2 ]; then
      spc="$(sed 's/chk=/chk=min0/;s@level2@'$level2'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/sp_template)"
      ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo"'\n\n'"$spc")"
   fi
   echo -e "insert or ignore into gaussian values (NULL,'min0','$ig09');\n.quit" | sqlite3 ${tsdirhl}/inputs.db
#   names[$m]="min0"
fi
#Perform m parallel calculations
echo Performing a total of $m ts opt calculations
if [ $m -gt 0 ]; then
#   doparallel "runTS.sh $tsdirhl" "$(echo ${names[@]})"
#echo doparallel "runTS.sh {1} $tsdirhl" "$(seq $m)"
   doparallel "runTS.sh {1} $tsdirhl" "$(seq $m)"
fi

