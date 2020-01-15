#!/bin/bash
#default sbatch FT2 resources
#SBATCH --time=08:00:00
#SBATCH -n 4
#SBATCH --output=PRODs-%j.log
#_remove_this_in_ft_SBATCH --partition=cola-corta,thinnodes
#SBATCH --ntasks-per-node=2
#SBATCH -c 12
#
# SBATCH -p shared --qos=shared
# SBATCH --ntasks-per-node=2
# SBATCH -c 10

exe="PRODs.sh"
sharedir=${AMK}/share
source utils.sh
#remove tmp files
tmp_files=(tmp* ff*)
trap cleanup EXIT INT

#current working dir
cwd=$PWD

if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
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

#####
read_input
##

##Initialize here the total number of frags
tnf=0
##Creating dir to make ab initio calcs
dir=$tsdirhl/PRODs/CALC
##Creating working dir to compare frags
working=$tsdirhl/PRODs/CALC/working
if [ ! -d "$dir" ]; then
   echo "$dir does not exist. It will be created"
   mkdir $dir
   mkdir $working
   echo "Screening" > $working/fraglist_screened
else
   echo "$dir already exists. Gathering info to working dir to avoid repeating structures"
   rm -rf $working
   mkdir $working
   echo "Screening" > $working/fraglist_screened
   for file in $(ls $dir/*-*.log)
   do
      ((tnf=tnf+1))
      get_geom_g09.sh $file > tmp_geom
      nn=$(basename $file .log)
      compare_frags.sh tmp_geom frag${tnf}_$nn $working
   done
fi
sqlite3 $dir/inputs.db "drop table if exists gaussian; create table gaussian (id INTEGER PRIMARY KEY,name TEXT, input TEXT);"

number=0
m=0
echo "PR list with frags" > $tsdirhl/PRODs/PRlist_frag
for name in $(sqlite3 ${tsdirhl}/PRODs/prodhl.db "select name from prodhl")
do
   ((number=number+1))
   line[$number]=$(awk '{if($2=="'$number'") print $0}' $tsdirhl/PRODs/PRlist)
   name0=$(echo $name | sed 's@_min@ min@g' | awk '{print $2}')
   formula="$(sqlite3 $tsdirhl/PRODs/prodhl.db "select natom,geom from prodhl where name='$name'" | sed 's@|@\n\n@g' | FormulaPROD.sh)"
   nfrag=$(awk '{print $1}' tmp_nf)
   echo "Number: $number Name: $name # of frags: $nfrag"
   rm -f ff$number
   for j in $(seq 1 $nfrag)
   do 
      charge=$(awk '{if(NR == FNR) {n[NR]=$1;++naf}}
      /Mulliken charges:/{getline
      i=1
      while(i<='$natom'){
        getline
        j=1
        while(j<=naf){
          if($1==n[j]) {q+=$3}
          j++
          }
        i++
        }
      }
      END{print int(q+0.5)}'  tmp_Frag$j $tsdirhl/IRC/${name0}.log)
      noe=$(noe.sh tmp_frag$j.xyz $charge)
      if [ $noe -eq 1 ]; then
         mult=1
      else
         mult=2
      fi
      name="$(awk '{if(NR==2) print $1}' tmp_frag$j.log)"

      sqlnamep=${name}.q${charge}.m${mult}
      nisql="$(sqlite3 $dir/inputs.db "select name from gaussian where name like '%$sqlnamep%'")"
      if [ -f ${dir}/${sqlnamep}-1.log ] && [ -z "$nisql" ]; then
         ni=$(ls ${dir}/${sqlnamep}-*.log | wc -l | awk '{print $1+1}')
      else
         ni="$(echo "$nisql" | awk 'BEGIN{FS="-"};{name=$1;n=$2};END{print n+1}')"
      fi

      nn=${name}.q${charge}.m${mult}-$ni
      echo $nn >> ff$number
      ((tnf=tnf+1))
      awk '{if(NF==4) print $0}' tmp_frag$j.xyz >tmp_geom
      compare_frags.sh tmp_geom frag${tnf}_$nn $working
      nl=$(awk 'END{print NF}' $working/fraglist)
     ((m=m+1))
##calc only if the frag is not repeated
      if [ $nl -eq 2 ]; then
         nnc=${name}_q${charge}_m${mult}_$ni
         chk="$(echo "%chk=$nnc")"
         cal="$(sed 's/ts,noeigentest,//;s/tkmc/'$temperature'/;s@level1@'$level1'@;s/charge/'$charge'/;s/mult/'$mult'/;s@iop@'$iop'@'  $sharedir/hl_input_template)"
         geo="$(awk '{if(NF==4) print $0};END{print ""}' tmp_frag$j.xyz)"
         ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo")"
         if [ $noHLcalc -eq 2 ]; then
            spc="$(sed 's/chk=/chk='$nnc'/;s@level2@'$level2'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/sp_template)"
            ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo"'\n\n'"$spc")"
         fi
      else
         ig09="$(echo salir)"
      fi
      echo -e "insert into gaussian values (NULL,'$nn','$ig09');\n.quit" | sqlite3 ${dir}/inputs.db
   done
####
#Make PRlist_frag file
   dumm=$(awk '{f[NR]=$1}
   END{
   i=1
   while(i<=NR){
     if(i<NR) printf "%s + ",f[i]
     if(i==NR) printf "%s",f[i]
     i++
     }
   print ""
   }' ff$number)
   echo "${line[$number]}" "$dumm" >> $tsdirhl/PRODs/PRlist_frag
####
done
#submit the calcs
echo Performing a total of $m opt calculations
if [ $m -gt 0 ]; then
   doparallel "runTS.sh {1} $dir" "$(seq $m)"
fi


