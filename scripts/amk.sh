#!/bin/bash
source utils.sh
#On exit remove tmp files
tmp_files=(ConnMat tmp* ScalMat *.arc *.mop fort.* partial_opt ts_opt *_backup rotate.dat minn black_list*)
trap 'err_report2 $LINENO $gauss_line' ERR
trap cleanup EXIT INT

##Defining paths and names
cwd=$PWD
sharedir=${AMK}/share
exe=$(basename $0)

# Printing the references of the method
print_ref
##
##Define input file and create symbolic link-->amk.dat
if [ $# -eq 0 ]; then usages "One argument is required" ; fi
inputfile=$1
if [ ! -f $inputfile ]; then
   echo "The file $inputfile does not exist"
   exit
fi
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

###Reading stuff from inputfile
read_input
###keywords check: molecule and rate
keywords_check
### check files and write some stuff 
xyzfiles_check
###Peforming some calcs for the various samplings
sampling_calcs
### print method and screening sections
amkscript=1
print_method_screening
#################################
##  Starting the calculations
#################################
echo ""
echo "CALCULATIONS START HERE"
#####for association and vdw get the association complexes
if [ $sampling -ge 30 ]; then
   if [ ${xyz_exists} -eq 0 ]; then
      echo "Selecting a ${frA}-${frB} structure"
      exec_assoc
      echo "Screening the structures"
      screening_assoc.sh $inputfile
      if [ ! -f ${molecule}.xyz ]; then
         echo "A structure for the ${frA}-${frB} complex could not be found"
         exit 1
      fi
   else
      echo "${frA}-${frB} structure detected in the working directory"
   fi
fi
###for association stop here
if [ $sampling -eq 30 ]; then
   echo ""
   echo "END OF THE CALCULATIONS"
   echo "Check your ${frA}-${frB} structure in file ${molecule}.xyz" 
   echo ""
   exit 
fi
##For MD-based methods continue 
##select starting structure
sel_mol.sh $inputfile $multiple_minima
##
frag_check
##lift MD-constraint in subsequent iterations 
kmcfilell=${tsdirll}/KMC/RXNet_long.cg_groupedprods
minfilell=${tsdirll}/MINs/SORTED/MINlist_sorted
if [ -f $kmcfilell ] && [ -f $minfilell ] && [ $mdc -ge 1 ] && [ $ndis -eq 1 ]; then mdc=0 ; fi
##template for the dynamics
generate_dynamics_template
##make temporary folders
rm -rf partial_opt  && mkdir partial_opt
rm -rf ts_opt   && mkdir ts_opt
if [ $sampling -ne 4 ]; then
   rm -rf coordir && mkdir coordir
fi
###Opt the starting structure
echo "Optimizing the starting structure"
echo "$min_template"                    > ${molecule}_freq.mop
awk  'NF==4{print $0}' ${molecule}.xyz >> ${molecule}_freq.mop
if [ $md -eq 1 ]; then 
echo "In MD-micro, a frequency calculation is performed as well"
   echo "$freq_template" >> ${molecule}_freq.mop
fi
mopac ${molecule}_freq.mop 2>/dev/null
###Check the optimization
geo_min=$(get_geom_mopac.sh ${molecule}_freq.out)
if [ "$geo_min" = "Error" ];then
   echo "The input structure could not be optimized. Check your XYZ file"  
   exit
else
   echo "$geo_min"  > opt_start.xyz
fi

if [ $postp_alg -ge 1 ]; then
###########min0 is the reference minimum. If it exists, take its energy
   if [ -f $tsdirll/MINs/min.db ]; then
      e0=$(sqlite3 ${tsdirll}/MINs/min.db "select energy from min where name='min0_0'")
   else
      if [ "$program_opt" = "mopac" ]; then
         e0=$(awk '/FINAL HEAT OF FORMATION =/{e0=$6};END{print e0}' ${molecule}_freq.out )
      else
###if LL program is gaussian, perform frequency with gaussian to get e0
         nameg=${molecule}_freq_gauss
         chk="$(echo %chk=$nameg)"
         cal="$(sed 's/ts,noeigentest,//;s/tkmc/'$temperature'/;s@level1@'$method_opt'@;s/charge/'$charge'/;s/mult/'$mult'/;s@iop@@' $sharedir/hl_input_template)"
##we take geo from XYZ file
         geo="$(awk 'NF==4{print $0};END{print ""}' ${molecule}.xyz)"
         ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo")"
         echo -e "$ig09\n\n" > ${nameg}.dat
         g09 <${nameg}.dat >${nameg}.log
         ok=$(awk 'BEGIN{ok=0};/Frequencies/{++nfreq;if($3>0 && $4>0 && nfreq==1) ok=1};END{print ok}' ${nameg}.log)
         if [ $ok -eq 0 ];then
            echo "The input structure could not be optimized with gaussian. Check your XYZ file"  
            exit
         fi
         e0=$(get_energy_g09_${LLcalc}.sh ${nameg}.log 1)         
         emaxts=$(echo "scale=6; $emaxts/627.51" | bc | awk '{printf "%14.6f",$1}')
      fi
   fi
###########3
   emaxts=$(echo "scale=6; $emaxts+$e0" | bc | awk '{printf "%14.6f",$1}')
fi
if [ $mdc -ge 1 ]; then itrajn=5 ; fi
###Loop over the trajectories
for i in $(seq 1 $itrajn) 
do 
  named=${molecule}_dyn${i}
  echo ""
##Empty temporary folders
  rm -rf partial_opt/* ts_opt/* 
####
#This is only for internal dynamics (MOPAC)
####
  if [ $mdc -ge 1 ] && [ $i -gt 1 ]; then
     echo "Searching for reactive events with:"
  else
     echo ""
     echo "+-+-+-+-+-+-+-+-+-          Trajectory ${i}          +-+-+-+-+-+-+-+-+-"
     if [ $md -eq 0 ]; then
       echo "Performing BXDE MD"
       bxde.py $inputfile &>  ${named}.log
       if [ ! -f traj.xyz ]; then
          echo "traj.xyz does not exist"
          continue 
       else
          mv traj.xyz coordir/${named}.xyz
          mv bond_order.txt coordir/${named}.bo
       fi
     elif [ $md -eq 1 ]; then 
        echo "Performing standard MD"
        echo "$dytem1"     > ${named}.mop
        initialqv_mopac_samp1.sh ${molecule}_freq.out $seed $excite $nlms $lstnm | nm.exe >> ${named}.mop
        mopac ${named}.mop &> ${named}.log 
        if [ ! -f ${named}.xyz ]; then
           echo "${named}.xyz does not exist"
           continue 
        else
           mv ${named}.xyz coordir
        fi
     elif [ $md -eq 2 ]; then
        echo "Performing standard MD"
        echo "$dytem1"     > ${named}.mop
        initialqv_mopac_samp2.sh ${molecule}_freq.out $excite $nlms $lstnm $thmass | termo.exe | sed 's/ 1.d69/1.d69/g' >> ${named}.mop
        mopac ${named}.mop &> ${named}.log
        if [ ! -f ${named}.xyz ]; then
           echo "${named}.xyz does not exist"
           continue 
        else
           mv ${named}.xyz coordir
        fi
     else
       echo "Reading external dynamics results from coordir"
     fi
  fi

  if [ $postp_alg -eq 0 ]; then
     echo "End of traj "$i
     echo "Only trajs. No post-processing algorithm applied to search for TSs"
     break
  elif [ $postp_alg -eq 1 ]; then
     postp_file=bbfs     
  elif [ $postp_alg -eq 2 ]; then
     postp_file=bots     
  fi
###########
#From here everything is common for internal and external dynamics
###########
  if [ $i -eq 1 ]; then
     echo "  *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*"  > ${postp_file}.out
     echo "                $postp_file algorithm results          " >> ${postp_file}.out
     echo "  *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*" >> ${postp_file}.out
     echo ""                                                        >> ${postp_file}.out
  fi
  echo "  Trajectory $i" >> ${postp_file}.out

  if [ $postp_alg -eq 1 ]; then
     if [ $mdc -eq 0 ]; then
        snapshots_mopac.sh coordir/${named}.xyz  $irange | bbfs.exe >> ${postp_file}.out
     elif [ $mdc -ge 1 ]; then
        namedc=${molecule}_dyn1
        irange=$((20-4*(i - 1) ))
        irangeo2=$(echo "scale=0; $irange/2" | bc )
        echo "Time window (fs) = $irange "
        snapshots_mopac.sh coordir/${namedc}.xyz  $irange | bbfs.exe >> ${postp_file}.out
     fi
  elif [ $postp_alg -eq 2 ]; then
     bots.py $natom $cutoff $stdf ${named} >> ${postp_file}.out
  fi

  path=$(awk '/Number of paths/{np=$4};END{print np}' ${postp_file}.out )
  if [ $path -eq 0 ]; then
     echo "This traj has no paths "
     continue
  fi

  echo "Npaths=" $path
  chapath[0]=0
  for ip in $(seq $path)
  do
# Find the highest energy point
    if [ $postp_alg -eq 1 ]; then
       ijc=$(awk '/Joint path=/{if($2=='$ip') ijc=$5};END{print ijc}' ${postp_file}.out)
    else
       ijc=0
    fi
    jp=$((ip - 1))
##If previous path was multiple, continue 
    chapath[$ip]=$ijc
    if [ ${chapath[$jp]} -eq 1 ]; then continue ; fi
##
    if [ $ijc -eq 0 ]; then
       echo "Path" $ip" (Single): $nppp attempt(s)  to locate the ts" 
       ll=$((wrkmode-1))
       dlt=$((wrkmode+1))
       ul=1
    elif [ $ijc -eq 1 ]; then 
       echo "Path" $ip" (Multiple): several attempts to locate the ts" 
       ll=$((1 - irangeo2))
       dlt=$(echo "scale=2; $irange/6" | bc | awk '{print int($1+0.5)}')
       ul=$irangeo2  
    fi
    npo=0
    for itspt in $(seq $ll $dlt $ul)
    do 
       npo=$((npo + 1))
       ctspt=$((100*ip + irangeo2 + itspt))
       echo "$min_template"         > partial_opt/pes$ctspt
       if [ $postp_alg -eq 1 ]; then
          cat partial_opt/fort.$ctspt >> partial_opt/pes$ctspt
       else
          cat partial_opt/fort.$ip >> partial_opt/pes$ctspt
       fi
       mopac partial_opt/pes$ctspt  2>/dev/null
       geo_pes=$(get_geom_mopac.sh partial_opt/pes${ctspt}.out)
       if [ "$geo_pes" = "Error" ]; then continue ; fi

       name=ts${i}_${ip}_${ctspt}
       if [ "$program_opt" = "mopac" ]; then
          echo "$ts_template"                      > ts_opt/$name
          echo "$geo_pes" | awk 'NF==4{print $0}' >> ts_opt/$name
          echo "$freq_template"                   >> ts_opt/$name
          mopac ts_opt/$name 2>/dev/null
          file=ts_opt/${name}.out
#If too many variables, run ts int
          if [ $(awk 'BEGIN{f=0};/Too many variables/{f=1};END{print f}' $file) -eq 1 ]; then
             sed -i 's/ts /ts int /g' ts_opt/$name
             mopac ts_opt/$name 2>/dev/null
          fi
###
          prog=1
       else
#construct g09 input file
          chk="$(echo %chk=ts_opt/ts$name)"
          cal="$(sed 's/tkmc/'$temperature'/;s@iop@@;s@level1@'$method_opt'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/hl_input_template)"
          geo="$(echo "$geo_pes" | awk 'NF==4{print $0};END{print ""}')"
          ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo")"
          echo -e "$ig09\n\n" > ts_opt/${name}.dat
          g09 <ts_opt/${name}.dat >ts_opt/${name}.log && gauss_line=$(echo $LINENO)
          file=ts_opt/${name}.log
          ok=$(awk 'BEGIN{ok=0};/Frequencies/{++nfreq;if($3<0 && $4>0 && nfreq==1) ok=1};END{print ok}' $file)
          if [ $ok -eq 1 ]; then
             get_energy_g09_${LLcalc}.sh $file 1   > tmp_gauss
             get_freq_g09.sh $file >> tmp_gauss
             prog=2
          else
             printf "     Pt%2s: failed-->EF algorithm was unable to optimize a TS\n" $npo
	     continue    
          fi
       fi

       fe="$(mopac_freq_ts.sh $file $prog)"
       fi="$(echo "$fe" | awk '{printf "%10.0f",$1}')"
       ei="$(echo "$fe" | awk '{printf "%14.6f",$2}')"
       if [[ ("$fi" -eq -1) ]]; then
          printf "     Pt%2s: failed-->Lowest real freq is negative\n" $npo
          continue
       elif [[ ("$fi" -eq -2) ]]; then
          printf "     Pt%2s: failed-->Sum of 2 lowest real freqs < 10cm-1\n" $npo
          continue
       elif [[ ("$fi" -eq -3) ]]; then
          printf "     Pt%2s: failed-->Stationary point is a minimum\n" $npo
          continue
       elif [[ ("$fi" -eq -4) ]]; then
          printf "     Pt%2s: failed-->EF algorithm was unable to optimize a TS\n" $npo
          continue
       elif (( $(echo "$ei > $emaxts" |bc -l) )); then
          printf "     Pt%2s: TS optimized but not added-->E=%4s kcal/mol > %4s kcal/mol\n" $npo $ei $emaxts
          continue
       fi
       if [[ ("$fi" -ge "$imag") ]]; then
          string="$(echo "$fe" | awk '{printf "%10.0f %10.4f %10.0f %10.0f %10.0f %10.0f",$1,$2,$3,$4,$5,$6}')"
# GLB added lock to tslist so that duplicate numbers are not created
          (
          flock -x 200 || exit 1
          if [ -f "$tslistll" ]; then
             ok=$(diff.sh $string $tslistll $prog)
             if [[ ("$ok" -eq "-1") ]]; then
                nt=$(awk '{nt=$2};END{print nt + 1}' $tslistll )
                name=ts${nt}_${nb}
                printf "ts%5s%18s%70s traj= %4s Path= %10s\n" $nt $name "$string" $i $nb  >> $tslistll
                cp ${file} $tsdirll/${name}.out
                printf "     Pt%2s: TS optimized and added to ts list\n" $npo
             else
                printf "     Pt%2s: TS optimized but not added-->redundant with ts %4s\n" $npo $ok
             fi
          else
             nt=1
             name=ts${nt}_${nb}
             printf "ts%5s%18s%70s traj= %4s Path= %10s\n" $nt $name "$string" $i $nb  >> $tslistll
             cp ${file} $tsdirll/${name}.out
             printf "     Pt%2s: TS optimized and added to ts list\n" $npo
          fi
          ) 200>>${tslistll}.lock
          if [ $mdc -ge 1 ]; then exit ; fi
          break
       else
          printf "     Pt%2s: TS optimized but not added-->imag=%4si cm-1 < %4si cm-1\n" $npo $fi $imag
       fi
    done
  done
done


if [ $sampling -ne 30 ]; then 
   echo ""
   echo "END OF THE CALCULATIONS" 
   echo ""
fi

