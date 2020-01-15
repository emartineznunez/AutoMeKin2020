#!/bin/bash
source utils.sh
#On exit remove tmp files
tmp_files=(ConnMat deg* fort.* intern* mingeom ScatMat ts_tors* ScalMat *_opt.* tors.* geom*)
trap 'err_report2 $LINENO $gauss_line' ERR
trap cleanup EXIT INT

#Paths and names
cwd=$PWD
sharedir=${AMK}/share
exe=$(basename $0)

if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi

###Reading stuff from inputfile
read_input
###

if [ $sampling -ge 30 ];then
   echo "No torsion calculations for this sampling"
   exit
fi

###check whether we do a scan
kmcfilell=$tsdirll/KMC/RXNet_long.cg_groupedprods
minfilell=$tsdirll/MINs/SORTED/MINlist_sorted
confilell=$tsdirll/working/conf_isomer.out
factor=1.5

if [ -f $minfilell ] && [ -f $kmcfilell ]; then
   echo "Select the relevant minima to find conformational isomers"
   minn=$(awk '/min0/{print $2}' $minfilell)
   if [ ! -f $confilell ]; then
      echo "File $confilell does not exist and we cannot proceed"
      exit
   fi
   minok=$(awk 'BEGIN{min='$minn'}
   {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
   j=1
   while(j<=iso[NR]){
      if('$minn'==m[NR,j]) min=m[NR,1]
      j++
      }
   }
   END{print min}' $confilell )
else
   echo "Files $minfilell and/or $kmcfilell do not exist and we cannot proceed"
   exit
fi
###

en=$(awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile )

#Compute the absolute values of emax and emin 
if [ -f $tsdirll/MINs/min.db ]; then
   e0=$(sqlite3 ${tsdirll}/MINs/min.db "select energy from min where name='min0_0'")
else
   if [ "$program_opt" = "mopac" ]; then
      e0=$(awk '/FINAL HEAT OF FORMATION =/{e0=$6};END{print e0}' ${molecule}_freq.out )
   else
      name=$molecule"_freq_gauss"
      e0=$(get_energy_g09_${LLcalc}.sh ${name}.log 1)
      emaxts=$(echo "scale=6; $emaxts/627.51" | bc | awk '{printf "%14.6f",$1}')
   fi
fi

echo The value of e0=$e0
echo The relative value of emaxts=$emaxts
emaxts=$(echo "scale=6; $emaxts+$e0" | bc | awk '{printf "%14.6f",$1}')
echo The absolute value of emaxts=$emaxts

flag="tors"
for min in $(get_minn.sh $kmcfilell $minok $en $factor | awk 'NF==1{print $1}')
do
   echo "Running the tors calculations for  MIN $min"
   names="MIN"$min
   echo "$bo_template"            > tors_opt.mop
   sqlite3 $tsdirll/MINs/SORTED/mins.db "select geom from mins where name='$names'" >> tors_opt.mop
   echo "$freq_template"         >> tors_opt.mop
   mopac tors_opt.mop 2>/dev/null
   get_geom_mopac.sh tors_opt.out  >mingeom 
   createMat.py mingeom 1 $nA
   deg="$(awk '{si=0;for(i=1;i<=NF;i++) {si+=$i;if(i==NF) print si,$0 }}' ConnMat)" 
   bo="$(awk '/BOND ORDERS AND VALENCIES/{getline;getline;getline;getline
   for(i=1;i<='$natom';i++) {getline; for(j=3;j<=i+1;j++) {if($j<2.0) print i,j-2} } }' tors_opt.out)"
   echo "$deg" > deg_bo
   echo "$bo" >> deg_bo
   awk '{
   if(NR<='$natom'){
     deg[NR]=$1
     l=0
     if(deg[NR]>1) {for(i=2;i<=NF;i++) {if($i==1) {++l;jatom[NR,l]=i-1 } }}
     }
   else {
     if(deg[$1]>1 && deg[$2]>1) {
       ok=0
       j=1
       while( j<=deg[$1] ){
         if(jatom[$1,j] != $2) k=jatom[$1,j]
         else ok=1
         ++j
         }
       j=1
       while( j<=deg[$2] ){
         if(jatom[$2,j] != $1) {l=jatom[$2,j];break}
         ++j
         }
       if(ok==1) print k,$0,l
       }
     }
   }' deg_bo  > deg_bo.out
   
   ntor=$(wc -l deg_bo.out | awk '{print $1}')
   echo "Number of torsions $ntor"
   if [ $ntor -eq 0 ]; then exit ; fi
   for itor in $(awk '{print NR}' deg_bo.out)
   do
      echo "Running the TS search for tors $itor of molecule tors"
      lr=$(awk 'NR=='$itor'{print $1;exit}' deg_bo.out )
      l1=$(awk 'NR=='$itor'{print $2;exit}' deg_bo.out )
      l2=$(awk 'NR=='$itor'{print $3;exit}' deg_bo.out )
      l3=$(awk 'NR=='$itor'{print $4;exit}' deg_bo.out )
      if [ $l1 -gt $lr ]; then ((l1=l1-1)) ; fi
      if [ $l2 -gt $lr ]; then ((l2=l2-1)) ; fi
      if [ $l3 -gt $lr ]; then ((l3=l3-1)) ; fi
      get_geom_mopac.sh tors_opt.out >intern.dat
      awk 'NR=='$itor'{print $0}' deg_bo.out >>intern.dat
      intern.exe <intern.dat>intern.out
      fok=$(awk 'BEGIN{fok=1};/Abort/{fok=0};END{print fok}' intern.out)
      if [ $fok -eq 0 ]; then continue ; fi
      internlastatom="$(awk '{print $1,$2,$3,$4,$5,$6,$7,'$l1','$l2','$l3'}' intern.out)"
      sed 's/method/'"$method"' charge='$charge'/g' $sharedir/path_template >tors.mop
      get_geom_mopac.sh tors_opt.out | awk 'NF==4{print $0}' | awk '{if(NR!= '$lr') print $0}' >> tors.mop 
      echo "$internlastatom" >> tors.mop
      mopac tors.mop 2>/dev/null
      awk '/VARIABLE        FUNCTION/{++i;getline
      e[i]=$3
      getline
      getline
      getline
      j=1
      while(j<='$natom'){
        getline
        l[i,j]=$0
        j++
        }
      }
      END{
      i0=2
      while(i0<=i){
        if(e[i0]>e[i0-1] && e[i0]>e[i0+1]) {
           k=1
           proc=1
           while(k<=nmax){
             diff0=e[i0]-emax[k]
             diff=sqrt(diff0*diff0)
             if(diff<=0.01) proc=0
             k++
             } 
           if(proc==0) {i0++;continue}
           ++nmax
           emax[nmax]=e[i0]  
           j=1
           while(j<='$natom'){
             print l[i0,j]
             j++
             }
           }
        i0++
        } 
      }' tors.out>geomts_tors0
      nomax="$(awk 'END{print NR/'$natom'}' geomts_tors0)"
      echo "$nomax possible candidate(s)"
      for inmax in $(seq 1 $nomax)
      do
         echo "$inmax th candidate(s)"
         awk 'BEGIN{nr0='$natom'*('$inmax'-1)+1;nrf='$natom'*'$inmax'};{if(NR>=nr0 && NR<=nrf) print $0}' geomts_tors0 >geomts_tors
         name="ts_tors"$itor
         echo "$ts_template"     > ${name}.mop
         awk '$7=1' geomts_tors >> ${name}.mop
         echo "$freq_template"  >> ${name}.mop
         mopac ${name}.mop 2>/dev/null
         if [ "$program_opt" = "mopac" ]; then
            file=${name}.out
            prog=1
         else
#construct g09 input file
            chk="$(echo %chk=$name)"
            cal="$(sed 's/tkmc/'$temperature'/;s@iop@@;s@level1@'$method_opt'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/hl_input_template)"
            geo="$(get_geom_mopac.sh ${name}.out | awk '{if(NF==4) print $0};END{print ""}')"
            ig09="$(echo -e "$chk"'\n'"$cal"'\n'"$geo")"
            echo -e "$ig09\n\n" > ${name}.dat
            g09 <${name}.dat >${name}.log && gauss_line=$(echo $LINENO)
            file=${name}.log
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

   #check the ts
         fe="$(mopac_freq_ts.sh $file $prog)" 
         fi="$(echo "$fe" | awk '{printf "%10.0f",$1}')"
         ei="$(echo "$fe" | awk '{printf "%14.6f",$2}')"
         if [[ ("$fi" -eq -1) ]]; then
            echo "TS opt failed. Lowest real freq < 0 cm-1"
            continue
         elif [[ ("$fi" -eq -2) ]]; then
            echo "TS opt failed. Sum of two Lowest real freqs < 50 cm-1"
            continue
         elif [[ ("$fi" -eq -3) ]]; then
            echo "TS opt failed. Stationary point is a minimum"
            continue
         elif [[ ("$fi" -eq -4) ]]; then
            echo "TS opt failed. Stationary point search crashed"
            continue
         elif (( $(echo "$ei > $emaxts" |bc -l) )); then
            echo "TS has an energy $ei, which is greater than" $emaxts
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
                  printf "ts%5s%18s%70s traj= %4s Path= %10s\n" $nt $name "$string" $flag $nb  >> $tslistll
                  cp ${file}  $tsdirll/${name}.out 
                  printf "     Pt%2s: TS optimized and added to ts list\n" $itor
                  continue
               else
                  printf "     Pt%2s: TS optimized but not added-->redundant with ts %4s\n" $itor $ok
                  continue
               fi
            else
               nt=1
               name=ts${nt}_${nb}
               printf "ts%5s%18s%70s traj= %4s Path= %10s\n" $nt $name "$string" $flag $nb  >> $tslistll
               cp ${file}  $tsdirll/${name}.out
               printf "     Pt%2s: TS optimized and added to ts list\n" $itor
            fi
            ) 200>>${tslistll}.lock
         else
            printf "     Pt%2s: TS optimized but not added-->imag=%4si cm-1, which is lower than %4si cm-1\n" $itor $fi $imag
         fi
      done
   done
done

