dir=barrierless_diss
file=$dir/barrierless.dat
catsum=catalysis_res/KMC/RXNet_catalysis
ratsum=catalysis_res/KMC/Rates
correlsum=catalysis_res/KMC/Correl
cp $ratsum ratsum
cp $catsum catsum
cp $correlsum correlsum
lastkmc=`awk 'BEGIN{max=0};{if(NF==5) {if($2>max)max=$2;if($3>max)max=$3;if($4>max)max=$4;if($5>max)max=$5}};END{print max}' ratsum  `
if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "read catalysis file is missing"
   exit
fi
eta=`awk '{if($1=="eta") print $2}' amk.dat`
tkmc=`awk '{if($1=="TKMC") print $2}' amk.dat`
noHLcalc=`awk '/HLcalc/{print $2}' $inputfile`
if [ $noHLcalc -eq 1 ]; then
   HLcalc1=`awk '/HLcalc/{print $3}' $inputfile`
   level1=`awk '{if($1=="level1") print $2}' $inputfile`
   basis1=`awk '{if($1=="basis1") print $2}' $inputfile`
   echo "HL using one level only for energies and frequencies" $HLcalc1
   HLcalc=$HLcalc1
elif [ $noHLcalc -eq 2 ]; then
   HLcalc1=`awk '/HLcalc/{print $3}' $inputfile`
   HLcalc2=`awk '/HLcalc/{print $4}' $inputfile`
   level2=`awk '{if($1=="level2") print $2}' $inputfile`
   basis2=`awk '{if($1=="basis2") print $2}' $inputfile`
   echo "HL using two levels" $HLcalc1 "and" $HLcalc2
   HLcalc=$HLcalc2
else
   echo "check the input file "
fi

awk '/leading/{print $0}' $file > dumigb
if [ ! -f $catsum ]; then
   echo "Please, run catalysis.py first"
   exit
fi
sed 's/MIN//g' dumigb | sed 's/_/ /g'  >dumigb.log
nt=0
set `awk '{print NR}' dumigb` 
for i
do 
   ((nt=nt+1))
   subs=`awk 'NR=='$i',NR=='$i'{print $4}' dumigb`
###
   natom=`get_n_atoms.sh $dir/diss$nt.log`
###
   natoma=`get_n_atoms.sh $dir/fragtoop$nt.log`
   ((natomb=natom-natoma))
   cn=$nt
   min=`awk 'NR=='$nt',NR=='$nt'{print $1}' dumigb.log`
   if [ ! -f $dir/down$cn.log ];then
      echo "Please, run Barrierless_CAT.sh first"
      exit 
   fi
   flb=`get_dist.sh $dir/down$cn.log $natoma`
#only the barrierless channels ($flb -eq 1)  are analyzed
   if [ $flb -eq 1 ]; then
      echo $natoma > geoma
      echo "" >> geoma
      get_lastgeomg09.sh  $dir/fragtoop$cn".log" >>geoma
      formul_a=`FormulaMOL.sh geoma`
      where=`awk '/Subsystem/{if($4=="'$formul_a'") print $3}' $catsum`
      formul_b=`awk 'NR=='$nt',NR=='$nt'{print $NF}' dumigb`
      formul_ab=`echo "$formul_b+$formul_a"`
      subi=`awk '/Subsystem/{if($3=="'$subs'") {print $2;exit}}' $catsum` 
      minkmc=`awk '/Subsystem/{if($3=="'$subs'") {subi=$2;ok=1}}
      {for(i=1;i<=NF;i++) if($i=="'$min'_"subi && ok==1 && i>2) {print $(i+1);exit}}' $catsum`
      gmin=`awk '/Subsystem/{if($3=="'$subs'") {subi=$2;ok=1}}
      {for(i=1;i<=NF;i++) if($i=="'$min'_"subi && ok==1 && i>2) {print $(i+2);exit}}' $catsum`

      if [ ! -z $where ]; then
         echo "Where= $where" > matchmin2.dat
         echo $HLcalc $noHLcalc >> matchmin2.dat
         fragmin=`matchmin2.sh matchmin2.dat $dir/fragtoop$cn".log"`
         if [ $fragmin -gt 0 ]; then
            fragminkmc=`awk '/Subsystem/{if($3=="'$where'") {subi=$2;ok=1}}
            {for(i=1;i<=NF;i++) if($i=="'$fragmin'_"subi && ok==1 && i>2) {print $(i+1);exit}}' $catsum`
###the min might not be in RXNet    because it might be connected by a higher enery barrier
            if [ -z $fragminkmc ]; then
               ((fragminkmc=lastkmc+1))
               lastkmc=$fragminkmc
            fi
         else
            ((fragminkmc=lastkmc+1))
            lastkmc=$fragminkmc
         fi
      else
######emilio
         ((fragminkmc=lastkmc+1))
         lastkmc=$fragminkmc
      fi
      get_energy$HLcalc.sh $dir/fragtoop$cn".log" $noHLcalc >gfragmin
      get_gcorr.sh $dir/fragtoop$cn".log" >>gfragmin
      nif=`get_NIF.sh $dir/fragtoop$cn".log" `
      if [ $nif -gt 0 ]; then
         echo "Calc. number $cn barrierless dissociation. One fragment has $nif imaginary frequencies...skip it" 
      fi
      awk '/Reference energy:/{printf "%20.10f\n",$3}' $catsum >gpr
      awk '/Subsystem/{if($2=='$subi') printf "%20.10f\n",$NF}' catalysis_res/KMC/Energies >>gpr
      awk '{e[NR]=$1};END{printf "%20.10f\n",e[1]+e[2]}' gfragmin >>gpr
      awk '{if($1=="'$formul_b'") printf "%20.10f\n",$2}' prods.out >>gpr
      gpr=`awk '{e[NR]=$1};END{print (e[3]+e[4]+e[2]-e[1])*627.51}' gpr`
#this is the minimum}     
      frb=$formul_b
      frbkmc=`awk '{if(NR>2 && NF==3 && $3=="'$frb'") {print $1;exit} }' $catsum`
      if [ -z $frbkmc ]; then
         ((frbkmc=lastkmc+1))
         lastkmc=$frbkmc   
      fi
      echo $gpr >dif
      echo $gmin >>dif
      deltag=`awk 'BEGIN{dif=0};{e[NR]=$1};END{printf "%8.3f\n",e[1]-e[2]}' dif`
#      dif=`awk 'BEGIN{dif=0};{e[NR]=$1};END{diff=e[1]-e[2];if(diff>0) dif=1;print dif}' dif`
#      if [ $dif -eq 0 ]; then
#         echo "Calc. number $cn barrierless dissociation. DeltaG is negative ( $deltag  )...skip it" 
#      fi
#      awk '/Reference energy:/{printf "%20.10f\n",$3}' $catsum >gpr
#      if [ $dif -eq 1 ] && [ $nif -eq 0 ] ; then
      if [ $nif -eq 0 ] ; then
         echo "Calc. number $cn barrierless dissociation from min: (orig: $min of $subs) $minkmc $gmin <--> $fragminkmc + $frbkmc  $gfragmin " 
         echo $eta $deltag $tkmc > diff_and_diss.dat
         diff_and_diss.exe <diff_and_diss.dat>diff_and_diss.out     
         kdiff=`awk 'NR==1,NR==1{print $1}' diff_and_diss.out`
         kdiss=`awk 'NR==2,NR==2{print $1}' diff_and_diss.out`
         awk '{if($1=="Subsystem:" && $2=="'$subi'") {
            print $0
	    printf "%30s %4s %4s %4s    0\n",'$kdiff','$fragminkmc','$frbkmc','$minkmc' 
	    printf "%30s %4s    0 %4s %4s\n",'$kdiss','$minkmc','$fragminkmc','$frbkmc'
            }
         else
	    print $0
         }' ratsum > ratsum.log
         cp ratsum.log ratsum
###
         awk '{if($1=="Subsystem:" && $2=="'$subi'") {
            print $0
	    printf "%30s %4s %4s %4s    0 Gprod=%10s\n",'$kdiff','$fragminkmc','$frbkmc','$minkmc','$gpr' 
	    printf "%30s %4s    0 %4s %4s Gprod=%10s\n",'$kdiss','$minkmc','$fragminkmc','$frbkmc','$gpr'
            }
         else
	    print $0
         }' correlsum > correlsum.log
         cp correlsum.log correlsum
###
         awk '{if($1=="Subsystem" && $2=="'$subi'") add=1}
         {if($1=="TSlabel" && add==1)  {
            print $0
	    printf "TS   null        null ==>    MIN    '$min'_'$subi'    '$minkmc'       '$gmin'    <-->   '$formul_a' + '$formul_b'            '$frbkmc' + '$fragminkmc'         '$gpr'\n"
            add=0
            }
         else
	    print $0
         }' catsum > catsum.log
         cp catsum.log catsum
      fi
   else
      echo "Calc. number $cn barrierless dissociation presumably has a barrier...skip it"
   fi
done
reunite_prods.sh
cp dratsum $ratsum"_with_barrierless"
cp dcatsum $catsum"_with_barrierless"
cp dcorrelsum $correlsum"_with_barrierless"
rm -f ratsum* catsum* dratsum ddratsum dcatsum ddcatsum correlsum*
