###########
#This script finds reactions that proceed without a barrier:
# for instance:
# Here, we look for shallow minima (of any type that can dissociate)
#it puts both fragments some distance appart and then it runs downhill calculation
#This script performs g09 calculations and the second part does the analysis
##########
if [ $# -eq 0 ]; then
   ci=1
elif [ $1 -eq 0 ]; then
   ci=0
else
  ci=1
fi

bd=barrierless_diss
if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk.dat file is missing"
   exit
fi
if [ -d $bd ];then
   echo "$bd already exists"
   exit
else
   mkdir $bd
fi



catsum=catalysis_res/KMC/RXNet_catalysis
awk '{if(NR>2){
  if(NF==3) {
    print $3
    }
   else
    exit
  }
}' $catsum >prods0.log
vdw=`awk 'BEGIN{vdw=0};{if($1=="vdw") vdw=1};END{print vdw}' $inputfile `
if [ $vdw -eq 1 ]; then
   cp thdist thdist_backup
   nvdw=`awk '{if($1=="vdw") print $2}' $inputfile `
   echo $nvdw "Van der Waals distances to be taken into account"
   for i in $(seq 1 $nvdw)
   do
      at1=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $1; i++}  }' $inputfile `
      at2=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $2; i++}  }' $inputfile `
      dis=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i') print $3; i++}  }' $inputfile `
      echo "Distance $i between atoms ${at1} and ${at2} is ${dis}"
      awk '{if($1=="'${at1}'" && $2=="'${at2}'") {print $1,$2,"'${dis}'"}
      else if($2=="'${at1}'" && $1=="'${at2}'") {print $1,$2,"'${dis}'"}
      else {print $0}
      } ' thdist >thdist_vdw
      cp thdist_vdw thdist
   done
fi


echo "Summary of the calcs" > $bd/barrierless.dat
echo "********************" >> $bd/barrierless.dat
TKMC=`awk '{if($1 == "TKMC") {print $2;nend=1}};END{if(nend==0) print "298"}' $inputfile`
charge=`awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' $inputfile `
mult=`awk 'BEGIN{mu=1};{if($1=="mult") mu=$2};END{print mu}' $inputfile `
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

sed 's/level1/'$level1'\/'$basis1'/g' hl_input_template | sed 's/charge/'$charge'/g' | sed 's/mult/'$mult'/g'  > g09inp0
sed 's/tkmc/'$TKMC' calc/g' g09inp0 | sed 's/opt=(ts,noeigentest,calcall,noraman)//g'  >g09inp

sed "s/'/ /g" atsymb >atsdum1
sed "s/,/ /g" atsdum1 > atsdum2
sed "s/+00/+00 /g" atsdum2 > atsdum3
awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' atsdum2  > atsdum2.out
awk '/ams=/{amlab=1}
/character/{++nt}
{for(i=1;i<=(NF-1);i++) {if($1 != "real" && amlab==1 && nt==0) {++j;m[j]=$i} }}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) {++k;print $i,m[k]}
} ' atsdum3 > atsdum3.out


awk '{if(NF==4) print $0}' cat.xyz >mingeom0
met_label=`awk '{if( NR == FNR) {l[NR]=$1;tne=NR}}
{if(FNR > 1 && NR >FNR ) {
   IGNORECASE = 1
   i=1
   while(i<=tne){
      if( $1 == l[i] && i>=21 && i<=30) print $1 
      if( $1 == l[i] && i>=39 && i<=48) print $1 
      if( $1 == l[i] && i>=72 && i<=80) print $1 
      i++
      }
  }
}' atsdum2.out mingeom0`
awk '{if($1=="'$met_label'") {print NR-2;exit}
}' cat.xyz >atoms_to_remove
nmet_label=`awk '{print $1}' atoms_to_remove`
rm -f minn minn.log
subsystems.py >subsystems.out
#Getting other possibilities from tsdirHL
n=0
tc=0
#Search for possible products besides the prods0.log
if [ $ci -eq 1 ]; then
   echo "possible prods"> prods
   set `awk '{print $0}' subs`
   for i
   do
      prodir="tsdirHL_"$i"/PRODs/CALC"
      if [ -d $prodir ]; then
         set `ls $prodir/*.log | awk '{print $1}' ` 
         for j
         do
           pr=`basename $j`
           echo $pr >pr
           sed 's/\./ /g' pr | sed 's/\-/ /g' | sed 's/m//g'  >pr.log
           prod=`awk '{if($1 !~/'$met_label'/) print $1}' pr.log`
           multi=`awk '{print $3}' pr.log`
           if [ $multi -eq 1 ] && [ ! -z $prod ] ; then 	
              get_energy$HLcalc.sh $j $noHLcalc >en
              get_gcorr.sh $j >>en
              en=`awk '{e[NR]=$1};END{printf "%20.10f",e[1]+e[2]}' en`
              echo $prod $en >> prods
           fi
         done
      fi
   done
   awk '{if(NR>1) {++i;p[i]=$1;tot=i;pe[i]=$0}}
   END{i=1
   while(i<=tot){ok=1
      j=1
      while(j<i){
        if(p[i]==p[j]) {ok=0;break} 
        j++
        }
      if(ok==1) print pe[i]
      i++
      }
   }' prods >prods.out
   paste prods0.log prods.out > prod_comp
   awk '{if(NF==3) {++i;ni[i]=$1;tot=i;p[NR]=$2;e[NR]=$3}}
   {if(NF==2) p[NR]=$1;e[NR]=$2}
   END{i=1
   while(i<=NR){ok=1
      j=1
      while(j<=tot){
        if(p[i]==ni[j]) {ok=0;break}
        j++
        }
      if(ok==1) print p[i],e[i]
      i++
      }
   }' prod_comp>prods.log
fi

nprod=`wc -l prods.log | awk '{print $1}'`
nprod0=`wc -l prods0.log | awk '{print $1}'`
nc=0
tc=0
set `awk '{print $0}' subs`
for i
do
   ((nc=nc+1))
   echo "" >> $bd/barrierless.dat
   echo "System: $i" >> $bd/barrierless.dat
   echo "*******" >> $bd/barrierless.dat
   ((n=n+1)) 
# ${frag[$j]}
   minfile="tsdirHL_"$i"/MINs/SORTED/MINlist_sorted"
   mindir="tsdirHL_"$i"/MINs/SORTED"
   confile="tsdirHL_"$i"/working/conf_isomer.out"
   minn=`awk '/min0/{print $2}' $minfile` 

   res=`awk 'BEGIN{res=0}
   {for(i=1;i<=NF;i++) if($i=='$minn') res=$1}
   END{
   print res
   }' $confile`
   if [ $res -gt 0 ]; then
      minn=$res
   fi
   get_minn_fromoverall.sh $catsum $minn $nc 20
   set `awk '{print $0}' minn`
   for j
# Loop over the structures (most important minima) of a given subsystem
   do
      min=$mindir/MIN$j"_*.rxyz"
      min_name=`basename $min`
      awk '{if(NF==4) print $0}' $min >mingeom0
      natom=`wc -l mingeom0 | awk '{print $1}'`
      echo $natom > filemin
      echo "" >> filemin
      cat mingeom0 >> filemin
#remove only the metal
      cat atoms_to_remove mingeom0 > pmg
      awk '{if(NF==1) {++i;lab[i]=$1;tot=i}}
      {if(NF==4){ok=1;++k
        i=1
        while(i<=tot){
          if(k==lab[i]) {ok=0;break}
          i++
          }
        if(ok==1) print $0
        }
      }' pmg >mingeom
#mingeom does not contain catalyzer now

      createMat.sh
#      ConnMat Look for fragments and choose appropriate ones
# now identify two or more fragments (if there is more than one)
      awk 'END{print NR"\n"}' mingeom> geom
      cat mingeom >> geom
      FormulaPROD_ef.sh geom > molec
      nmolec=`wc -l molec | awk '{print $1}'` 
      for inmp0 in $(seq 1 $nprod0)
      do 
          m[$inmp0]=1
      done
      for nm in $(seq 1 $nmolec)
      do
         molec=`awk 'NR=='$nm',NR=='$nm'{print $1}' molec`
######First, the initial prods (prods0)
         for nmp0 in $(seq 1 $nprod0)         
         do 
            prodi0=`awk 'NR=='$nmp0',NR=='$nmp0'{print $1}' prods0.log`
            if [ "$molec" == "$prodi0" ] && [ ${m[$nmp0]} -eq 1 ] ; then
               m[$nmp0]=0
	       ((tc=tc+1))
	       natomB=`awk 'NR==1,NR==1{print $1}' frag$nm.xyz`
               ((natomA=natom-natomB))
	       echo "$min_name of subsystem $i leading to $prodi0"
	       echo "$min_name of subsystem $i leading to $prodi0" >>$bd/barrierless.dat
               if [ $noHLcalc -eq 1 ]; then
	          balecalc.sh $natom $natomA $tc $bd $nm $nmet_label 1 $charge $mult
               elif [ $noHLcalc -eq 2 ]; then
	          balecalc.sh $natom $natomA $tc $bd $nm $nmet_label 2 $charge $mult $level2 $basis2
               fi
               queueHLoptdown_sub.sh $bd  diss$tc lanza$tc $tc 
               queueHL_sub.sh $bd  fragtoop$tc lanzaop$tc $tc
            fi
         done
######Now, other prods (prods)
         for nmp in $(seq 1 $nprod)         
         do 
            prodi=`awk 'NR=='$nmp',NR=='$nmp'{print $1}' prods.log`
            if [ "$molec" == "$prodi" ]; then
	       ((tc=tc+1))
	       natomB=`awk 'NR==1,NR==1{print $1}' frag$nm.xyz`
               ((natomA=natom-natomB))
	       echo "$min_name of subsystem $i leading to $prodi"
	       echo "$min_name of subsystem $i leading to $prodi" >>$bd/barrierless.dat
               if [ $noHLcalc -eq 1 ]; then
	          balecalc.sh $natom $natomA $tc $bd $nm $nmet_label 1 $charge $mult
               elif [ $noHLcalc -eq 2 ]; then
	          balecalc.sh $natom $natomA $tc $bd $nm $nmet_label 2 $charge $mult $level2 $basis2
               fi
               queueHLoptdown_sub.sh $bd  diss$tc lanza$tc $tc 
               queueHL_sub.sh $bd  fragtoop$tc lanzaop$tc $tc
            fi
         done

      done 
   done
done



