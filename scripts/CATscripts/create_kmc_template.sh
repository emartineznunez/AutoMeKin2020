inputfile=amk.dat
kmc_files=catalysis_res/KMC/kmc_files
if [ -d $kmc_files ]; then
   echo "$kmc_files folder exists"
else
   echo "Create $kmc_files folder"
   mkdir $kmc_files
fi

kmcfile=catalysis_res/KMC/kmc_template.dat
rxnfile=catalysis_res/KMC/RXNet_catalysis
ratsum=catalysis_res/KMC/Rates_with_barrierless
nmol=`awk '{if($1=="nmol") print $2}' $inputfile`
nespconst=`awk '{if($1=="nespconst") print $2}' $inputfile`
awk '{if($1!="Subsystem:")print $0}' $ratsum > catalysis_res/KMC/rates_with_barrierless
echo "kmc calculation" > $kmcfile
nproc=`wc -l catalysis_res/KMC/rates_with_barrierless | awk '{print $1}'`
nspec=`awk 'BEGIN{l=0};{if($2>l) l=$2;if($3>l) l=$3;if($4>l) l=$4;if($5>l)l=$5};END{print l}' catalysis_res/KMC/rates_with_barrierless`
echo "$nproc $nspec 3" >> $kmcfile
nic=`awk '{if($1=="species") print NF-1}' $inputfile`
echo $nic
for i in $(seq 1 $nic)
do
    iic[$i]=`awk '{if($1=="c0") {++i;if(i=='$i') {print $2;exit} }}' $inputfile`
    awk '{if($1=="c0") {++i;if(i=='$i') {print $0;exit} }}' $inputfile >dumc0
    c00[$i]=`awk '{if(NF==4) print $3;if(NF==5) print $3*$5}' dumc0`
    c0[$i]=`awk '{for(i=3;i<=NF;i++) print $i}' dumc0`
    if [ -d tsdirHL_${iic[$i]} ]; then
       indx0=`awk '/min0/{print $2}' tsdirHL_${iic[$i]}/MINs/SORTED/MINlist_sorted `
       isu=`awk '{if($1=="Subsystem" && $3=="'${iic[$i]}'") print $2}' $rxnfile `
       flag=$indx0"_"$isu
       indxf[$i]=`awk '{if($6=="'$flag'") {print $7;exit};if($11=="'$flag'") {print $12;exit}}' $rxnfile `
    else
       indxf[$i]=`awk '{if(NR>2 && $2=="'${iic[$i]}'") {print $1;exit}}' $rxnfile`
    fi 
done
min_conc=${c00[1]}
for i in $(seq 1 $nic)
do
   if [ 1 -eq "$(echo "${c00[$i]} < $min_conc" | bc)" ]; then
       min_conc=${c00[$i]}
   fi
done
echo $min_conc >min_conc
vol=`awk '{print '$nmol'/$1/6.022e23}' min_conc`
echo $vol >> $kmcfile
echo $nespconst >> $kmcfile
for i in $(seq 1 $nespconst)
do
    ec[$i]=`awk '{if($1=="nespconst") print $(2+'$i')}' $inputfile`
    echo $i ${ec[$i]}
    awk '{if(NR>2 && $2=="'${ec[$i]}'") {print $1;exit}}' $rxnfile>> $kmcfile
done
cat catalysis_res/KMC/rates_with_barrierless >> $kmcfile

n=0
for i in $(seq 1 $nspec)
do
   ok=0
   for j in $(seq 1 $nic)
   do
     if [ $i -eq ${indxf[$j]} ]; then
        ok=1
        jn=$j
     fi
   done
   if [ $ok -eq 0 ]; then
      echo "0" >> $kmcfile
   else
      ((n=n+1))
      echo "C0_$n "${c0[$jn]} >> $kmcfile
   fi
done
timec=`awk '{if($1=="Maxtime") print $2}' $inputfile`
dtc=`awk '{if($1=="Stepsize") print $2}' $inputfile`
echo $timec  $dtc >>$kmcfile

