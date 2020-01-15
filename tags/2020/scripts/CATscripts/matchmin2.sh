file=$1
path2=`awk '/Where=/{print "tsdirHL_"$2}' $file` 
HLcalc=`awk 'NR==2,NR==2{print $1}' $file`
noHLcalc=`awk 'NR==2,NR==2{print $2}' $file`
fraglog=$2
awk '{if(NR>2) {
	j=0
	p=1
	while(j<=i){
            if($8==min[j]) p=0		
            j++		
            } 
	if(p==1) {++i;min[i]=$8}
	if($10=="MIN") {
		j=0
		p=1
		while(j<=i){
                   if($11==min[j]) p=0		
		   j++
		   }	
		if(p==1) {++i;min[i]=$11}
		}
	}
}
END{tot=i
j=1
while(j<=tot){
   print min[j]
   j++
   }
}' $path2/KMC/RXNet_long.cg_groupedprods >minlist
energy=`get_energy$HLcalc.sh $fraglog $noHLcalc`
ls $path2/MINs/SORTED/MIN*.rxyz > minlisttot
sed 's/MIN/ /g' minlisttot | sed 's/_min/ /g' > dumi 
n=0
minok=0
set `awk '{print $0}' minlisttot `
for i
do
    ((n=n+1))
    echo $energy >comp
    awk 'NR==2,NR==2{print $2}' $i >>comp
    ediff=`awk 'BEGIN{pm=0};{e[NR]=$1};END{diff=(e[2]-e[1])*627.51;diff=sqrt(diff*diff);if(diff<0.005) pm=1;print pm}' comp`
    min=`awk 'NR=='$n',NR=='$n'{print $(NF-1)}' dumi`
    if [ $ediff -eq 1 ]; then
# see if the min is in minlist
	minok=`awk 'BEGIN{min=0}
        {if('$min'==$1) min=$1
        }
        END{
        print min
        }' minlist`
# check the isomers
        if [ $minok -eq 0 ]; then
           f=$path2/working/conf_isomer.out
	   minok=`awk 'BEGIN{min=0}
           {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
           j=1
	   while(j<=iso[NR]){
              if('$min'==m[NR,j]) min=m[NR,1]
              j++ 
              }
           }
           END{print min}' $f `
        fi
        echo $minok
	exit
    fi
done
echo $minok


