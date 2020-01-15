correlfile=catalysis_res/KMC/Correl_with_barrierless
sed 's/=/= /g' $correlfile > correlfile
rxnetfile=catalysis_res/KMC/RXNet_catalysis_with_barrierless
rm -f la_react_pro
ls catalysis_res/KMC/kmc_files/kmc*.out > kmcf
set `awk '{print $1}' kmcf`
for i
do
   echo $i
   awk '/counts per process/{ok=1}
   {if(ok==1 && NF==2 && $NF>0) {++tot;suma+=$2;lal[tot]=$1;lac[tot]=$2}
   }
   END{if(ok==0) exit
   i=1
   while(i<=tot){
      p=lac[i]/suma*100
      if(p>0) print lal[i],lac[i]
      i++
      }
   }' $i >> la_react_pro0
done
#remove duplicates in la_react_pro
awk '{++i;lal[i]=$1;lac[i]=$2;tot=i}
END{i=1
while(i<=tot){
   ok=1
   j=1
   while(j<i){
     if(lal[i]==lal[j]) {ok=0;break}
     j++
     }
   if(ok==1) print lal[i],lac[i]
   i++
   }
}' la_react_pro0 >la_react_pro
cat la_react_pro correlfile > ccf
awk '{if(NF==2) {++j;la[j]=$1;tot=j}
if($1!="Subsystem:" && NF>2) {
  ++i
  k=1
  while(k<=tot){
    if(i==la[k]) print $(NF-1),$NF
    k++
    }
  }
}' ccf >ccf.log
cat ccf.log $rxnetfile > rrcf
awk '{if($1=="TS" && NF==2) {++ts;tsl[ts]=$NF}
if($1=="Gprod=" && NF==2) {++pr;prl[pr]=$NF}
}
/Subsystem/{print $0}
{if($1=="TS" && NF>2){
  ok=0
  i=1  
  while(i<=ts){
     if($2==tsl[i]) {ok=1;break}
     i++
     }
  j=1
  while(j<=pr){
     if($NF==prl[j]) {ok=1;break}
     j++
     } 
  if(ok==1) print $0
  }
}' rrcf > catalysis_res/KMC/RXNet_active
