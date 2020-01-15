#!/bin/bash
kmcfile=$1
minn=$2
en=$3
awk 'BEGIN{m[0]='$minn';n=0;k=0}
{if($10=="MIN" && $5<'$en'){
  i=0
  while(i<=n){
    if( $8==m[i] ) {++k;m[k]=$11;break}
    if($11==m[i] ) {++k;m[k]=$8;break}
    i++
    }
  n=k
  }
}
END{i=0
while(i<=n){
   j=0
   ok=1
   while(j<i){
      if(m[j]==m[i]) ok=0
      j++
      }
   if(ok==1) print m[i]
   i++
   }
}' $kmcfile >tmp01
cnom[0]=$(wc -l tmp01 | awk '{print $1}')
for i in $(seq 1 10)
do
   output="$(awk '{if(NR==FNR) m[NR]=$1;n0=NR;n=n0}
   {if($10=="MIN" && $5<'$en'){
     i=1
     while(i<=n){
       if( $8==m[i] ) {++kk;k=n0+kk;m[k]=$11}
       if($11==m[i] ) {++kk;k=n0+kk;m[k]=$8}
       i++
       }
     n=n0+k
     }
   }
   END{i=1
   while(i<=n){
      j=0
      ok=1
      while(j<i){
         if(m[j]==m[i]) ok=0
         j++
         }
      if(ok==1) print m[i]
      i++
      }
   }' tmp01 $kmcfile)"
   echo "$output" > tmp01
   cnom[$i]=$(wc -l tmp01 | awk '{print $1}')
   if [ ${cnom[$i]} -eq ${cnom[$i-1]} ]; then 
#      echo "convergence after $i iterations"
      break
   fi
done
cat tmp01 $kmcfile |  awk '{if(NF==1) {++k;m[k]=$1;n=k}}
{if(NF>1)  {
   ++nol
   if(nol<=2) print $0
   if(nol>2){
      i=1
      while(i<=n){
        if($8==m[i] && $5<'$en') {print $0;break}
        if($10=="MIN" && $11==m[i] && $5<'$en') {print $0;break}
        i++ 
        }
      }
  } 
}'  

##

