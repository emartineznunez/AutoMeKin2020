kmcfile=$1
minn=$2
subs=$3
ene=$4
sed 's/_/ /g' $kmcfile > kmcfile
en=`awk '{if($7=="'$minn'" && $8=="'$subs'") e0=$10}
{if($13=="'$minn'" && $14=="'$subs'") e0=$16}
END{print '$ene'+e0} ' kmcfile`

awk 'BEGIN{m[0]='$minn';n=0;k=0}
/Subsystem/{if($2=='$subs') ok=1}
/Subsystem/{if($2!='$subs') ok=0}
{if($12=="MIN" && ok==1 && $4<'$en'){
  i=0
  while(i<=n){
    if( $7==m[i] ) {++k;m[k]=$13;break}
    if($13==m[i] ) {++k;m[k]=$7;break}
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
}' kmcfile >minn0
cp minn0 minn
cnom[0]=`wc -l minn | awk '{print $1}'`
for i in $(seq 1 10)
do
   awk '{if(NR==FNR) m[NR]=$1;n0=NR;n=n0}
   /Subsystem/{if($2=='$subs') ok=1}
   /Subsystem/{if($2!='$subs') ok=0}
   {if($12=="MIN" && ok==1 && $4<'$en'){
     i=1
     while(i<=n){
       if( $7==m[i] ) {++kk;k=n0+kk;m[k]=$13;break}
       if($13==m[i] ) {++kk;k=n0+kk;m[k]=$7;break}
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
   }' minn kmcfile >minn$i
   cp minn$i minn
   cnom[$i]=`wc -l minn | awk '{print $1}'`
   if [ ${cnom[$i]} -eq ${cnom[$i-1]} ]; then
      break
   fi
done
