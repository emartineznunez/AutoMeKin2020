sed "s/'/ /g" atsymb >atsdum1
sed "s/,/ /g" atsdum1 > atsdum2
awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' atsdum2  > atsdum2.out
awk 'BEGIN{huge=1000000}
{if( NR == FNR) l[NR]=$1}
/Charge =/{p1=NR}
/orientation:/{ getline
getline
getline
getline
i=1
while(i<=huge){
  getline
  if(NF==1) break
  la[i]=l[$2]
  x[i]=$4
  y[i]=$5
  z[i]=$6
  natom=i
  i++
  }
}
END{
i=1
while(i<=natom){
   print la[i],x[i],y[i],z[i]
   i++
   }
print ""
}' atsdum2.out $1
