g09file=$1
energy=$2
zpe=$3
gcorr=$4
label=$5
awk 'BEGIN{huge=1000000;zero=0}
{if( NR == FNR) l[NR]=$1}
/orientation:/{ getline
getline
getline
getline
i=1
while(i<=huge){
  getline
  if(NF==1) break
  n[i]=$2
  x[i]=$4
  y[i]=$5
  z[i]=$6
  natom=i
  i++
  }
}
/Frequencies/{++j;for(i=3;i<=NF;i++) {ll=i-2;k=(j-1)*3+ll;freq[k]=$i}}
END{
print natom
printf "%2s %15s %4s %8s %5s %14s\n","E=","'$energy'","ZPE=","'$zpe'","Gcorr=","'$gcorr'"
i=1
while(i<=natom){
  print l[n[i]],x[i],y[i],z[i]
  print l[n[i]],x[i]+6*('$label'-1),y[i],z[i]  >> "catalysis_res/geome_"'$label'
  i++
  }
nfreq=k
k=1
while(k<=nfreq) {
  if(freq[k]>0) printf "%5.0f\n",freq[k]
  if(freq[k]>0) printf "%5.0f\n",freq[k]  >> "catalysis_res/frequ_"'$label'
  k++
  }
}' atsdum2.out catalysis_res/$g09file".log" > catalysis_res/$g09file".rxyz"


