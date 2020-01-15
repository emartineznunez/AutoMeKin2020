awk 'BEGIN{nif=0};/Frequencies/{++j;for(i=3;i<=NF;i++) {l=i-2;k=(j-1)*3+l;freq[k]=$i}}
END{nfreq=k
k=1
while(k<=nfreq) {
  if(freq[k]<0) ++nif
  k++
  }
print nif
}' $1

