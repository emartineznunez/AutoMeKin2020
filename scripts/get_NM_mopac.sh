file=$1
output=$2
echo "[Molden Format]" > ${output}.freq
echo "[FREQ]"         >> ${output}.freq
##echo the freqs here
get_freq_mopac.sh $file >> ${output}.freq
echo "       " >> ${output}.freq
echo "[FR-COORD]" >> ${output}.freq

awk 'BEGIN{ifreq=0;atobohr=1.889726}
/Empirical Formula/{natom=$(NF-1)}
/ORIENTATION OF MOLECULE IN FORCE CALCULATION/{
getline
getline
getline
i=1
while(i<=natom){
 getline
 x[i]=$3;y[i]=$4;z[i]=$5
 print $2,$3*atobohr,$4*atobohr,$5*atobohr  >>  "'${output}'.freq"
 if(i>1) {
   dx=x[i]-x[i-1]
   dy=y[i]-y[i-1]
   dz=z[i]-z[i-1]
   mod=sqrt(dx*dx+dy*dy+dz*dz)
   ddx[i]=dx/mod;ddy[i]=dy/mod;ddz[i]=dz/mod
 }
 if(i>2) {
   prod=ddx[i]*ddx[i-1]+ddy[i]*ddy[i-1]+ddz[i]*ddz[i-1] 
   sp[i]=sqrt( prod*prod )
   sump+=sp[i]
   diff=natom-2-sump
 }
 i++ 
}
print "" >> "'${output}'.freq"
print "[FR-NORM-COORD]" >> "'${output}'.freq"
if(natom==2 || diff <0.001 ) 
  nlin=1
else
  nlin=0
ntimes=int((3*natom-6)/8+1)
}
/NORMAL COORDINATE ANALYSIS/{
itimes=1
while(itimes<=ntimes){
  i=1
  while(i<=7){
    getline
    ++i
    }
  for(inf=1;inf<=NF;++inf) {++ifreq;freq[ifreq]=$inf} 
  getline
  i=1
  while(i<=3*natom){
   getline
   for(inf=2;inf<=NF;++inf) {nfreq=ifreq-NF+inf;nm[i,nfreq]=$inf} 
   ++i
   }
 ++itimes
 }
}
END{i=1
while(i<=ifreq){
 print "Vibration "i >> "'${output}.freq'" 
   j=1
   while (j<=natom) {
   print nm[3*j-2,i]*atobohr,nm[3*j-1,i]*atobohr,nm[3*j,i]*atobohr >> "'${output}'.freq"
   ++j
   }
 ++i
 }
}' $file 


