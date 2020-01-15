file=$1
natoma=$2
natom=`get_n_atoms.sh $file`
awk 'BEGIN{huge=10^10;natom='$natom';natoma='$natoma';natomb=natom-natoma;i=0;j=0;watot=0;wbtot=0}
{if( NR == FNR) l[NR]=$1;w[NR]=$2}
/ orientation:/{++ngm;getline
getline
getline
getline
i=1
k=0
while(i<=huge){
  getline
  if(NF==6 && i<=natoma) {
      xa[i,ngm]=$4;ya[i,ngm]=$5;za[i,ngm]=$6
      if(ngm==1) {na[i]=$2;wa[i]=w[$2];watot+=wa[i]}}
  if(NF==6 && i>natoma)  {
      ++k;xb[k,ngm]=$4;yb[k,ngm]=$5;zb[k,ngm]=$6
      if(ngm==1) {nb[k]=$2;wb[k]=w[$2];wbtot+=wb[k]}}
  if(NF==1) break
  i++
  }
}
END{
igm=1
while(igm<=ngm){
   i=1
   xcma=0
   ycma=0
   zcma=0
   while(i<=natoma){
     xcma+=wa[i]*xa[i,igm]
     ycma+=wa[i]*ya[i,igm]
     zcma+=wa[i]*za[i,igm]
     i++
     }
   xcma=xcma/watot
   ycma=ycma/watot
   zcma=zcma/watot
   i=1
   xcmb=0
   ycmb=0
   zcmb=0
   while(i<=natomb){
     xcmb+=wb[i]*xb[i,igm]
     ycmb+=wb[i]*yb[i,igm]
     zcmb+=wb[i]*zb[i,igm]
     i++
     }
   xcmb=xcmb/wbtot
   ycmb=ycmb/wbtot
   zcmb=zcmb/wbtot
   dx=xcma-xcmb
   dy=ycma-ycmb
   dz=zcma-zcmb
   dcm[igm]=sqrt(dx^2+dy^2+dz^2)
   igm++
   }
diff=dcm[1]-dcm[ngm]
if(diff>0) print "1"
if(diff<0) print "0"
}' atsdum3.out $file
