file=mingeom0
natoma=$1
ncm=$2
natom=`wc -l $file | awk '{print $1}'`
awk 'BEGIN{natom='$natom';natoma='$natoma';natomb=natom-natoma;i=0;j=0;count=0;watot=0;wbtot=0;ncm='$ncm'}
{if( NR == FNR) l[NR]=$1;w[l[NR]]=$2}
{if(NF==4) ++count}
{if(count>0 && count<=natoma) {j++
   na[j]=$1;xa[j]=$2;ya[j]=$3;za[j]=$4;wa[j]=w[$1];watot+=wa[j]}
else if(count>natoma) {i++
   nb[i]=$1;xb[i]=$2;yb[i]=$3;zb[i]=$4;wb[i]=w[$1];wbtot+=wb[i]}
}
END{i=1
xcma=0
ycma=0
zcma=0
while(i<=natoma){
  xcma+=wa[i]*xa[i]
  ycma+=wa[i]*ya[i]
  zcma+=wa[i]*za[i]
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
  xcmb+=wb[i]*xb[i]
  ycmb+=wb[i]*yb[i]
  zcmb+=wb[i]*zb[i]
  i++
  }
xcmb=xcmb/wbtot
ycmb=ycmb/wbtot
zcmb=zcmb/wbtot
dx=xcma-xcmb
dy=ycma-ycmb
dz=zcma-zcmb
dcm=sqrt(dx^2+dy^2+dz^2)
dx=dx/dcm
dy=dy/dcm
dz=dz/dcm
i=1
while(i<=natoma){
   xa[i]=xa[i]-xcma
   ya[i]=ya[i]-ycma
   za[i]=za[i]-zcma
   printf "%4s %20.10f %20.10f %20.10f\n",na[i],xa[i],ya[i],za[i]
   printf "%4s %20.10f %20.10f %20.10f\n",na[i],xa[i],ya[i],za[i] >"fragtoop"
   i++
   }
i=1
while(i<=natomb){
   xb[i]=xb[i]-xcmb-ncm*dx*dcm
   yb[i]=yb[i]-ycmb-ncm*dy*dcm
   zb[i]=zb[i]-zcmb-ncm*dz*dcm
   printf "%4s %20.10f %20.10f %20.10f\n",nb[i],xb[i],yb[i],zb[i]
   i++
   }
}' atsdum3.out $file
