ind=$1
jnd=$2
if [ -f "ScalMat" ]; then rm ScalMat ; fi
awk '
NR==FNR{++i;label[i]=$1;x[i]=$2;y[i]=$3;z[i]=$4}
NR>FNR{++j;label1[j]=$1;label2[j]=$2;rsc[j]=$3}
END{
IGNORECASE = 1
i0=1
while(i0<=i){
 cm[i0,i0]=0
 ii=i0+1
 while(ii<=i){
  cm[i0,ii]=0
  xij=x[i0]-x[ii]
  yij=y[i0]-y[ii]
  zij=z[i0]-z[ii]
  r[i0,ii]=sqrt(xij*xij+yij*yij+zij*zij)
  k=1
  l=0
  while(k<=j){
   dum1=label[i0]
   dum2=label[ii]
   if(dum1==label1[k] && dum2==label2[k]) {l=k}
   if(dum1==label2[k] && dum2==label1[k]) {l=k}
   k++
   }
   print rsc[l] >> "ScalMat"
   if(r[i0,ii]<1.1*rsc[l]) cm[i0,ii]=1
   if(i0=='$ind' && ii=='$jnd') cm[i0,ii]=0
   cm[ii,i0]=cm[i0,ii]
  ii++
  }
 i0++
 }
i0=1
while(i0<=i){
 for(j0=1;j0<=i;j0++) { printf "%s ",cm[i0,j0]}; print ""
 i0++
 }
}' mingeom thdist > ConnMat

