#!/bin/bash
source utils.sh

sharedir=${AMK}/share
cwd=$PWD
inputfile='amk.dat'

#on exit remove tmp files
tmp_files=(tmp* filedum frag* atsdum* code* cop ddf* dum* Frag* isteps ConnMat labels mingeom nf ScalMat sprint* )
trap cleanup EXIT INT

exe=$(basename $0)

rm -f isteps
awk '{if(NR!=2) print $0};{if(NR==2) print ""}' $1 > filedum 
awk 'BEGIN{print "Labels"};NF==4{print $1}' $1 > labels
file=filedum
##
read_input
##

elements=${sharedir}/elements
natom=$(awk 'NR==1,NR==1{print $1}' $file )

createMat.py $file 1 $nA

echo "1 $natom" | cat - ConnMat | cat - ConnMat | sprint2.exe > sprint.out

awk '/Results for the Laplacian/{
getline
sumt=$6+$7+$8
sumd=$6+$7
sumu=$6
if(sumt<=0.001) 
  print "3"
else if (sumd<=0.001)
  print "2"
else if (sumu<=0.001)
  print "1"
}' sprint.out > nf
nfrag=$(awk '{print $1}' nf)
##############   
awk '/Natom/{natom=$2};/Adjacency matrix/{getline
i=1
while(i<=natom){
  print $0
  getline
  i++
  }
}' sprint.out >ddf

for i in $(seq 1 10)
do
   ((j=i-1))
   awk '{for(i=1;i<=NF;i++) {v[NR,i]=$i;if($i>0) {++l[NR];ik[NR,l[NR]]=i}}}
   END{
   j=1
   code=0
   idum=0
   while(j<=NF){
     k=1
     while(k<=l[j]){
       for(i=1;i<=NF;i++) v[j,i]=v[j,i]+v[ik[j,k],i]
       k++
       }
     i=1
     while(i<=NF){
       if(i>j) ++idum
       iexp=NF*(NF-1)-idum
       if(v[j,i]>1) v[j,i]=1
       code+=v[j,i]^iexp
       i++
       }
     for(i=1;i<=NF;i++) {printf "%2.0f ",v[j,i]};print ""
     j++
     }
   print code > "code"'$i'
   }' ddf >ddf.out
   code[$i]=`awk '{print $1 }' code$i`
   if [ -f code$j ]; then
      if [ ${code[$i]} -eq  ${code[$j]} ] ;then
        echo "$i steps" >> isteps
        break
      fi
   fi
   cp ddf.out ddf
done

awk '{for(i=1;i<=NF;i++) {if($i>0) printf "%3.0f ",i};print ""
}' ddf >ddf.log
awk '{for(i=1;i<=NF;i++) {noanr[NR]=NF;v[NR,i]=$i}}
END{
print "number of fragments=",'$nfrag'
i=1
while(i<=NR){
  j=1
  p=1
  while(j<=i-1){
    if(v[i,1]==v[j,1]) p=0
    j++
    }  
  if(p==1) ++nn
  if(p==1 && nn<'$nfrag') {print "Number of atoms of this fragment=",noanr[i]
  for(j=1;j<=noanr[i];j++) {printf "%3.0f ",v[i,j]};print "" }
  i++
  }
}' ddf.log >dum.log
##############   

awk '/Natom/{natom=$2;print natom}
/Adja/{getline
while(i<=natom){
  if(NF==natom) print $0
  getline
  i++
  }
}' sprint.out >dum
   
cat dum labels dum.log >dummi

awk '{if(NR == FNR) {z[NR]=NR;l[NR]=$1;n[NR]=NR/10+1;tne=NR}}
{if(NR>FNR && FNR==1) natom=$1}
{if(NR>FNR && NF==natom)  {
  ++i
  for(j=1;j<=natom;j++)  a[i,j]=$j
  }
}
/of atoms/{
++nf
nat=$7
noa[nf]=nat
getline
for(j=1;j<=nat;j++) {at[nf,j]=$j}
}
/Labels/{
i=1
while(i<=natom) {
   getline
   IGNORECASE = 1
   j=1
   while(j<=tne){
     if($1 == l[j]) z[i]=j 
     j++  
     } 
     print i,$1,z[i],n[z[i]]
   i++
   }
}
END{
for(inf=1;inf<=nf;inf++){
for(j=1;j<=noa[inf];j++){
  kk=1 
  while(kk<=natom) {
     sum+=n[z[at[inf,j]]]*a[at[inf,j],kk]*n[z[kk]]
     kk++
     }
  }
}
  print "Code=",sum*'$nfrag'
  print "Code=",sum*'$nfrag' >> "cop"
} ' $elements dummi >>dum.log



awk 'BEGIN{nrfh=1d20}
/number of fragments/{nfrag=$4}
/Number of atoms of this fragment/{++i;natom[i]=$7;getline
for(j=1;j<=NF;j++) at[i,j]=$j
nrfh=NR
}
{
if(NR>nrfh) {
k=1
while(k<=1000){
   pr=1
   if($1== "Code=") break
   l=1
   while(l<=nfrag-1){
     m=1
     while(m<=natom[l]){
       if($1 == at[l,m]) pr=0
       m++
       }
     l++
     }
   if(pr==1) {++natom[nfrag];at[nfrag,natom[nfrag]]=$1}
   k++
   getline
   }
 }
}
END{
i=1
while(i<=nfrag){
  j=1
  while(j<=natom[i]){
    print at[i,j] > "Frag"i
    j++
    }
  i++
  }
}' dum.log
rm -f frag.log
for i in $(seq 1 $nfrag)
do
   wc -l Frag$i | awk '{print $1}' > frag$i.xyz
   echo "" >> frag$i.xyz
   cat Frag$i > dum
   echo "break" >>dum
   cat $file >> dum
   awk '/break/{fl=1}
   {if(fl==0) {++n;at[n]=$1}}
   {if(fl==1 && NF==4) {
     ++nat
     j=1
     while(j<=n){
       if(nat==at[j]) print $0
       j++
       }
     }
   }' dum >> frag$i.xyz
   chemical_formula.sh frag$i.xyz >frag$i.log 
   cat frag$i.log >>frag.log
done

awk '/natom=/{natom=$2;getline;print natom,$1}' frag.log >dum
sort -n dum >dumo.log
awk '{f[NR]=$2}
END{
i=1
while(i<=NR){
  if(i<NR) printf "%s + ",f[i]
  if(i==NR) printf "%s",f[i]
  i++
  }
print ""
}' dumo.log
