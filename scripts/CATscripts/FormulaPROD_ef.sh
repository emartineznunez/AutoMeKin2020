rm -r frag*.xyz frag*.log
rm -f isteps
awk '{if(NR!=2) print $0};{if(NR==2) print ""}' $1 > filedum 
file=filedum
inputfile=amk.dat
vdw=`awk 'BEGIN{vdw=0};{if($1=="vdw") vdw=1};END{print vdw}' $inputfile `
if [ $vdw -eq 1 ]; then
   cp thdist thdist_backup
   nvdw=`awk '{if($1=="vdw") print $2}' $inputfile `
#   echo $nvdw "Van der Waals distances to be taken into account"
   for i in $(seq 1 $nvdw)
   do
      at1[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $1; i++}  }' $inputfile `
      at2[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $2; i++}  }' $inputfile `
      dis[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i') print $3; i++}  }' $inputfile `
#      echo "Distance $i between atoms ${at1[$i]} and ${at2[$i]} is ${dis[$i]}"
      awk '{if($1=="'${at1[$i]}'" && $2=="'${at2[$i]}'") {print $1,$2,"'${dis[$i]}'"}
      else if($2=="'${at1[$i]}'" && $1=="'${at2[$i]}'") {print $1,$2,"'${dis[$i]}'"}
      else {print $0}
      } ' thdist >thdist_vdw
      cp thdist_vdw thdist
   done
fi


#echo $file

sed "s/'/ /g" atsymb >atsdum1
sed "s/,/ /g" atsdum1 > atsdum2
awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' atsdum2  > atsdum2.out

echo "Labels" >labels
awk '{if(NF==4) print $0 }' $file > mingeom
awk '{if(NF==4) print $1 }' $file >> labels


natom=`awk 'NR==1,NR==1{print $1}' $file `
echo "1" $natom > sprint.dat

 

createMat.sh
cat ConnMat >> sprint.dat
createMat.sh
cat ConnMat >> sprint.dat
sprint2.exe <sprint.dat >sprint.out

### if vdw=1 restore thdist
if [ $vdw -eq 1 ]; then
   cp thdist_backup thdist
fi
###


awk '/Results for the Laplacian/{
getline
sum6=$6+$7+$8+$9+$10+$11
sum5=$6+$7+$8+$9+$10
sum4=$6+$7+$8+$9
sum3=$6+$7+$8
sum2=$6+$7
if(sum6<=0.001) 
  print "6"
else if(sum5<=0.001) 
  print "5"
else if(sum4<=0.001) 
  print "4"
else if(sum3<=0.001) 
  print "3"
else if(sum2<=0.001)
  print "2"
else
  print "1"
}' sprint.out > nf
nfrag=`awk '{print $1}' nf`

##############   
awk '/Natom/{natom=$2};/Adjacency matrix/{getline
i=1
while(i<=natom){
  for(j=1;j<=natom;j++) {
   if(i==j) 
     el=1
   else
     el=$j
   printf "%s ",el }
  printf "\n"
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
} ' atsdum2.out dummi >>dum.log



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

awk '/natom=/{natom=$2;getline;print $1}' frag.log 
