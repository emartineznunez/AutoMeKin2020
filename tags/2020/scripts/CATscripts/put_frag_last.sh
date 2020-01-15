nm=$1
cat filemin frag$nm".xyz" > compose
awk 'NR==1,NR==1{
natom=$1
getline
i=1
while(i<=natom) {
  getline
  atom_old[i]=$0
  i++
  }
getline
natomf=$1
getline
i=1
while(i<=natomf){
   getline
   atomf_new[i]=$0
   jnew=natom-natomf+i
   i++
   }
}
END{i=1
sfrga=natom-natomf
while(i<=natom){
   ok=1
   j=1
   while(j<=natomf){
     if(atom_old[i]==atomf_new[j]) {ok=0;break}
     j++
     }
   if(ok==1) atom_new[i]=atom_old[i]
   if(ok==0 && i>sfrga) atom_new[i]=atom_old[i] 
   if(ok==0 && i<=sfrga) {++i1;la1[i1]=i} 
   if(ok==1 && i>sfrga) {++i2;la2[i2]=i} 
   i++
   }
i=1
while(i<=i1){
  atom_new[la1[i]]=atom_old[la2[i]]
  atom_new[la2[i]]=atom_old[la1[i]]
  i++ 
  }
i=1
while(i<=natom){
   print atom_new[i]
   i++
   }
}' compose
