#swap over hydrogens h1 and h2
h1=$1
h2=$2
awk '{++nat
if(NR=='$h1') 
   p1=$0
else if(NR=='$h2')
   p2=$0
else {
   p[nat]=$0
   }
}
END{
k=1
while(k<=nat){
  if(k=='$h1') 
     print p2
  else if(k=='$h2') 
     print p1
  else
     print p[k]
  k++
  }
}' mingeom0 >dum.out
cp dum.out mingeom0
