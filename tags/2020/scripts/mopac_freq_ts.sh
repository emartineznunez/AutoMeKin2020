file=$1
if [ $2 -eq 1 ]; then
   awk 'BEGIN{zero=0}
   /Empirical Formula/{natom=$7;natmone=natom-1}
   /FINAL HEAT OF FORMATION =/{e=$6;f=-4}
   /MASS-WEIGHTED COORDINATE ANALYSIS/{
   getline; getline; getline; getline; getline; getline; getline;
   imagf=$1
   f=sqrt(imagf*imagf)
   f1=$2
   f2=$3
   f3=$4
   f4=$5
   s12=f1+f2
   if(f1<0) f=-1
   if(s12<10) f=-2
   if(imagf>0) f=-3
   }
   END{
   if(f==0) f=zero
   if(e==0) {e=zero;f=-4}
   print f,e,f1,f2,f3,f4}' $file
else
   awk 'BEGIN{zero=0}
   NR==1{ e=$1}
   NR>1{++i;freq[i]=$1
   }
   END{
   imagf=freq[1]
   f=sqrt(imagf*imagf)
   f1=freq[2]
   f2=freq[3]
   f3=freq[4]
   f4=freq[5]
   s12=f1+f2
   if(f1<0) f=-1
   if(s12<10) f=-2
   if(imagf>0) f=-3
   if(f==0) f=zero
   if(e==0) {e=zero;f=-4}
   print f,e,f1,f2,f3,f4}' tmp_gauss
fi

