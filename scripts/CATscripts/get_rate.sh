file=$1
cn=$2
maxtime=$3
code_prod=$4
maxtime=2
rm -f k
awk 'BEGIN{t0=1}
/Calculation number/{cn=$3}
/Time=/{t=$2
i=1
while(i<=100){
   getline
   if($1=="'$code_prod'" && cn=='$cn' && t>=t0 && t<='$maxtime') print t,$2
   if(NF==0) break
   i++
   }
}' $file > k
rmv=`awk 'BEGIN{rmv=0}
{t=$1}
END{
if(t<0.95*'$maxtime') rmv=1
print rmv
}' k`
if [ $rmv -eq 1 ]; then
   rm k
fi
