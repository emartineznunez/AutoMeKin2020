changes=$1
input=$2
output=$3
awk '{if( NR == FNR) {c1[NR]=$1;c2[NR]=$2;tc=NR}}
{if(NR >FNR) {
 ch=0
 for(i=1;i<=tc;i++){if($1==c1[i]){ch=1;c=c2[i]} }
 if(ch==0) print $0
 if(ch==1) print c 
 }
}' $changes $input > $output
