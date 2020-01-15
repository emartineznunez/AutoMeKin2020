awk 'BEGIN{huge=10^10}
/Charge =/{i=1
while(i<=huge){
   getline
   if(NF==4) ++natom
   if(NF==0) break
   i++
   }
}
END{print natom}' $1
