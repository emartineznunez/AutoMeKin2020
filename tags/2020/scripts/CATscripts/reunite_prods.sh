prfile=prods.log
cp catsum.log dcatsum
cp ratsum.log dratsum
cp correlsum.log dcorrelsum
nprod=`wc -l $prfile | awk '{print $1}'`
for i in $(seq 1 $nprod)
do
   prodi=`awk 'NR=='$i',NR=='$i'{print $1}' $prfile`
   echo "checking $prodi"
   gl=`awk '{
   if($10=="'$prodi'" || $12=="'$prodi'") {
     if($13>$15) 
        print $13
     else
        print $15
     exit}
   }' dcatsum `
   if [ -z $gl ]; then 
      echo "No products like this"
      continue  
   fi
   awk '{
   if($10=="'$prodi'" || $12=="'$prodi'") {
     ++i
     tot=i
     if($13>$15) 
       la[i]=$13
     else
       la[i]=$15}
   }
   END{i=2
   while(i<=tot){
    print la[i]
    i++
    }
   }' dcatsum >bl 
   echo $gl
   nsubs=`wc -l bl | awk '{print $1}'`
   for j in $(seq 1 $nsubs)
   do
      bl=`awk 'NR=='$j',NR=='$j'{print $1}' bl`
      sed 's/ '$bl'/ '$gl'/g' dcatsum > ddcatsum
      sed 's/ '$bl'/ '$gl'/g' dratsum > ddratsum
      sed 's/ '$bl'/ '$gl'/g' dcorrelsum > ddcorrelsum
      cp ddcatsum dcatsum 
      cp ddratsum dratsum 
      cp ddcorrelsum dcorrelsum 
   done
done 
