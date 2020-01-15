rm rotate.out
for i in $(seq 1 100) 
do 
   echo "18" >>rotate.out
   echo "" >>rotate.out
   rotate.exe <rotate.dat>>rotate.out 
done
