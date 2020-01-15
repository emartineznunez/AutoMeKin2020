d=1.8
natom=$1
natomA=$2
tc=$3
bd=$4
nm=$5
nmet_label=$6
noHLcalc=$7
charge=$8
mult=$9
echo "%chk=calc"$tc > $bd/diss$tc.dat
echo "%chk=calcop"$tc > $bd/fragtoop$tc.dat
sed 's/calc/opt=modredundant int(grid=ultrafine)/g' g09inp >> $bd/diss$tc.dat
sed 's/calc/opt=(calcall,noraman,cartesian,maxcycle=100)/g' g09inp >> $bd/fragtoop$tc.dat
put_frag_last.sh $nm >mingeom0
anchat=`awk 'BEGIN{dmin=10^10}
    {if(NR>'$natomA'){++natb
    x[natb]=$2;y[natb]=$3;z[natb]=$4}
    }
    {if(NR=='$nmet_label') {x[0]=$2;y[0]=$3;z[0]=$4}}
END{
i=1
while(i<=natb){
   dx=x[i]-x[0]
   dy=y[i]-y[0]
   dz=z[i]-z[0]
   d=sqrt(dx^2+dy^2+dz^2)
   if(d<dmin) {dmin=d;anch=i}
   i++
   }
print anch+'$natomA'
}' mingeom0`
sep_frag.sh $natomA $d >> $bd/diss$tc.dat
echo "$natom" > $bd/diss$tc.xyz
echo "" >> $bd/diss$tc.xyz
sep_frag.sh $natomA $d >> $bd/diss$tc.xyz
echo "" >> $bd/diss$tc.dat
echo $nmet_label $anchat "f" >> $bd/diss$tc.dat
echo "" >> $bd/diss$tc.dat
cat fragtoop >> $bd/fragtoop$tc.dat
echo "" >> $bd/fragtoop$tc.dat
if [ $noHLcalc -eq 2 ]; then
  level2=$10
  basis2=$11
  sed 's/chk=/chk=min0/g' sp_template | sed 's/level2/'$level2'\/'$basis2'/g' | sed 's/charge/'$charge'/g' | sed 's/mult/'$mult'/g' >> $bd/fragtoop$tc.dat
fi

