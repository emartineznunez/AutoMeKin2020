rm -f tswrk/geom0 tswrk/geom1 
file=$1
natom=$2
awk 'BEGIN{lok=0}
/FINAL HE/{lok=1}
/CARTESIAN COORDINATES/{
getline
if(lok==0) getline
if(lok==0) getline
i=1
while(i<='$natom'){
 getline
 if(lok==0 && i==1) print $3,$4,$5 > "tswrk/geom0"
 if(lok==0 && i>1) print $3,$4,$5 >> "tswrk/geom0"
 if(lok==1 && i==1) print $3,$4,$5 > "tswrk/geom1"
 if(lok==1 && i>1) print $3,$4,$5 >> "tswrk/geom1"
 i++
 }
}' $file
