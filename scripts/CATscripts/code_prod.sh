file=$1
prodi=$2
gl=`awk '{
if($10=="'$prodi'" || $12=="'$prodi'") {
  if($13>$15)
     print $13
  else
     print $15
  exit}
}' $file `
echo $gl
