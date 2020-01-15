file=$1
min=$2
awk '{if($2=='$min') print $4}' $file

