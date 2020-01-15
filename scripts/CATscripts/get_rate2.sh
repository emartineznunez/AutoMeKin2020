file="fit.log"
vol=$1
awk 'BEGIN{avog=6.022e23}
/b               =/{rate=$3/avog/'$vol'}
END{print rate}' $file
