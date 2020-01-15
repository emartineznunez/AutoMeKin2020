#!/bin/bash
#if [ -n "$LUSTRE_SCRATCH" ]
#then
#  mkdir $LUSTRE_SCRATCH/$2
#  export LUSTRE_SCRATCH=$LUSTRE_SCRATCH/$2
#fi
#if [ -n "$TMPDIR" ]
#then
#  mkdir $TMPDIR/$2
#  export TMPDIR="$TMPDIR/$2"
#fi

cd "$2"
name="$(sqlite3 inputs.db "select name from gaussian where id=$1")"
inp="$(sqlite3 inputs.db "select input from gaussian where id=$1")"
echo -e "$inp\n\n" > ${name}.dat
if [ "$inp" != "salir" ]; then
   g09 <${name}.dat >${name}.log
   t=$(awk 'BEGIN{t=0};/Error termination via/{t=1};END{print t}' ${name}.log)
   if [ $t -eq 1 ]; then
#Constructing new input file
     echo -e "$inp\n\n" | sed 's/calcall,noraman/cartesian,maxcycle=100,calcall,noraman/g' > ${name}.dat
     g09 <${name}.dat> ${name}.log
   fi
fi
#rm ${name}.dat
