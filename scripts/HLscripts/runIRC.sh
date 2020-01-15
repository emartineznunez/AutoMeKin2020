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
cd "$2"/IRC
name="$(sqlite3 inputs.db "select name from gaussian where id=$1")"
echo -e "$(sqlite3 inputs.db "select input from gaussian where id=$1")\n\n" > ${name}.dat
g09 <${name}.dat >${name}.log
rm ${name}.dat
