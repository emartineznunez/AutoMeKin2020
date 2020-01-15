#!/bin/bash
name="$(sqlite3 $2/IRC/inputs.db "select name from mopac where id=$1")"
mopac ${2}/IRC/minf_${name}.mop 2>/dev/null
mopac ${2}/IRC/minr_${name}.mop 2>/dev/null

