#!/bin/bash
sharedir=${AMK}/share
elements=${sharedir}/elements
file=$1
q=$2

awk 'BEGIN{p=0}
{if(NR == FNR) n[$1]=NR}
{if(NF==4) noe+=n[$1]
}
END{res=noe-'$q';reso2=res/2;if(reso2==int(reso2)) p=1;print p}' $elements $file
