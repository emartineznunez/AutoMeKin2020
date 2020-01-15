#!/bin/bash
sharedir=${AMK}/share
source utils.sh
#On exit remove tmp files
tmp_files=(tmp*)
trap 'err_report $LINENO' ERR
trap cleanup2 EXIT INT

exe=$(basename $0)
cwd=$PWD
inputfile=amk.dat
##reading HL stuff
read_input
##

sed 's/ts,noeigentest,//;s/tkmc/'$temperature'/;s@level1@'$level1'@;s/charge/'$charge'/;s/mult/'$mult'/;s@iop@'$iop'@' $sharedir/hl_input_template > tmp_min_opt
if [ -z $iop ]; then
   sed 's@level1@'$level1' sp@;s/charge/'$charge'/;s/mult/'$mult'/;/opt/,/temp/d' $sharedir/hl_input_template > tmp_min_sp_prod
else
   sed 's/opt=(ts,noeigentest,calcall,noraman)//;s@level1@'$level1' sp@;s/charge/'$charge'/;s/mult/'$mult'/;/temp/d;s@iop@'$iop'@' $sharedir/hl_input_template > tmp_min_sp_prod
fi
i=$1
frm1=$(awk 'BEGIN{npt=0};/Pt /{npt=$2};END{print npt}' $tsdirhl/IRC/ircf_$i.log )
frm2=$(awk 'BEGIN{npt=0};/Pt /{npt=$2};END{print npt}' $tsdirhl/IRC/ircr_$i.log )
if [ $frm1 -le 2 ] || [ $frm2 -le 2 ]; then
   mode=1
else
   mode=2
fi

chkf="$(echo "%chk=ircf_"$i)"
chkr="$(echo "%chk=ircr_"$i)"

if [ $mode -eq 1 ]; then
   get_NM_g09.sh $tsdirhl/$i.log  1 > tmp_geomf_$i
   get_NM_g09.sh $tsdirhl/$i.log -1 > tmp_geomr_$i
else
   get_geom_irc_g09.sh $tsdirhl/IRC/ircf_$i.log > tmp_geomf_$i
   get_geom_irc_g09.sh $tsdirhl/IRC/ircr_$i.log > tmp_geomr_$i
fi
nf=$(nfrag.sh tmp_geomf_$i ${nfrag_th})
if [ $nf -eq 1 ]; then
   calf="$(cat tmp_min_opt)"
else
   calf="$(cat tmp_min_sp_prod)"
fi

nf=$(nfrag.sh tmp_geomr_$i ${nfrag_th})
if [ $nf -eq 1 ]; then
   calr="$(cat tmp_min_opt)"
else
   calr="$(cat tmp_min_sp_prod)"
fi

geof="$(cat tmp_geomf_$i)"
geor="$(cat tmp_geomr_$i)"
ig09f="$(echo -e "$chkf"'\n'"$calf"'\n'"$geof")"
ig09r="$(echo -e "$chkr"'\n'"$calr"'\n'"$geor")"
if [ $noHLcalc -eq 2 ]; then
   spcf="$(sed 's/chk=/chk=ircf_'$i'/;s@level2@'$level2'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/sp_template)"
   spcr="$(sed 's/chk=/chk=ircr_'$i'/;s@level2@'$level2'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/sp_template)"
   ig09f="$(echo -e "$chkf"'\n'"$calf"'\n'"$geof"'\n\n'"$spcf")"
   ig09r="$(echo -e "$chkr"'\n'"$calr"'\n'"$geor"'\n\n'"$spcr")"
fi
echo -e "insert or ignore into gaussian values (NULL,'minf_$i','$ig09f');\n.quit" | sqlite3 ${tsdirhl}/IRC/inputs.db
echo -e "insert or ignore into gaussian values (NULL,'minr_$i','$ig09r');\n.quit" | sqlite3 ${tsdirhl}/IRC/inputs.db
###
