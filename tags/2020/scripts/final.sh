#!/bin/bash
sharedir=${AMK}/share

source utils.sh
#On exit remove tmp files
tmp_files=(tmp* convergence.*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)

if [ -f amk.dat ];then
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
cwd=$PWD
###reading input file
read_input
###
if [ $rate -eq 0 ]; then
   postb="T$temperature"
   units="s"
elif [ $rate -eq 1 ]; then
   postb="E$energy"
   units="ps"
fi
##

final=FINAL_LL_${molecule}
####
rm -rf $final 
mkdir $final 
mdir=${final}/normal_modes
mkdir $mdir

##copy the products
rm -rf tmp_pe
cp ${tsdirll}/PRODs/prod.db $final
sqlite3 $final/prod.db "select energy,formula from prod" | awk '{for (i=1;i<=NF;i++) printf "%s",$i;printf "\n"}' | sed 's@|@ @g' >tmp_pf
###
##prods within 1 kcal/mol are considered the same
awk '{e[NR]=$1;f[NR]=$2}
END{i=1
while(i<=NR){
   np[i,0]=i
   j=i+1
   while(j<=NR){
     diff=e[i]-e[j]
     diff=sqrt(diff*diff)
     if(e[i]==0 && e[j]==0) diff=10
     if(diff<1 && f[i]==f[j]) {n[i]+=1;np[i,n[i]]=j;p[j]=1}
     j++
     }
   if(p[i]!=1 && n[i]>=1) {for(k=0;k<n[i];k++) printf np[i,k] " ";print np[i,n[i]]}
   i++
   }
}' tmp_pf > tmp_ep
##Make sure there is an empty line
echo "" >> tmp_ep
##remove the same products
for i in $(awk '{for (i=2; i<=NF; i++) print $i}' tmp_ep)
do
   sqlite3 $final/prod.db "delete from prod where id=$i"
done
cat tmp_ep > tmp_RXNet
cat tmp_ep > tmp_RXNet.cg
##Making final RXNet and RXNet.cg files
awk 'BEGIN{if('$rate'==0) fl="DG"; if('$rate'==1) fl="DE"}
NR==1{printf " TS #    %2s(kcal/mol)    -------Path info--------\n",fl}
NR>1{printf "%5.0f %12.3f       %4s %4s %4s %4s %4s\n",$2,$5,$7,$8,$9,$10,$11}' ${tsdirll}/KMC/RXNet >> tmp_RXNet
if [ -f ${tsdirll}/PRODs/PRlist_kmc ];then
   awk 'NR>1{for(i=1;i<NF;i++) if(i!=3) printf $i " ";print $NF}' ${tsdirll}/PRODs/PRlist_kmc >> tmp_RXNet 
fi

awk 'BEGIN{np=0;fl=0}
{
if(fl==0) {
   ++jwiso
   niso[jwiso]=NF
   for(i=1;i<=NF;++i) {n[jwiso,i]=$i;niso2[$i]=NF}
   }
}
{if(NF==0) fl=1}
{
if(NF>1 && fl==1) {lp=0;lp1=0;lp2=0
 if($1 !="PROD" ) {lp=1}
 for(i=1;i<=jwiso;++i) {
   for(j=2;j<=niso[i];++j) {
       if($3 =="PROD" &&  $4==n[i,j]) {lp1=1;pr1=n[i,1];++np;pr[np]=pr1}
       if($6 =="PROD" &&  $7==n[i,j]) {lp2=1;pr2=n[i,1];++np;pr[np]=pr2}
       if($1 =="PROD" &&  $2==n[i,j]) {lp=1}
     }
   }
 if($3 =="PROD" && lp1==0) {++np;pr[np]=$4}
 if($6 =="PROD" && lp2==0) {++np;pr[np]=$7}
 if($1=="TS") print $0
 if($1!="PROD" && $1!="TS" && lp1==0 && lp2==0) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,$4,$5,$6,$7
 if($1!="PROD" && $1!="TS" && lp1==1 && lp2==0) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,pr1,$5,$6,$7
 if($1!="PROD" && $1!="TS" && lp1==0 && lp2==1) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,$4,$5,$6,pr2
 if($1!="PROD" && $1!="TS" && lp1==1 && lp2==1) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,pr1,$5,$6,pr2
 if(lp==0) {lpp=0
    for(i=1;i<=np;i++) {if($2==pr[i]) lpp=1}
    if(lpp==1) print $0}
 }
} ' tmp_RXNet  > ${final}/RXNet


awk 'BEGIN{if('$rate'==0) fl="DG"; if('$rate'==1) fl="DE"}
NR==1{printf " TS #    %2s(kcal/mol)    -------Path info--------\n",fl}
NR>2{printf "%5.0f %12.3f       %4s %4s %4s %4s %4s\n",$2,$5,$7,$8,$9,$10,$11}' ${tsdirll}/KMC/RXNet.cg >> tmp_RXNet.cg
if [ -f ${tsdirll}/PRODs/PRlist_kmc ];then
   awk 'NR>1{for(i=1;i<NF;i++) if(i!=3) printf $i " ";print $NF}' ${tsdirll}/PRODs/PRlist_kmc >> tmp_RXNet.cg 
fi

awk 'BEGIN{np=0;fl=0}
{
if(fl==0) {
   ++jwiso
   niso[jwiso]=NF
   for(i=1;i<=NF;++i) {n[jwiso,i]=$i;niso2[$i]=NF}
   }
}
{if(NF==0) fl=1}
{
if(NF>1 && fl==1) {lp=0;lp1=0;lp2=0
 if($1 !="PROD" ) {lp=1}
 for(i=1;i<=jwiso;++i) {
   for(j=2;j<=niso[i];++j) {
       if($3 =="PROD" &&  $4==n[i,j]) {lp1=1;pr1=n[i,1];++np;pr[np]=pr1}
       if($6 =="PROD" &&  $7==n[i,j]) {lp2=1;pr2=n[i,1];++np;pr[np]=pr2}
       if($1 =="PROD" &&  $2==n[i,j]) {lp=1}
     }
   }
 if($3 =="PROD" && lp1==0) {++np;pr[np]=$4}
 if($6 =="PROD" && lp2==0) {++np;pr[np]=$7}
 if($1=="TS") print $0
 if($1!="PROD" && $1!="TS" && lp1==0 && lp2==0) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,$4,$5,$6,$7
 if($1!="PROD" && $1!="TS" && lp1==1 && lp2==0) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,pr1,$5,$6,$7
 if($1!="PROD" && $1!="TS" && lp1==0 && lp2==1) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,$4,$5,$6,pr2
 if($1!="PROD" && $1!="TS" && lp1==1 && lp2==1) printf "%5s %12s       %4s %4s %4s %4s %4s\n",$1,$2,$3,pr1,$5,$6,pr2
 if(lp==0) {lpp=0
    for(i=1;i<=np;i++) {if($2==pr[i]) lpp=1}
    if(lpp==1) print $0}
 }
} ' tmp_RXNet.cg  > ${final}/RXNet.cg

##Add CONN or DISCONN in the last column of RXNet.cg 
file=${final}/RXNet.cg
awk '{if($1=="TS") print $1,$2}' ${tsdirll}/KMC/RXNet_long.cg_groupedprods >tmp_conn_proc
awk 'NR==FNR{++nts;a[nts]=$2}
NR>FNR{
   if($3=="MIN") {
     ts=$1
     flag="   DISCONN"
     for(i=1;i<=nts;i++) if(ts==a[i]) flag="   CONN"
     print $0,flag
   }
  else
     print $0
}' tmp_conn_proc $file> tmp_RXNet.cg
mv tmp_RXNet.cg ${final}/RXNet.cg
###
##Add CONN or DISCONN in the last column of RXNet.cg 

##copy the minima 
cp ${tsdirll}/MINs/SORTED/mins.db $final
sqlite3 $final/mins.db "alter table mins rename to min"
mv $final/mins.db $final/min.db

##copy the TSs 
cp ${tsdirll}/TSs/SORTED/tss.db $final
sqlite3 $final/tss.db "alter table tss rename to ts"
mv $final/tss.db $final/ts.db


##Making final MINinfo file
awk 'BEGIN{if('$rate'==0) fl="DG"; if('$rate'==1) fl="DE"}
NR==1{printf "MIN #    %2s(kcal/mol)\n",fl}
NR>=1{printf "%5.0f %12.3f\n",$2,$4}' ${tsdirll}/MINs/SORTED/MINlist_sorted > ${final}/MINinfo
if [ -f ${tsdirll}/working/conf_isomer.out ]; then
   echo "Conformational isomers are listed in the same line:" >> ${final}/MINinfo
   cat ${tsdirll}/working/conf_isomer.out >> ${final}/MINinfo
fi

##Making final TSinfo file
awk 'BEGIN{if('$rate'==0) fl="DG"; if('$rate'==1) fl="DE"}
NR==1{printf "TS  #    %2s(kcal/mol)\n",fl}
NR>=1{printf "%5.0f %12.3f\n",$2,$4}' ${tsdirll}/TSs/SORTED/TSlist_sorted > ${final}/TSinfo
if [ -f ${tsdirll}/working/conf_isomer_ts.out ]; then
   echo "Conformational isomers are listed in the same line:" >> ${final}/TSinfo
   cat ${tsdirll}/working/conf_isomer_ts.out >> ${final}/TSinfo
fi

##Making molden files to visualize freqs of TSs
n=0
for file in $(awk '{print $3}' ${tsdirll}/TSs/SORTED/TSlist_sorted) 
do
    f="$(basename $file .rxyz)"
    ((n=n+1)) 
    number="$(printf %04d ${n%})"
    get_NM_mopac.sh ${tsdirll}/${f}.out  $mdir/TS$number
done

##Making molden files to visualize freqs of MINs 
n=0
for file in $(sed 's/_min/ min/g;s/_0//' ${tsdirll}/MINs/SORTED/MINlist_sorted | awk '{print $4}') 
do
    f="$(basename $file .rxyz)"
    ((n=n+1)) 
    number="$(printf %04d ${n%})"
    if [ "$f" == "min0" ]; then
       get_NM_mopac.sh ${tsdirll}/MINs/${f}.out  $mdir/MIN$number
    else
       get_NM_mopac.sh ${tsdirll}/IRC/${f}.out  $mdir/MIN$number
    fi
done



# kinetics file
kmcfile=${tsdirll}/KMC/kmc${postb}.out 
brafile=${tsdirll}/KMC/branching${postb}.out

if [ -f $kmcfile ]; then
   cat $brafile >$final/kinetics$postb
   echo "" >> $final/kinetics$postb
   echo "Population of each species as a function of time" >>$final/kinetics$postb
   echo "++++++++++++++++++++++++++++++++++++++++++++++++" >>$final/kinetics$postb
   awk '{line[NR]=$0};/Calculation number/{fdata=NR+1};/Population/{ldata=NR}
   END{for(i=fdata;i<ldata;i++) print line[i]}' $kmcfile >tmp_kmc
   npro=$(awk 'BEGIN{npro=0};NR>1{++npro};END{print npro}' $brafile)
   pro="$(awk 'NR>1{j=0;for(i=NF-'$npro'+1;i<=NF;i++) {++j;if($i>0) l[j]=1}};END{for(i=1;i<='$npro';i++) if(l[i]==1) printf "%8s ",i}' tmp_kmc)"
   min="$(awk 'NR>1{nmin=NF-'$npro'-1;j=0;for(i=2;i<=NF-'$npro';i++) {++j;if($i>0) l[j]=1}};END{for(i=1;i<=nmin;i++) if(l[i]==1) printf "%8s ",i}' tmp_kmc)"
   echo  $npro "$pro" > tmp_kmc0
   for i in $(echo "$pro")
   do
      tmp="$(awk 'NR>1{++npro};{if(npro=='$i') for(i=2;i<=NF;i++) print $i}' $brafile | sed 's@ + @+@' )"
      echo $tmp | sed 's@ + @+@g' >>tmp_kmc0
   done
   echo "$min" >> tmp_kmc0
   awk 'BEGIN{units="'$units'"}
   {if(FNR==NR) if(FNR==1) {ntpro=$1;for(i=2;i<=NF;i++) {++npro;ipro[npro]=$i}
   for(i=1;i<=npro;i++) {getline;pro[ipro[i]]=$0}
   getline;for(i=1;i<=NF;i++) {++nmin;imin[i]=$i}
   }
   }
   /Time/{n=FNR
   for(i=1;i<=nmin;i++) min[i]="MIN"imin[i]
   printf "  Time(%2s) ",units
   for(i=1;i<=nmin;i++) printf "%8s ",min[i]
   for(i=1;i<=npro;i++) {printf "%20s ",pro[ipro[i]]}
   print ""}
   {if(FNR>n && NR>FNR) {
   printf "%10s ",$1;for(i=1;i<=nmin;i++) printf "%8s ",$(imin[i]+1)
   for(i=1;i<=npro;i++) {printf "%20s ",$(NF-ntpro+ipro[i])}
   print ""
   }
   }' tmp_kmc0 tmp_kmc >> $final/kinetics$postb
   ###Making plot
   awk '{l[NR]=$0};/Time/{ldata=NR+1};END{for(i=ldata;i<=NR;i++) print l[i]}' ${final}/kinetics$postb > ${final}/pop_data_$postb
   xhigh=$(awk '{xhigh=$1};END{print xhigh}' ${final}/pop_data_$postb)
   yhigh=$(awk 'BEGIN{yhigh=0};NR==1{for(i=2;i<=NF;i++) if($i>yhigh) yhigh=$i;print yhigh;exit}' ${final}/pop_data_$postb)
   nplot=$(awk 'NR==1{print NF-1;exit}' ${final}/pop_data_$postb)
   sed 's@xhigh@'$xhigh'@;s@yhigh@'$yhigh'@' ${sharedir}/pop_template.gnu > ${final}/population${postb}.gnu
   xtics=$(echo $xhigh | awk '{printf "%3.2e",$1/4}')
   echo set xtics $xtics >> ${final}/population${postb}.gnu
   echo "set xlabel 'Time ($units)' font 'Times-Roman, 18' " >> ${final}/population${postb}.gnu
   echo set key noreverse top  >> ${final}/population${postb}.gnu
   echo "set ylabel 'Population' font 'Times-Roman, 18'" >> ${final}/population${postb}.gnu
   echo "set multiplot" >> ${final}/population${postb}.gnu
   for i in $(seq $nplot)
   do
     title=$(awk '/Time/{print $(NF-'$nplot'+'$i')}' ${final}/kinetics$postb)
     if [ $i -gt 1 ]; then pre=re; fi
     ((j=i+1))
     echo "${pre}plot 'pop_data_$postb' u 1:$j w l title '$title'"  >>${final}/population${postb}.gnu
   done
   echo pause -1  >> ${final}/population${postb}.gnu
   plot_relevant.sh
   if [ -f diagram.gnu ]; then
      mv diagram.gnu ${final}/Energy_profile.gnu
   fi
fi
nx.sh LL
rxn.py LL ${molecule} 

if [ -f ${tsdirll}/KMC/RXNet.relevant ]; then
###create RXNet.rel
   file=${final}/RXNet.cg
   awk '{if($1=="TS") print $1,$2}' ${tsdirll}/KMC/RXNet.relevant >tmp_rel
   awk 'NR==FNR{++nts;a[nts]=$2}
   NR>FNR{
      if($3=="MIN") {
        ts=$1
        p=0
        for(i=1;i<=nts;i++) if(ts==a[i]) p=1
        if(p==1) print $0
      }
     else
        print $0
   }' tmp_rel $file | sed 's@CONN@@g' > ${final}/RXNet.rel
fi
###Add report.pdf to the final folder
######################################################v
formula=$(FormulaMOL.sh ${molecule}.xyz | awk 'NR==1{print $0}')
track.py ${formula}
if [ -f "convergence.txt" ]; then
   cp convergence.txt ${final}
fi
######################################################^
cd ${final}
######################################################V
if [ -f RXNet.rel ]; then
   gnuplot <population${postb}.gnu>population${postb}.pdf
   gnuplot <Energy_profile.gnu>Energy_profile.pdf
fi

#####################################################^
rm -rf population${postb}.gnu Energy_profile.gnu
#####################################################v
if ! [ -z "$AMK_REPORT" ]
then
   enscript -q --margins=60::: --header='$n|%W|Page $% of $=' -p MINinfo.ps MINinfo  
   enscript -q --margins=60::: --header='$n|%W|Page $% of $=' -p TSinfo.ps TSinfo  
   enscript -q --margins=60::: --header='$n|%W|Page $% of $=' -p RXNet.ps RXNet  
   enscript -q --margins=60::: --header='$n|%W|Page $% of $=' -p RXNet.cg.ps RXNet.cg  
   enscript -q --margins=60::: --header='$n|%W|Page $% of $=' -p convergence.ps convergence.txt  
   ps2pdf convergence.ps convergence.pdf
   ps2pdf MINinfo.ps MINinfo.pdf
   ps2pdf TSinfo.ps TSinfo.pdf
   ps2pdf RXNet.ps RXNet.pdf
   ps2pdf RXNet.cg.ps RXNet.cg.pdf
   pdftk ${sharedir}/header.pdf convergence.pdf MINinfo.pdf TSinfo.pdf RXNet.pdf RXNet.cg.pdf graph_all.pdf graph_kin.pdf population${postb}.pdf cat output report_in.pdf
   cpdf -scale-to-fit a4portrait report_in.pdf -o report.pdf
fi
#####################################################^
rm -rf MINinfo.* RXNet.cg.* RXNet.pdf RXNet.ps RXNet.rel.* TSinfo.* report_in.pdf convergence.ps convergence.pdf
cd ..
