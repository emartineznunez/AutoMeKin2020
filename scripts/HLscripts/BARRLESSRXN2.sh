if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing"
   exit
fi

molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
charge=`awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' $inputfile `
mult=`awk 'BEGIN{mu=1};{if($1=="mult") mu=$2};END{print mu}' $inputfile `
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirHL_'$molecule'"}' $inputfile`

if [ ! -d "$tsdirhl/MINs/SCAN" ]; then
   echo "$tsdirhl/MINs/SCAN does not exist. It will be created"
   mkdir $tsdirhl/MINs/SCAN
else
   echo "$tsdirhl/MINs/SCAN already exists"
fi

n=0
ls $tsdirhl/MINs/BO/*_bo.log > dumi
set `awk '{print $0}' dumi`
for i
do 
  echo "Scan Coordnates" > scan_coor.out
  echo "Sprint all scan coords" > Sprint
  ((n=n+1))
  echo $i
  name2=`basename $i _bo.log`
  awk 'BEGIN{huge=1000000;read=0}
  /orientation:/{ getline
  getline
  getline
  getline
  i=1
  while(i<=huge){
    getline
    if(NF==1) break
    natom=i
    i++
    }
  }
  /Wiberg bond index matrix/{read=1}
  {
  {if($1=="Atom") nn=$2}
  if(read==1 && NF>2 && $1!="Atom" && $1!="----") {
      i=3
      while(i<=NF){
        a1=int($1)
        a2=i-2+(nn-1)
        if($i>0.5 && $i <1.5 && a1<a2 )  print a1,a2
        i++
        }
      }
  }
  /Wiberg bond index, Totals by atom/{exit}' $i > scan_coor
  echo "Labels" >labels
  awk '{if(NR>=3 && NF==4) print $0 }' $tsdirhl/MINs/SORTED/MIN$name2"_"*.rxyz > mingeom
  awk '{if(NR>=3 && NF==4) print $1 }' $tsdirhl/MINs/SORTED/MIN$name2"_"*.rxyz >> labels

  awk '{if(NR<=2) print $0 };{if(NR>=3 && NF==4) print $0 }' $tsdirhl/MINs/SORTED/MIN$name2"_"*.rxyz > $name2".xyz"

  natom=`wc -l mingeom| awk '{print $1}' `

  createMat.sh
  mv ConnMat ConnMat0
 
  set `awk '{print NR}' scan_coor ` 
  for j
  do
    ind=`awk 'NR=='$j',NR=='$j'{print $1}' scan_coor`
    jnd=`awk 'NR=='$j',NR=='$j'{print $2}' scan_coor`
    echo "1" $natom > sprint.dat
    createMat2.sh
    cat ConnMat >> sprint.dat
    createMatnbo.sh $ind $jnd
    cat ConnMat >> sprint.dat
    sprint2.exe <sprint.dat >sprint.out_$ind$jnd
    ok=`awk 'BEGIN{ll=0};/Results for the La/{getline; if($6<0.001 && $7<0.001) ll=1};END{print ll}' sprint.out_$ind$jnd `
    awk '/Results for the La/{getline; s = ""; for (i = 6; i <= NF; i++) s = s $i " "; print s }' sprint.out_$ind$jnd >> Sprint 
    echo "Coord" $j "indexes" $ind $jnd "OK" $ok
    if [ $ok -eq 1 ]; then
       echo "COORD" $ind $jnd >> scan_coor.out 
    else
       echo "This coordinate does not lead to 2 products"
    fi
  done
# check for redundant coordinates

  paste scan_coor.out Sprint > scosp
  awk '/COORD/{print $0}' scosp >scosp.out
  awk '/COORD/{for(j=6; j<=NF; j++) mat[j,NR]=$j
  i=1
  while(i<=NR-1){
    rms=0
    j=6
    while(j<=NF){
      diff=(mat[j,NR]-mat[j,i])
      rms+=(diff*diff)
      j++
      }
    if(rms==0) print NR,i
    i++
    }
  }' scosp.out >redcoor
  cat redcoor scosp.out  >dumi
  awk '{if(NF==2) ln[NR]=$1}
  /COORD/{tot=NR-1;++j
  i=1
  nop=0
  while(i<=tot){
    if(j==ln[i]) nop=1
    i++
    }
  if(nop==0) print $2,$3
  }' dumi> scan_coor.log
  cp scan_coor.log scan_coor$name2".log"

  casscf_scans.py $name2 

done
