#!/bin/bash
# Run this script to test the program for the formic acid example
# You must have loaded amk module before:
# module load amk/2020

# run_test.sh --prefix=path_to_program
# path_to_program: is the path to the installation directory ($HOME/amk-2020 by default)
# ntasks: is number of parallel tasks
# niter: is the number of iterations
cwd=$PWD

path_to_program=$HOME/amk-2020/examples
tests0=(assoc rdiels_bias diels_bias FA_biasH2 FA_biasH2O FA_bxde FA_singletraj FA FAthermo FA_programopt vdW)

if [ $# -ge 1 ]; then
   args=${@}
   path_to_program=$(echo $args | sed  's/=/ /g' | awk 'BEGIN{f=0};{for(i=1;i<=NF;++i) if($i=="--prefix") f=i}
   END{if(f>0)
         print $(f+1)"/examples"
       else
         print "'$path_to_program'" }' | sed 's/\/\//\//')
   tests="$(echo $args | sed  's/=/ /g;s/,/ /g' | awk 'BEGIN{f=0};{for(i=1;i<=NF;++i) if($i=="--tests") f=i}
   END{if(f>0)
         for(i=f+1;i<=NF;i++) {if($i~/--/) exit;printf "%s ",$i}
       else
         print "all"
        }' )"
   if [[ "$tests" == *"all"* ]]; then
      tests=${tests0[@]}
   fi
else
   tests=${tests0[@]}
fi

if [ ! -f ${path_to_program}/FA.dat ]; then
   echo "The path to amk-2020 is not correct"
   exit
fi

ptp=$path_to_program
files=($ptp/*.dat $ptp/*.xyz)
ntasks=10
niter=2
runningtasks=$ntasks

for i in $(echo ${tests[@]})
do
   echo Running $i test. Directory: $i 
   if [ ! -f ${ptp}/${i}.dat  ];then
      echo "this test does not exist" 
      exit
   fi
   rm -rf $i && mkdir $i
   cd $i   
   cp ${files[@]} .
   if [ "$i" == "FA" ] || [ "$i" == "FA_bxde" ] || [ "$i" == "FAthermo" ] || [ "$i" == "FA_programopt" ] || [ "$i" == "vdW" ] || [ "$i" == "diels_bias" ]; then
      echo "LL calculations"
      time llcalcs.sh ${i}.dat $ntasks $niter $runningtasks > llcalcs.log 
   else
      time amk.sh ${i}.dat > ${i}.log
   fi
   if [ "$i" == "FA_programopt" ]; then
      echo "HL calculations"
      time hlcalcs.sh ${i}.dat $runningtasks > hlcalcs.log 
   fi
   if [ "$i" == "FA" ]; then
      cp -r FINAL_LL_FA FINAL_LL_FA_coarse_grained
      tsn=$(awk 'NR==2{print $1}' FINAL_LL_FA_coarse_grained/RXNet.cg)
      echo "Removing ts number $tsn"
      remove_ts.sh $tsn  > removets.log
      mv FINAL_LL_FA FINAL_LL_FA_coarse_grained_TS${tsn}_removed
      echo "Considering all states"
      rxn_network.sh allstates > allstates.log
      kmc.sh >> allstates.log
      final.sh >> allstates.log
      mv FINAL_LL_FA FINAL_LL_FA_all_states
      echo "One-level HL calculations"
      time hlcalcs.sh ${i}.dat $runningtasks > hlcalcs_onelevel.log 
      mv FINAL_HL_FA FINAL_HL_FA_onelevel
      cp $ptp/FA_2level.dat FA.dat 
      mv tsdirHL_FA tsdirHL_FA_onelevel
      echo "Two-level HL calculations"
      time hlcalcs.sh ${i}.dat $runningtasks > hlcalcs_twolevel.log 
      mv FINAL_HL_FA FINAL_HL_FA_twolevel
      mv tsdirHL_FA tsdirHL_FA_twolevel
   fi
   if [ "$i" == "FAthermo" ]; then
      mv FINAL_LL_FA FINAL_LL_FA_T300
      htemp=5000 
      kinetics.sh ${htemp} ll > htemp.log
      mv FINAL_LL_FA FINAL_LL_FA_T${htemp}
   fi
   if [ "$i" == "assoc" ] || [ "$i" == "vdW" ]; then
      cp Bz-N2.xyz Bz-N2.xyz_tmp
   fi
   if [ "$i" == "FA_singletraj" ]; then 
      tsll_view.sh > tslist.log
   fi
   rm -rf *.dat *.xyz
   if [ "$i" == "assoc" ] || [ "$i" == "vdW" ]; then
      mv Bz-N2.xyz_tmp Bz-N2.xyz
   fi
   echo ""
   echo $i test Done
   echo ""
   cd $cwd
done

echo "###########################"
echo "###Summary of all tests:###"
echo "###########################"
echo
echo "assoc test--> association complex geometry: Bz-N2.xyz"
echo
echo "rdiels test--> Stochastic. Possible TSs for reaction in folder tsdirLL_rdiels"
echo
echo "diels test--> Stochastic. Possible TSs for reaction in folder FINAL_LL_diels_bias"
echo
echo "FA_biasH2 test--> Stochastic. Possible TSs for reaction in folder tsdirLL_FA"
echo
echo "FA_biasH2O test--> Stochastic. Possible TSs for reaction in folder tsdirLL_FA"
echo
echo "FA_bxde test--> Stochastic. Possible TSs for reaction in folder FINAL_LL_FA"
echo
echo "FA_singletraj test--> 2 TSs for reaction in folder tsdirLL_FA. Results of tsll_view.sh in tslist.log"
echo
echo "FA test--> Stochastic. Check results in FINAL_LL_FA_coarse grained, FINAL_LL_FA_coarse_grained_TS1_removed, FINAL_LL_FA_all_states, FINAL_HL_FA_onelevel and FINAL_HL_FA_twolevel"
echo
echo "FAthermo test--> Stochastic. Check results in FINAL_LL_FA_T300 and FINAL_LL_FA_T5000"
echo
echo "FA_programopt test--> Stochastic. Check results in FINAL_LL_FA_programopt and FINAL_HL_FA_programopt"
echo
echo "vdW test--> Stochastic. Check results in FINAL_LL_Bz-N2"

