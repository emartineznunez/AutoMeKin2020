dir=$1
name=$2
lanza=$3
tc=$4
echo "module load g09" > $lanza
echo "g09<$dir/$name.dat>$dir/$name.log" >> $lanza
echo "echo \"%chk=down$tc\" > $dir/down$tc.dat" >> $lanza
echo "sed 's/calc/irc=(downhill,maxpoints=500,stepsize=50,gradientonly) iop(1\/108=-1)/g' g09inp >> $dir/down$tc.dat" >> $lanza
echo "get_lastgeomg09.sh $dir/$name.log >> $dir/down$tc.dat" >> $lanza
echo "g09<$dir/down$tc.dat>$dir/down$tc.log" >> $lanza
qsub -l num_proc=4,s_rt=10:00:00,s_vmem=10G,h_fsize=20G,h_stack=256M $lanza
