lanza=$1
qsub -l num_proc=1,s_rt=300:00:00,s_vmem=1G,h_fsize=2G,h_stack=256M $lanza
