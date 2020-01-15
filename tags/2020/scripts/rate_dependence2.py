#!/usr/bin/python

import sys
import os
import shutil
import math
import subprocess
import numpy as np

#Functions
def system_call(command):
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read()


def fileread(name):
        with open(name) as f:
                line=f.read().splitlines()
        return line

nsteps=50
j=0
#
dc0_max={}
sol={}
originpf='amk.dat'
imt=fileread(originpf)
for i in range(len(imt)):
	col=imt[i].split()
	if len(col)>1 and col[0]=='Maxtime':
		maxtime=col[1]
	if len(col)>1 and col[0]=='Product':
		prod=col[1]
code_prod=int(system_call('code_prod.sh '+'catalysis_res/KMC/RXNet_catalysis_with_barrierless '+ prod))
inputfile='catalysis_res/KMC/kmc_template.dat'
dirkmcs='catalysis_res/KMC/kmc_files'
icon=fileread(inputfile)
for i in range(len(icon)):
        col=icon[i].split()
	if i==2:
		vol=col[0]
        if 'C0' in col[0]:
                j+=1
                dc0_max[j]=float(col[2])
                if (len(col))==4:
                        sol=float(col[3])
                        dc0_max[j]=dc0_max[j]*sol
c0_max=dc0_max.values()
cmd='cp '+'${AMK}/share/CAT/param .'
os.system(cmd)
for ie in range(len(c0_max)):
	rate_ie=open('catalysis_res/KMC/kmc_files/rate_'+str(ie),'w')
	rate_ie_stat=open('catalysis_res/KMC/kmc_files/rate_'+str(ie)+'_stat','w')
	for i in range(nsteps):
		delta=float(c0_max[ie]/nsteps)
		c0_i=delta*(i+1)
		file=dirkmcs+'/kmc_'+str(ie)+'_'+str(i)+'.out'
		rate={}
		if os.path.isfile(file):
			for j in range(1,4):
				cmd='get_rate.sh '+file+' '+str(j)+' '+str(maxtime)+' '+str(code_prod)
       				os.system(cmd)
				if not os.path.isfile('k'):
					print "Calc. ",j," of molecule ",ie," step ",i," has not finished yet"
				else:
					cmd='gnuplot '+'${AMK}/share/CAT/fit.gp 2>error.log'
       					os.system(cmd)
					rate[j]=float(system_call('get_rate2.sh '+vol))
					cmd='rm error.log fit.log'
       					os.system(cmd)
			media=np.mean(rate.values())
			std=np.std(rate.values())
			print i,c0_i,media,std
			rate_ie.write(str(c0_i)+' '+str(media)+' '+str(std)+'\n' )
			rate_ie_stat.write(str(rate.values())+'\n' )
		else:
			print "File ",file," does not exist"
	rate_ie.close()
	rate_ie_stat.close()
cmd='rm param'
os.system(cmd)
