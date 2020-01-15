#!/usr/bin/python

import sys
import os
import shutil
import math
import subprocess

#Functions
def system_call(command):
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read()


def fileread(name):
        with open(name) as f:
                line=f.read().splitlines()
        return line

#use os.path.isfile to check whether a file exists

cmd='create_kmc_template.sh'
os.system(cmd)

nsteps=50
j=0
#
dc0={}
dc0_max={}
sol={}
inputfile='catalysis_res/KMC/kmc_template.dat'
makekmc='makekmci.sh'
makekmci='catalysis_res/KMC/makekmcinput'
icon=fileread(inputfile)
for i in range(len(icon)):
       	col=icon[i].split()
       	if 'C0' in col[0]:
		j+=1
		dc0[j]=float(col[1])
		dc0_max[j]=float(col[2])
		if (len(col))==4:
			sol=float(col[3])
			dc0[j]=dc0[j]*sol
			dc0_max[j]=dc0_max[j]*sol

c0i=dc0.values()
c0_max=dc0_max.values()
dflag={}
dc0f={}

for ie in range(len(c0i)):
        dflag[ie]='C0_'+str(ie+1)
flag=dflag.values()

n=0
for ie in range(len(c0i)):
	print "Varying concentration of species ",ie+1
	for i in range(nsteps):
		n+=1
		delta=float(c0_max[ie]/nsteps)
		c0x=delta*(i+1)
		for iei in range(len(c0i)):
			if iei==ie:
				dc0f[iei]=c0x
			else:
				dc0f[iei]=c0i[iei]
		c0f=dc0f.values()
		makekmcinput=open(makekmci,'w')
		for jk in range(len(c0i)):
			makekmcinput.write('C0_'+str(jk+1)+'\t'+str(c0f[jk])+'\n')
		makekmcinput.close()
		outputfile='catalysis_res/KMC/kmc_files/kmc_'+str(ie)+'_'+str(i)+'.dat'
		cmd=makekmc+'\t'+makekmci+'\t'+inputfile+'\t'+outputfile
 		os.system(cmd)
 		lanza='catalysis_res/KMC/kmc_files/lanza_'+str(n)+'.sh'	
		with open(lanza,'w') as file:
			file.write('cd '+'catalysis_res/KMC/kmc_files\n')
			file.write('kmc_bim.exe<kmc_'+str(ie)+'_'+str(i)+'.dat>kmc_'+str(ie)+'_'+str(i)+'.out\n')
		cmd='queue_kmc.sh '+lanza
      	 	os.system(cmd)
		
