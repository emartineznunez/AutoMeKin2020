#!/usr/bin/python

import sys
import os
import shutil
import math
import subprocess

#Functions
def fileread(name):
        with open(name) as f:
                line=f.read().splitlines()
        return line

#
if os.path.isfile('amk.dat'):
        inputfile='amk.dat'
else:
        inputfile='amk.dat'

print "Path to inputfile:",inputfile

min_name={}

xyz=fileread(inputfile)
for i in range(len(xyz)):
        col=xyz[i].split()
        if len(col)>0 and col[0]=='species':
                nmolecule=len(col)-1
                min_name=col

#
print "\nNumber of molecules=",nmolecule
if nmolecule>4:
	print "right now the maximum number of molecules is 4"
	sys.exit(0)	
for i in range(1,nmolecule+1):
	print "Molecule",i,"=",min_name[i]

if not min_name[1]=="cat":
	print "The first molecule should be cat"
        sys.exit(0)

#Now we calculate the number of possible combinations of nmolecules 
subsystem=['cat']
pairs={}
trios={}
fours={}
i=1
#First the pairs
print "\nCombinations with all molecules"
print "Single"
print 1,subsystem[0]
print "Pairs"
p=0
for i in range(2,nmolecule+1):
 	p+=1
	pairs[p]='cat_'+min_name[i]
	print p,pairs[p]
subsystem=subsystem+pairs.values()
print "Trios"
t=0
for i in range(1,p+1):
	for j in range(i+2,nmolecule+1):
                t+=1
		trios[t]=pairs[i]+'_'+min_name[j]
		print t,trios[t]
subsystem=subsystem+trios.values()
print "Fours"
fo=1
fours[fo]=trios[1]+'_'+min_name[nmolecule]
print fo,fours[fo]
subsystem=subsystem+fours.values()
nsubsystem=1+p+t+fo
print "Number of Subsystems=",nsubsystem
subs=open('subs','w')
for i in range(nsubsystem):
	print i+1,subsystem[i]
	subs.write(str(subsystem[i])+'\n')
subs.close()
