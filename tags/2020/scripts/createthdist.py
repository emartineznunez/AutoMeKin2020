#!/usr/bin/env python3
import sys
from ase.io import read

fexit = int(sys.argv[1])
mol   = str(sys.argv[2])
na    = int(sys.argv[3])
amk   = str(sys.argv[4])

if fexit ==1:
   exit()

rmol  = read(mol+'.xyz')
natom = len(rmol)
label = rmol.get_chemical_symbols()

inp = open(amk+'/share/list_of_cov_radii','r')
r_cov={}
for line in inp:
   i=line.split()[0]
   r_cov[i]=float(line.split()[1])

if na < natom: 
   inp = open(amk+'/share/list_of_vdw_radii','r')
   r_vdw={}
   for line in inp:
      i=line.split()[0]
      r_vdw[i]=float(line.split()[1])
      outputfile2 = open('thdist_vdw', "w")

line=[]
liner=[]
outputfile = open('thdist', "w")
for i1 in range(natom):
   for i2 in range(i1+1,natom):
         if i2 < na or i1 >= na:
            rth = str (round(1.2 * ( r_cov[label[i1]] + r_cov[label[i2]] ) ,3))
            flag='cov'
         else:
            rth = str (round(1.2 * ( r_vdw[label[i1]] + r_vdw[label[i2]] ) ,3))
            flag='vdw'
         exp = label[i1]+" "+label[i2]+" "+rth+" "+flag
         expr= label[i2]+" "+label[i1]+" "+rth+" "+flag
         if not exp in line and not exp in liner:
            if flag == 'cov':
               outputfile.write(exp+'\n') 
            else:
               outputfile2.write(exp+'\n') 
         line.append(exp)   
         liner.append(expr)   


