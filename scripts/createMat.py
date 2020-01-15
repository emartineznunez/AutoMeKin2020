#!/usr/bin/env python3.6
import networkx as nx
import numpy as np
import sys
from ase.io import read

#weig = 1 (integer numbers)
#weig = 2 (float numbers)
syst = str(sys.argv[1])
weig = int(sys.argv[2])

rmol = read(syst)
#Setting up initial stuff
natom  = len(rmol)
d      = rmol.get_all_distances(mic=False, vector=False)
symb   = rmol.get_chemical_symbols()
###
if len(sys.argv) == 4:
   na=int(sys.argv[3])
else:
   na=natom
##Three files are written here: Labels, ConnMat, ScalMat
outputfile2 = open('ConnMat', "w")
outputfile3 = open('ScalMat', "w")
###
###if there is one atom A has one element (0). Do it and exit
if natom == 1:
   outputfile2.write('0\n')
   exit()
###
sx1,sx2,d_ref,bondtype = np.genfromtxt("thdist",dtype="|U5",unpack=True)
dict_ref = {}
for i in range(len(d_ref)):
    dict_ref[sx1[i]+sx2[i]] = float(d_ref[i])
    if sx1[i] != sx2[i]:
        dict_ref[sx2[i]+sx1[i]] = float(d_ref[i])

if na < natom:
   sx1,sx2,d_ref,bondtype = np.genfromtxt("thdist_vdw",dtype="|U5",unpack=True)
   dict_ref_vdw = {}
   for i in range(len(d_ref)):
       dict_ref_vdw[sx1[i]+sx2[i]] = float(d_ref[i])
       if sx1[i] != sx2[i]:
           dict_ref_vdw[sx2[i]+sx1[i]] = float(d_ref[i])

G   = nx.Graph()
for i in range(natom):
    for j in range(i+1,natom):
        if j < na or i >= na:
           outputfile3.write(str ( dict_ref[ symb[i] + symb[j] ] )+' \n' )
           ratio = d[i][j] / 1.1 / dict_ref[ symb[i] + symb[j] ]
           diff  = d[i][j] - 1.1 * dict_ref[ symb[i] + symb[j] ]
        else:
           outputfile3.write(str ( dict_ref_vdw[ symb[i] + symb[j] ] )+' \n' )
           ratio = d[i][j] / 1.1 / dict_ref_vdw[ symb[i] + symb[j] ]
           diff  = d[i][j] - 1.1 * dict_ref_vdw[ symb[i] + symb[j] ]
        if weig == 1:
           weight = int(np.heaviside(-diff,0))
        elif weig == 2:
           weight = round( (1 - ratio**6) / (1 - ratio**12),5 )
        G.add_edge(i,j,weight=weight)

#Printing connectivity matrix
blank=' '
A=nx.adjacency_matrix(G)
Ar=A.A
for row in Ar:
    for item in row:
        outputfile2.write(str(item)+blank)
    outputfile2.write('\n')
