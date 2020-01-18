#!/usr/bin/env python3.6
from networkx import Graph, adjacency_matrix
from numpy import genfromtxt, heaviside
from sys import argv
from ase.io import read

#weig = 1 (integer numbers)
#weig = 2 (float numbers)
#weig = 3 (both: float and integer)
syst = str(argv[1])
weig = int(argv[2])

rmol = read(syst)
#Setting up initial stuff
natom  = len(rmol)
d      = rmol.get_all_distances(mic=False, vector=False)
symb   = rmol.get_chemical_symbols()
###
if len(argv) == 4:
   na=int(argv[3])
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
sx1,sx2,d_ref,bondtype = genfromtxt("thdist",dtype="|U5",unpack=True)
dict_ref = {}
for i in range(len(d_ref)):
    dict_ref[sx1[i]+sx2[i]] = float(d_ref[i])
    if sx1[i] != sx2[i]:
        dict_ref[sx2[i]+sx1[i]] = float(d_ref[i])

if na < natom:
   sx1,sx2,d_ref,bondtype = genfromtxt("thdist_vdw",dtype="|U5",unpack=True)
   dict_ref_vdw = {}
   for i in range(len(d_ref)):
       dict_ref_vdw[sx1[i]+sx2[i]] = float(d_ref[i])
       if sx1[i] != sx2[i]:
           dict_ref_vdw[sx2[i]+sx1[i]] = float(d_ref[i])

G1   = Graph()
G2   = Graph()
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
        weight1 = int(heaviside(-diff,0))
        weight2 = round( (1 - ratio**6) / (1 - ratio**12),5 )
        G1.add_edge(i,j,weight=weight1)
        G2.add_edge(i,j,weight=weight2)

#Printing connectivity matrix
blank = ' '
A = adjacency_matrix(G1)
Ar1 = A.A
A = adjacency_matrix(G2)
Ar2 = A.A
if weig >=2:
   for row in Ar2:
       for item in row:
           outputfile2.write(str(item)+blank)
       outputfile2.write('\n')
if weig !=2:
   for row in Ar1:
       for item in row:
           outputfile2.write(str(item)+blank)
       outputfile2.write('\n')

