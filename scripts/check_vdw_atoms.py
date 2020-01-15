#!/usr/bin/env python3
import sys
from ase.io import read

mol   = str(sys.argv[1])
amk   = str(sys.argv[2])

rmol  = read(mol+'.xyz')
label = rmol.get_chemical_symbols()

inp = open(amk+'/share/list_of_vdw_radii','r')
ref_atoms = []
for line in inp:
   ref_atoms.append(line.split()[0])

ok = 1
for lab in label:
   if lab not in ref_atoms:
      ok = 0
print(ok)
