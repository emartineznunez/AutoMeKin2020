#!/usr/bin/env python3 

import numpy as np
import math
from ase import Atoms
import io
import os
import sys
from mopacamk import MOPACamk
from ase.constraints import ExternalForce
from ase.io import read, write
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution,Stationary, ZeroRotation)
from ase import units
from numpy.random import standard_normal
import Langevin
import re
import Connectivity
from ase.optimize import BFGS

def runTrajectory(geom, T, Fric, totaltime, dt, adapLimit,window):
#r0 is a list with the initial distances
    size = len(geom.get_positions())
    r0 = np.zeros((size,size))
    for i in range(0,size):
        for j in range(i+1,size):
            r0[i][j]=geom.get_distance(i,j)
           
    timeStep = dt * units.fs
    numberOfSteps = int(totaltime/dt)
    BXDElower = 0
    BXDEupper = 0

    #get starting geometry for reference
    start = geom.copy()

## We set an external force of 4 eV/Ang on atoms 14 and 57 (python starts in 0)
## The - indicates an attractive force

#    c = ExternalForce(13,56,-4)
#    geom.set_constraint(c) 
##

    reactionCountDown = 0

    vel = geom.get_velocities()

    # Open files for saving data
    bofile = open("bond_order.txt", "w")
    # Initialise constraints to zero
    e_del_phi = 0
    j=0
    BXDon = True
    Call_criteria = True 
    Rxn = False
    Frag = False

    # Get current potential energy and check if BXD energy constraints are implemented
    ene = geom.get_potential_energy()
    e0=ene
    forces = geom.get_forces()

    BXDElower = ene - 10.0
    BXDEupper = ene + 5000

    maxEne = ene
    # Then set up reaction criteria or connectivity map
    con = Connectivity.BBFS(geom)
    mdInt = Langevin.Langevin(units.kB * T, Fric, forces, geom.get_velocities(), geom, timeStep)

    #Get forces
    forces = geom.get_forces()

    # Run MD trajectory for specified number of steps
    for i in range(0,numberOfSteps):


        # Get forces and energy from designated potential
        ene = geom.get_potential_energy()
   
        # Check for boundary hit on lower BXDE
        if BXDon and ene < BXDElower:
            hits += 1
            eBounded = True
            geom.set_positions(mdInt.oldPos)
            e_del_phi = geom.get_forces()
        else:
            eBounded = False
            hits = 0

#hits is a new variable to escape from trappings
        if hits >= 10:
            BXDElower = ene - 0.1
            print("Getting out of a trapping state")
            print("New BXDElower",(BXDElower - e0)*23.06)
#hits is a new variable to escape from trappings

        #check if we have passed upper boundary
        if ene > BXDEupper:
            j = 0
            BXDElower = BXDEupper - 0.1
            BXDEupper = BXDElower + 50000
            print("readjust boundaries") 
            e = (BXDElower - e0)*23.06
            print(e)
        if j < adapLimit and eBounded == False:
            j += 1
            if ene > maxEne:
                maxEne = ene
        elif j == adapLimit:
            j += 1
            BXDEupper = maxEne
        elif j > adapLimit:
            j += 1
            if (j/adapLimit)%4==0:
               diff = BXDEupper - BXDElower
               BXDEupper = BXDEupper - diff / 2
               print("new BXDEupper") 

        if eBounded is True:
            mdInt.constrain(e_del_phi)

        mdInt.mdStepPos(forces,timeStep,geom)
        forces = geom.get_forces()
        mdInt.mdStepVel(forces,timeStep,geom)

        #Print current positions to file
        idt = int(1/dt)
        if i % idt == 0:
            write("traj.xyz",geom.copy(),format="xyz",append=True)
            e = (ene - e0)*23.06
            emax=(BXDEupper - e0)*23.06
            emin=(BXDElower - e0)*23.06
            print(i*dt,j,e,emin,emax,BXDon,Rxn,reactionCountDown)
#read bond orders and write them to file bond_order.txt
            for item in geom.calc.get_bond_order():
                 bofile.write(item+' ')
            bofile.write('\n')
#end of reading bond orders
        con.update(geom)
        if Call_criteria is True:
            Rxn = con.criteriaMet
#Test to see if there is breakage
        for ii in range(0,size):
            for jj in range(ii+1,size):
                dist = geom.get_distance(ii,jj)
                if dist >= 5*r0[ii][jj]:
                    print("Fragmentation")
                    print(ii+1,jj+1,dist,r0[ii][jj])
                    Frag = True
        if Frag is True:
            break

        if Rxn is True:
            if reactionCountDown == 0:
                print("Possible Reaction")
                print(i*dt,con.ind1,con.ind2,con.ind3)
                reactionCountDown = window
                indicies = con.transitionIndices
            elif reactionCountDown > 1:
                reactionCountDown -= 1
                if reactionCountDown >= window - 20*idt:
                    if i % idt == 0:
                        print(i*dt,con.ind1,con.ind2,con.ind3)
                if reactionCountDown < window - 20*idt:
                    if Call_criteria is True:
                        print("Reaction confirmed")

                    Call_criteria = False
        else:
            BXDon = True
            reactionCountDown = 0

#remove traj.xyz if it exists
os.system('rm -rf traj.xyz')

inputfile = str(sys.argv[1])
xyzfile="opt_start.xyz"
temp,fric,totaltime,dt,adap,window,method = 1000,0.5,5000,0.25,100,500,'pm7'
inp = open(inputfile, "r")
for line in inp:
    if re.search("temp ", line):
        temp = float(line.split()[1])
    if re.search("Friction", line):
        fric = float(line.split()[1])
    if re.search("fs", line):
        totaltime = int(line.split()[1])
    if re.search("Dt", line):
        dt = float(line.split()[1])
    if re.search("AdaptiveLimit", line):
        adap = int(line.split()[1])
    if re.search("Window", line):
        window = int(line.split()[1])
        window = int(window/dt)
    if re.search("LowLevel ", line):
        method = str(line.split(' ',1)[1]).rstrip()
rmol = read(xyzfile)
rmol.calc = MOPACamk(label='bxde', method=method+' PRTXYZ THREADS=1')
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
print("    BXDE input parameters    ")
print("    Temperature(K) = ",temp)
print("    Friction       = ",fric)
print("    Totaltime(fs)  = ",totaltime)
print("    Delta_t(fs)    = ",dt)
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

MaxwellBoltzmannDistribution(rmol, temp * units.kB)
runTrajectory(rmol,temp,fric,totaltime,dt,adap,window)

