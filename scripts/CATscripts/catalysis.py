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

#
if not os.path.exists('catalysis_res'):
	os.makedirs('catalysis_res')
	print "catalysis_res does not exist and it'll be created"

if not os.path.exists('catalysis_res/KMC'):
	os.makedirs('catalysis_res/KMC')
	print "catalysis_res/KMC does not exist and it'll be created"
#else:
#	shutil.rmtree('catalysis_res/KMC')
#	os.makedirs('catalysis_res/KMC')
#	print "catalysis_res/KMC already exists, but it'll be removed and created again"

#use os.path.isfile to check whether a file exists	

if os.path.isfile('amk.dat'):
        inputfile='amk.dat'
else:
        print "amk.dat is missing"
        sys.exit(0)

print "Path to inputfile:",inputfile

min_name={}
gap={}
min_form={}
min_form_subs={}
natom={}
energy={}
zpe={}
gcorr={}
TKMC='298'
noHLcalc=0

xyz=fileread(inputfile)
for i in range(len(xyz)):
        col=xyz[i].split()
        if len(col)>0 and col[0]=='species':
                nmolecule=len(col)-1
                min_name=col
        if len(col)>0 and col[0]=='TKMC':
                TKMC=col[1]
        if len(col)>0 and col[0]=='eta':
                eta=col[1]
        if len(col)>0 and col[0]=='HLcalc':
                noHLcalc=col[1]
                if int(noHLcalc)==0:
                        print "The level of theory should be specified"
                        sys.exit(0)
                elif int(noHLcalc)==1:
                        HLcalc1=col[2]
                        HLcalc=HLcalc1
                        print "\nHL using one level only for energies and frequencies ",HLcalc1
                elif int(noHLcalc)==2:
                        HLcalc1=col[2]
                        HLcalc2=col[3]
                        HLcalc=HLcalc2
                        print "\nHL using two levels ",HLcalc1," and ",HLcalc2
        if len(col)>0 and col[0]=='level1':
                        level1=col[1]
        if len(col)>0 and col[0]=='level2':
                        level2=col[1]

#Diffusion controlled rate constant calculation
diffbim=open('catalysis_res/diffbim.dat','w')
diffbim.write(str(eta)+' '+str(TKMC))
diffbim.close()
kdiff=float(system_call('diff_bim.exe<catalysis_res/diffbim.dat'))

#

if int(noHLcalc) == 1:
 	print level1
elif int(noHLcalc)==2:
	print level1,"//",level2
print "\nNumber of species=",nmolecule
if nmolecule>4:
	print "right now the maximum number of species is 4"
	sys.exit(0)	
for i in range(1,nmolecule+1):
	print "Molecule",i,"=",min_name[i]
print "Temperature of the KMC calcs=",TKMC
print "Diffusion-controlled rate constant=",kdiff," M-1 s-1"

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
print "Number of Subsystems=",nsubsystem,"\n"
#We start off optimizing the minima
for i in range(1,nmolecule+1):
	name=min_name[i]
	if os.path.isfile('catalysis_res/'+name+'.log'):
 		print name,"alredy optimized"
	else:
		print name," not optimized"
		fg09=open('catalysis_res/'+name+'.dat','w')
 		chk='%chk='+name
		theory='#p '+level1+' opt=(calcall,noraman) temperature='+TKMC
##construct g09 input file
		fg09.write(chk+'\n')
		fg09.write(theory+'\n')
		fg09.write('\nmin opt')
		fg09.write('\n\n')
		fg09.write('0 1\n')
		fcoor = name+'.xyz'
		xyz=fileread(fcoor)
		for ii in range(2,len(xyz)):
			fg09.write(xyz[ii]+'\n')
		fg09.write('\n')
		if int(noHLcalc) ==2:
			fg09.write("--Link1--\n")
			fg09.write(chk+'\n')
		        theory='#p '+level2+' geom=check'
		        fg09.write(theory+'\n')
# set up cesga script
		cmd='queueHL_sub.sh catalysis_res'+'\t'+ name +' lanza_'+name+'.sh'
                os.system(cmd)
		sys.exit(0)
cmd='sed "s/\'/ /g" atsymb >atsdum1'
os.system(cmd)
cmd='sed "s/,/ /g" atsdum1 > atsdum2'
os.system(cmd)
cmd='atsdum2.sh'
os.system(cmd)

cmd='rm -f catalysis_res/geome_* catalysis_res/frequ_* catalysis_res/min0.rxyz'
os.system(cmd)
	
for i in range(1,nmolecule+1):
	cmd='get_energy'+HLcalc+'.sh catalysis_res/'+min_name[i]+'.log '+noHLcalc
        energy[i]=float(system_call(cmd))
        line=fileread('catalysis_res/'+min_name[i]+'.log')
	with open('catalysis_res/'+min_name[i]+'.log','r') as f:
		for line in f:
			if 'Zero-point vibrational energy' in line:
				zpe[i]=float(f.next().split()[0])
			if 'Thermal correction to Gibbs Free Energy' in line:
				gcorr[i]=float(line.split()[6])
	cmd='log_to_rxyz.sh '+min_name[i]+'\t'+str(energy[i])+'\t'+str(zpe[i])+'\t'+str(gcorr[i])+'\t'+str(i)
	os.system(cmd)
	xyz=fileread(min_name[i]+'.xyz')
	natom[i]=int((xyz[0]))
entot=sum(energy.values())
zptot=sum(zpe.values())
gitot=sum(gcorr.values())
natot=sum(natom.values())


fmin0=open('catalysis_res/min0.rxyz','w')
fmin0.write(str(natot)+'\n')
fmin0.write('E= '+str(entot)+' ZPE= '+str(zptot)+' Gcorr= '+str(gitot)+'\n')
fmin0.close()
for i in range(1,nmolecule+1):	
	cmd=('cat catalysis_res/geome_'+str(i)+' >> catalysis_res/min0.rxyz')
	os.system(cmd)
for i in range(1,nmolecule+1):	
	cmd=('cat catalysis_res/frequ_'+str(i)+' >> catalysis_res/min0.rxyz')
	os.system(cmd)

cmd='rm -f catalysis_res/geome_* catalysis_res/frequ_*'
os.system(cmd)
# e0 is the energy of the separated nmolecule fragments. This is the zero of energy
e0=entot+gitot

#

enadf=open('catalysis_res/KMC/Energies','w')
rxnet=open('catalysis_res/KMC/RXNet_catalysis','w')
rates=open('catalysis_res/KMC/Rates','w')
correl=open('catalysis_res/KMC/Correl','w')
rxnet.write('Reference system:\n')
rxnet.write('Klabel    Name   Formula\n')
for i in range(2,nmolecule+1):
	min_form[i]=system_call('FormulaMOL.sh '+min_name[i]+'.xyz')
	rxnet.write(' %2s %10s %10s' % (str(i-1),min_name[i],min_form[i]))
min_form_list=min_form.values()
for i in range(len(min_form_list)):
	min_form_list[i]=min_form_list[i].strip()
##
rxnet.write('Reference energy: '+str(e0)+' Hartree\n')
enadf.write('Reference energy: '+str(e0)+' Hartree\n')
rxnet.write('Klabel==label employed in the KMC simulations\n')
rxnet.write('Gibbs Free Energy (G) in kcal/mol\n')
rxnet.write('Diffusion-controlled rate constant= '+str(kdiff)+' M-1 s-1\n')
##
## cat is now third
##
lastmm=nmolecule-1

for i in range(nsubsystem):
	print 'Subsystem ',i+1,' is: ',subsystem[i]

for i in range(nsubsystem):
	gap[i]=lastmm
	rates.write('Subsystem: '+str(i+1)+'\t'+subsystem[i]+'\n')
	correl.write('Subsystem: '+str(i+1)+'\t'+subsystem[i]+'\n')
	tsdirll='tsdirLL_'+subsystem[i]
	tsdirhl='tsdirHL_'+subsystem[i]
	if not os.path.exists(tsdirhl):
		print tsdirhl,'does not exist\n'
		print 'Complete all calculations before proceeding'
		sys.exit(0)
	min_form_subs[i]=system_call('FormulaMOL.sh '+tsdirhl+'/MINs/min0.rxyz')
	min_form_subs_list=min_form_subs.values()
	for ii in range(len(min_form_subs_list)):
		min_form_subs_list[ii]=min_form_subs_list[ii].strip()
	print '\nSubsystem',i+1,subsystem[i]
	rxnet.write('\nSubsystem '+str(i+1)+' '+subsystem[i]+' '+min_form_subs[i]+'\n')
	rxnet.write('  TSlabel        G_TS          MINlabel Klabel       G_MIN              PRlabel              Klabel        G_PRO  DEG_PATH\n')
	print 'tsdirll= ',tsdirll
	k=0
	edum={}
 	for j in range(1,nmolecule+1):
		if not min_name[j] in subsystem[i]:
			ee=fileread('catalysis_res/'+min_name[j]+'.rxyz') 	
        		col=ee[1].split()
			k+=1
			edum[k]=float(col[1])+float(col[5])
	if k>0:
		enadd=sum(edum.values())
	else:
		enadd=0
	enadf.write('Subsystem '+str(i+1)+' '+subsystem[i]+' Enadd= '+str(enadd)+'\n')
	file=tsdirhl+'/MINs/min0.rxyz'
	if os.path.isfile(file):
		ee=fileread(file) 	
        	col=ee[1].split()
		enall=float(col[1])+float(col[5])+float(enadd)
	else:
		print "You should first optimize the global minimum of this subsystem"
		sys.exit(0)
##ediff is the energy that has to be added to every single min and ts of subsystem[i]
	ediff=(enall-e0)*627.51
	rxnetfile=fileread(tsdirhl+'/KMC/RXNet_long.cg_groupedprods')
	ratefile=fileread(tsdirhl+'/KMC/TST/rate'+TKMC+'.out')
	minn={}
	kk=0
	for n in range(2,len(rxnetfile)):
		coll=rxnetfile[n].split()
		kk+=1
		minn[kk]=int(coll[7])
 		if coll[9] =='MIN':
			kk+=1
			minn[kk]=int(coll[10])
	lastmin=max(minn.values())+gap[i]
	prfile=fileread(tsdirhl+'/PRODs/PRlist_kmc')
	klabel={}
	r1={}
	ts={}
	p1={}
	p2={}
	kc=0
	for line in range(len(rxnetfile)):
		col=rxnetfile[line].split()
		if line >=2:
			file=tsdirhl+'/MINs/SORTED/MINlist_sorted'
			cmd='get_e+g.sh '+file+' '+col[7]
			emin1=float(system_call(cmd))
			klabel[int(col[7])]=int(col[7])+gap[i]
			if klabel[int(col[7])] > lastmm:
				lastmm=klabel[int(col[7])]
			if col[9]=='MIN':
			        cmd='get_e+g.sh '+file+' '+col[10]
			        emin2=float(system_call(cmd))
				klabel[int(col[10])]=int(col[10])+gap[i]
				if klabel[int(col[10])] >lastmm:
					lastmm=klabel[int(col[10])]
				rxnet.write('TS%5s_%s    %8.3f ==>    MIN%5s_%s  %4s      %7.3f    <-->   MIN%5s_%s               %4s      %8.3f\n' % (col[1],i+1,float(col[4])+ediff,col[7],i+1,klabel[int(col[7])],emin1+ediff,col[10],i+1,klabel[int(col[10])],emin2+ediff))
				kc+=1
				ts[kc]=int(col[1])
				r1[kc]=klabel[int(col[7])]				
				p1[kc]=klabel[int(col[10])]
				p2[kc]=0
				kc+=1
				ts[kc]=int(col[1])
				r1[kc]=klabel[int(col[10])]
				p1[kc]=klabel[int(col[7])]				
				p2[kc]=0
			if col[9]=='PROD':
				deg=col[13]
				code=col[2]
				for line2 in range(len(prfile)):
					col2=prfile[line2].split()
					if code in col2[2]:
						prcode=col2[1]
						prname=fileread(tsdirhl+'/PRODs/CALC/PR'+col2[1]+'_'+col2[2])
                        			cole=prname[1].split()	
                        			eprod=float(cole[1])+float(cole[5])+float(enadd)
						epr=(eprod-e0)*627.51
						del col2[0:3]
						matches=list(set(col2)&set(min_form_list))
						matches2=list(set(col2)&set(min_form_subs_list))
						for ik in range(len(matches)):
							klabel[int(col[10]),ik]=min_form_list.index(matches[ik])+1
						if len(matches)==0:
							lastmin+=1
							klabel[int(col[10]),0]=lastmin
							lastmin+=1
							klabel[int(col[10]),1]=lastmin	
						elif len(matches)==1:
							lastmin+=1
							klabel[int(col[10]),1]=lastmin
						pr=' '.join(col2)
	
						if len(matches2)==1:
							indxs=min_form_subs_list.index(matches2[0])
							matchmin=open('catalysis_res/matchmin.dat','w')
							matchmin.write('Where= '+subsystem[indxs]+'\n')
							matchmin.write('From= '+subsystem[i]+'\n')
							matchmin.write(min_form_subs[indxs]+' '+str(prcode)+'\n')
							matchmin.write(HLcalc+' '+noHLcalc+'\n')
							matchmin.close()
						        minp=system_call('matchmin.sh catalysis_res/matchmin.dat')
							os.system('rm -f catalysis_res/matchmin.dat')
							if int(minp)>0:
         							klabel[int(col[10]),1]=int(minp)+gap[indxs]
								lastmin-=1
						if klabel[int(col[10]),1] > lastmm:
							lastmm=klabel[int(col[10]),1] 
						deltag=float(col[4])+ediff-epr
						tstbim=open('catalysis_res/tstbim.dat','w')
						tstbim.write(str(deltag)+' '+str(TKMC)+' '+str(deg))
						tstbim.close()
						krev=float(system_call('tst_bim.exe<catalysis_res/tstbim.dat'))
						if krev>kdiff:
							krev=kdiff
						os.system('rm -f catalysis_res/tstbim.dat')
						rates.write('%30s %4s %4s %4s    0\n' % (krev,klabel[int(col[10]),0],klabel[int(col[10]),1],klabel[int(col[7])]))
						correl.write('%30s %4s %4s %4s    0 TS%5s_%s\n' % (krev,klabel[int(col[10]),0],klabel[int(col[10]),1],klabel[int(col[7])],col[1],i+1))
						rxnet.write('TS%5s_%s    %8.3f ==>    MIN%5s_%s  %4s      %7.3f    <-->   %-18s %4s +%4s      %8.3f  %4s\n' % (col[1],i+1,float(col[4])+ediff,col[7],i+1,klabel[int(col[7])],emin1+ediff,pr,klabel[int(col[10]),0],klabel[int(col[10]),1],epr,deg))
						kc+=1
						ts[kc]=col[1]
						r1[kc]=klabel[int(col[7])]
						p1[kc]=klabel[int(col[10]),0]
						p2[kc]=klabel[int(col[10]),1]
	ikc=0
	for line2 in range(len(ratefile)):
		col=ratefile[line2].split()
		ikc+=1
		if klabel.has_key(int(col[2])):
			rates.write('%30s %4s    0 %4s    0\n' % (col[0],r1[ikc],p1[ikc]))
			correl.write('%30s %4s    0 %4s    0 TS%5s_%s\n' % (col[0],r1[ikc],p1[ikc],ts[ikc],i+1))
		else:
			rates.write('%30s %4s    0 %4s %4s\n' % (col[0],r1[ikc],p1[ikc],p2[ikc]))
			correl.write('%30s %4s    0 %4s %4s TS%5s_%s\n' % (col[0],r1[ikc],p1[ikc],p2[ikc],ts[ikc],i+1))
rxnet.close()	
rates.close()
correl.close()
enadf.close()
