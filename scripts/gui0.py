from Tkinter import *
import sys
import os 

#funciones de procesamiento
def myfunction(event):
    canvas.configure(scrollregion=canvas.bbox("all"),width=70,height=90)

def myfunction2(event):
    canvas2.configure(scrollregion=canvas2.bbox("all"),width=70,height=90)

def fileread(name):
        with open(name) as f:
                line=f.read().splitlines()
        return line

###########3
# load the listbox with data
def load_item():
# Get the number of atoms of the molecule
	listbox1.delete(0,END)
	if job.get()==1:
		name=molecule.get()
		xyz=fileread(name+'.xyz')
		natom.set(int(xyz[0]))
	if sampling.get()==1:
		for item in range(1,natom.get()+1):
    			listbox1.insert(END, item)

def load_item2():
# Get the number of atoms of the molecule
	listbox2.delete(0,END)
	if job.get()==1:
		name=molecule.get()
		xyz=fileread(name+'.xyz')
		nmode.set(3*int((xyz[0]))-6)
	if sampling.get()==2:
		for item in range(1,nmode.get()+1):
    			listbox2.insert(END, item)

def delete_item():
	try:
		index=listbox1.curselection()[0]
 		listbox1.delete(index)
    	except IndexError:
        	pass
 
def delete_item2():
	try:
		index=listbox2.curselection()[0]
 		listbox2.delete(index)
    	except IndexError:
        	pass
#######33k
def proc1():
	if sampling.get()==1:
		name=molecule.get()
		samp="canonical"
	elif sampling.get()==2:
		name=molecule.get()
		samp="microcanonical"
	else:
		name=fr1.get()+'_'+fr2.get()
		samp="association"
	fname=name+'_'+samp+'.dat'			
	inpf=open(fname,'w')
	inpf.write('--General section--\n')
	if job.get()==1:
		inpf.write('molecule '+str(name)+'\n')
		inpf.write('charge '+str(charge.get())+'\n')
		inpf.write('mult '+str(mult.get())+'\n')
	else:
		inpf.write('molecule '+str(name)+'\n')
		inpf.write('charge '+str(chcat.get())+'\n')
		inpf.write('mult '+str(mucat.get())+'\n')
	inpf.write('LLcalc '+str(llcalc.get())+'\n')
        sp2corr=spcorr.get()
	if (sp2corr=="CCSD(T)"):
		sp2corr="CCSDT"
	if (sp2corr=="NO"):
		inpf.write('HLcalc 1'+'\t'+str(hlcalc1.get())+'\n')
		inpf.write('level1'+'\t'+str(level1.get())+'/'+str(basis1.get())+'\n')
	else:
		inpf.write('HLcalc 2'+'\t'+str(hlcalc1.get())+'\t'+str(sp2corr)+'\n')
		inpf.write('level1'+'\t'+str(level1.get())+'/'+str(basis1.get())+'\n')
		inpf.write('level2'+'\t'+str(level2.get())+'/'+str(basis2.get())+'\n')
        if iop.get()!= "":
		inpf.write('iop'+'\t'+str(iop.get())+'\n')

	if job.get()==2:
		inpf.write('\n--Catalysis section--\n')
		inpf.write('species '+str(species1.get())+'\t'+str(species2.get())+'\t'+str(species3.get())+'\t'+str(species4.get())+'\n')
		inpf.write('eta '+str(eta.get())+'\n')

	inpf.write('\n--CDS section--\n')
	if sampling.get()<=2:
		if job.get()==1:	
			xyz=fileread(name+'.xyz')
			natom.set(int(xyz[0]))
        	nmode.set(3*natom.get()-6)
	if sampling.get() ==1:
		inpf.write('sampling canonical\n')
		inpf.write('temp'+'\t'+str(tdyn.get())+'\n')
###
		temp_list = list(listbox1.get(0, END))
		temp_list = [str(item) for item in temp_list]
 		inpf.write('atoms ')
		if natom.get()==len(temp_list) or len(temp_list)==0:
 			inpf.write('all\n')
		else:
 			inpf.write(str(len(temp_list))+' ')
			for i in range(len(temp_list)-1):
 				inpf.write(temp_list[i]+',')
 			inpf.write(temp_list[len(temp_list)-1]+'\n')
###
	if sampling.get() ==2:
		inpf.write('sampling microcanonical\n')
		inpf.write('etraj'+'\t'+str(edyn.get())+'\n')
###
		temp_list = list(listbox2.get(0, END))
		temp_list = [str(item) for item in temp_list]
 		inpf.write('modes ')
		if nmode.get()==len(temp_list) or len(temp_list)==0:
 			inpf.write('all\n')
		else:
 			inpf.write(str(len(temp_list))+' ')
			for i in range(len(temp_list)-1):
 				inpf.write(temp_list[i]+',')
 			inpf.write(temp_list[len(temp_list)-1]+'\n')
	if sampling.get() ==3:
		inpf.write('sampling association\n')
		inpf.write('A= '+str(fr1.get())+'\n')
		inpf.write('B= '+str(fr2.get())+'\n')
		inpf.write('rotate '+str(pivot1.get())+'\t'+str(pivot2.get())+'\t'+str(dbp.get())+'\t'+str(mind.get())+'\n')

	if sampling.get()<=2:
		inpf.write('ntraj '+str(ntraj.get())+'\n')
		inpf.write('fs '+str(time.get())+'\n')
		if seed.get() != "random":
			inpf.write('seed '+str(seed.get())+'\n')
		if sampling.get()==1 and thmass.get() != 0:
			inpf.write('thmass '+str(thmass.get())+'\n')
		inpf.write('\n--BBFS section--\n')
		inpf.write('emaxts '+str(emaxts.get())+'\n')
		inpf.write('emints '+str(emints.get())+'\n')
		inpf.write('freqmin '+str(freqmin.get())+'\n')
		if fastmode.get()==1:
			inpf.write('fastmode\n')
		if refdist.get()==1:
			inpf.write('NOcreatethdist\n')
			cmd='createthdist0.sh '+str(name)+'\t'+str(job.get())
			os.system(cmd)
	inpf.write('\n--Screening of the structures section--\n')
	inpf.write('avgerr 1e-4'+'\n')
	inpf.write('bigerr 1'+'\n')
	inpf.write('thdiss 0.1'+'\n')
	if nvdw.get()>0:
		inpf.write('vdw '+str(nvdw.get())+'\n')
		inpf.write(str(atp1.get())+'\t'+str(atpd1.get())+'\n')
		if nvdw.get()==2:
			inpf.write(str(atp2.get())+'\t'+str(atpd2.get())+'\n')
		if nvdw.get()==3:
			inpf.write(str(atp2.get())+'\t'+str(atpd2.get())+'\n')
			inpf.write(str(atp3.get())+'\t'+str(atpd3.get())+'\n')
		
	
	if sampling.get()<=2:
		inpf.write('\n--Kinetics section--\n')
		if samplingkmc.get()==1:
			inpf.write('Rate canonical\n')
			inpf.write('TKMC '+str(tkmc.get())+'\n')
		else:
			inpf.write('Rate microcanonical\n')
			inpf.write('EKMC '+str(ekmc.get())+'\n')
		inpf.write('MaxEn '+str(emaxkmc.get())+'\n')
		inpf.write('nmol '+str(nmol.get())+'\n')
		if min0.get()=="":
			inpf.write('imin min0'+'\n')
		else:	
			inpf.write('imin '+str(min0.get())+'\n')
		if job.get()==2:
			spc=[]	
			spc0={}
			nespconst=0	
			if c01.get()<0:
				nespconst+=1
				spc0[nespconst]=sp1.get()
				spc=spc+spc0.values()
			if c02.get()<0:
				nespconst+=1
				spc0[nespconst]=sp2.get()
				spc=spc+spc0.values()
			if c03.get()<0:
				nespconst+=1
				spc0[nespconst]=sp3.get()
				spc=spc+spc0.values()
			if c04.get()<0:
				nespconst+=1
				spc0[nespconst]=sp4.get()
				spc=spc+spc0.values()
			if nspecies.get()==1:
				if sol1.get()>0:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+' '+str(sol1.get())+'\n')
				else:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+'\n')
			if nspecies.get()==2:
				if sol1.get()>0:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+' '+str(sol1.get())+'\n')
				else:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+'\n')
				if sol2.get()>0:
					inpf.write('c0 '+str(sp2.get())+' '+str(c02.get())+' ' +str(c02m.get())+' '+str(sol2.get())+'\n')
				else:
					inpf.write('c0 '+str(sp2.get())+' '+str(c02.get())+' ' +str(c02m.get())+'\n')
			if nspecies.get()==3:
				if sol1.get()>0:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+' '+str(sol1.get())+'\n')
				else:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+'\n')
				if sol2.get()>0:
					inpf.write('c0 '+str(sp2.get())+' '+str(c02.get())+' ' +str(c02m.get())+' '+str(sol2.get())+'\n')
				else:
					inpf.write('c0 '+str(sp2.get())+' '+str(c02.get())+' ' +str(c02m.get())+'\n')
				if sol3.get()>0:
					inpf.write('c0 '+str(sp3.get())+' '+str(c03.get())+' ' +str(c03m.get())+' '+str(sol3.get())+'\n')
				else:
					inpf.write('c0 '+str(sp3.get())+' '+str(c03.get())+' ' +str(c03m.get())+'\n')
			if nspecies.get()==4:
				if sol1.get()>0:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+' '+str(sol1.get())+'\n')
				else:
					inpf.write('c0 '+str(sp1.get())+' '+str(c01.get())+' ' +str(c01m.get())+'\n')
				if sol2.get()>0:
					inpf.write('c0 '+str(sp2.get())+' '+str(c02.get())+' ' +str(c02m.get())+' '+str(sol2.get())+'\n')
				else:
					inpf.write('c0 '+str(sp2.get())+' '+str(c02.get())+' ' +str(c02m.get())+'\n')
				if sol3.get()>0:
					inpf.write('c0 '+str(sp3.get())+' '+str(c03.get())+' ' +str(c03m.get())+' '+str(sol3.get())+'\n')
				else:
					inpf.write('c0 '+str(sp3.get())+' '+str(c03.get())+' ' +str(c03m.get())+'\n')
				if sol4.get()>0:
					inpf.write('c0 '+str(sp4.get())+' '+str(c04.get())+' ' +str(c04m.get())+' '+str(sol4.get())+'\n')
				else:
					inpf.write('c0 '+str(sp4.get())+' '+str(c04.get())+' ' +str(c04m.get())+'\n')
			if nespconst>0:
				del spc[0]
				inpf.write('nespconst '+str(nespconst)+' '+' '.join(spc)+'\n')
		inpf.write('Maxtime '+str(timekmc.get())+'\n')
		inpf.write('Stepsize '+str(dt.get())+'\n')
		if job.get()==1:
			inpf.write('ImpPaths '+str(percentrel.get())+'\n')
		inpf.write('PathInfo All\n')
		if product.get()!="":
			inpf.write('Product '+str(product.get())+'\n')
			

	inpf.close()

def proc22kmc():
	if samplingkmc.get()==1 and sampling.get()<=2:
		tkmc.set("298")
		ekmc.set("")
		nmol.set("1000")
		timekmc.set("1e12")
		dt.set("1e11")
		emaxkmc.set("40")
	if samplingkmc.get()==2 and sampling.get()<=2:
		tkmc.set("")
		ekmc.set("150")
		nmol.set("1000")
		timekmc.set("10")
		dt.set("0.1")
		emaxkmc.set("200")
	if sampling.get()==3:
		tkmc.set("")
		ekmc.set("")
		nmol.set("")
		timekmc.set("")
		dt.set("")
		emaxkmc.set("")
		

def proc22():
	if sampling.get()==1:
		tdyn.set("10000")
		edyn.set("")
		fr1.set("")
		fr2.set("")
		pivot1.set("")
		pivot2.set("")
		dbp.set("")
		mind.set("")
		emaxts.set("200")
		emints.set("-200")
		freqmin.set("200")
		nmol.set("1000")
		ntraj.set("1")
		time.set("500")
		seed.set("random")
		thmass.set("0.0")
	elif sampling.get()==2:
		tdyn.set("")
		edyn.set("400")
		emaxts.set("400")
		fr1.set("")
		fr2.set("")
		pivot1.set("")
		pivot2.set("")
		dbp.set("")
		mind.set("")
		freqmin.set("200")
		emints.set("-200")
		nmol.set("1000")
		ntraj.set("1")
		time.set("500")
		seed.set("random")
		thmass.set("0.0")
	else:
		tdyn.set("")
		edyn.set("")
		fr1.set("name")
		fr2.set("name")
		pivot1.set("com")
		pivot2.set("com")
		dbp.set("2.0")
		mind.set("1.0")
 		emaxts.set("")
 		emints.set("")
		freqmin.set("")
		tkmc.set("")
		ekmc.set("")
		nmol.set("")
		timekmc.set("")
		dt.set("")
		emaxkmc.set("")
		ntraj.set("")
		time.set("")
		seed.set("")
		thmass.set("")
				

def proc2():
	if job.get()==1:
		min0.set("ask later")
		species1.set("")
		species2.set("")
		species3.set("")
		species4.set("")
		ss1.set("")
		ss2.set("")
		ss3.set("")
		ss4.set("")
		ss5.set("")
		ss6.set("")
		ss7.set("")
		ss8.set("")
		eta.set("")
		charge.set("0")
		charge1.set("")
		charge2.set("")
		charge3.set("")
		charge4.set("")
		mult.set("1")
		mult1.set("")
		mult2.set("")
		mult3.set("")
		mult4.set("")
	else:
		min0.set("")
                charge.set("")
                charge1.set("0")
                charge2.set("0")
                charge3.set("0")
                charge4.set("0")
                mult.set("")
		mult1.set("1")
		mult2.set("1")
		mult3.set("1")
		mult4.set("1")
		species1.set("cat")
		eta.set("0.01")
def proc3():
        sp1.set(species1.get())
        sp2.set(species2.get())
        sp3.set(species3.get())
        sp4.set(species4.get())
        if job.get()==1:
		return
	ss1.set("")
	ss2.set("")
	ss3.set("")
	ss4.set("")
	ss5.set("")
	ss6.set("")
	ss7.set("")
	ss8.set("")
	min_name={}
	nmolecule=0
	if not species1.get()=="":
		nmolecule+=1			     	 
		min_name[1]=species1.get()
	if not species2.get()=="":
		nmolecule+=1			     	 
		min_name[2]=species2.get()
	if not species3.get()=="":
		nmolecule+=1			     	 
		min_name[3]=species3.get()
	if not species4.get()=="":
		nmolecule+=1			     	 
		min_name[4]=species4.get()
        nspecies.set(nmolecule)
	if nmolecule==1:
		return
        subsystem=['cat']
        pairs={}
        i=1
        p=0
        for i in range(2,nmolecule+1):
                p+=1
                pairs[p]='cat_'+min_name[i]
        subsystem=subsystem+pairs.values()
	ss1.set(subsystem[0])
	ss2.set(subsystem[1])
        t=0
	if nmolecule>2:
        	trios={}
	        for i in range(1,p+1):
                	for j in range(i+2,nmolecule+1):
                        	t+=1
                        	trios[t]=pairs[i]+'_'+min_name[j]
        	subsystem=subsystem+trios.values()
		ss3.set(subsystem[2])
		ss4.set(subsystem[3])
		if nmolecule==4:
        		fours={}
        		fours[1]=trios[1]+'_'+min_name[nmolecule]
        		subsystem=subsystem+fours.values()
			ss5.set(subsystem[4])
			ss6.set(subsystem[5])
			ss7.set(subsystem[6])
			ss8.set(subsystem[7])
	i=0
	for txt in subsystem:
		i+=1
		Radiobutton(f02,text=txt,font=("Helvetica",11),indicatoron=0,variable=cs,value=i,width=25,command=proc4).grid(row=i,column=8,columnspan=2)

def proc4():
	if species2.get()=="":
		nos=1
	elif species3.get()=="":
		nos=2
	elif species4.get()=="":
		nos=3
	else:
		nos=4
	for i in range(1,9):
		if cs.get()==i:
			if i==1:
				ss.set(ss1.get())
				chcat.set(charge1.get())
				mucat.set(mult1.get())
				xyz1=fileread(species1.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=0
				natom.set(natom1)
				if sampling.get()==3:	
					fr1.set(species1.get())
					fr2.set("")
			if i==2:
				ss.set(ss2.get())
				chcat.set(charge1.get()+charge2.get())
				mucat.set(mult1.get()+mult2.get()-1)
				xyz1=fileread(species1.get()+'.xyz')
				xyz2=fileread(species2.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=(int(xyz2[0]))
				natom.set(natom1+natom2)
				if sampling.get()==3:	
					fr1.set(species1.get())
					fr2.set(species2.get())
			if i==3:
				ss.set(ss3.get())
				chcat.set(charge1.get()+charge3.get())
				mucat.set(mult1.get()+mult3.get()-1)
				xyz1=fileread(species1.get()+'.xyz')
				xyz2=fileread(species3.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=(int(xyz2[0]))
				natom.set(natom1+natom2)
				if sampling.get()==3:	
					fr1.set(species1.get())
					fr2.set(species3.get())
			if i==4:
				ss.set(ss4.get())
				if nos==3:
					chcat.set(charge1.get()+charge2.get()+charge3.get())
					mucat.set(mult1.get()+mult2.get()+mult3.get()-2)
					xyz1=fileread(species1.get()+'.xyz')
					xyz2=fileread(species2.get()+'.xyz')
					xyz3=fileread(species3.get()+'.xyz')
					natom1=(int(xyz1[0]))
					natom2=(int(xyz2[0]))
					natom3=(int(xyz3[0]))
					natom.set(natom1+natom2+natom3)
					if sampling.get()==3:	
						fr1.set(species1.get()+'_'+species2.get())
						fr2.set(species3.get())
				elif nos==4:
					chcat.set(charge1.get()+charge4.get())
					mucat.set(mult1.get()+mult4.get()-1)
					xyz1=fileread(species1.get()+'.xyz')
					xyz2=fileread(species4.get()+'.xyz')
					natom1=(int(xyz1[0]))
					natom2=(int(xyz2[0]))
					natom.set(natom1+natom2)
					if sampling.get()==3:	
						fr1.set(species1.get())
						fr2.set(species4.get())
			if i==5:
				ss.set(ss5.get())
				chcat.set(charge1.get()+charge2.get()+charge3.get())
				mucat.set(mult1.get()+mult2.get()+mult3.get()-2)
				xyz1=fileread(species1.get()+'.xyz')
				xyz2=fileread(species2.get()+'.xyz')
				xyz3=fileread(species3.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=(int(xyz2[0]))
				natom3=(int(xyz3[0]))
				natom.set(natom1+natom2+natom3)
				if sampling.get()==3:	
					fr1.set(species1.get()+'_'+species2.get())
					fr2.set(species3.get())
			if i==6:
				ss.set(ss6.get())
				chcat.set(charge1.get()+charge2.get()+charge4.get())
				mucat.set(mult1.get()+mult2.get()+mult4.get()-2)
				xyz1=fileread(species1.get()+'.xyz')
				xyz2=fileread(species2.get()+'.xyz')
				xyz3=fileread(species4.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=(int(xyz2[0]))
				natom3=(int(xyz3[0]))
				natom.set(natom1+natom2+natom3)
				if sampling.get()==3:	
					fr1.set(species1.get()+'_'+species2.get())
					fr2.set(species4.get())
			if i==7:
				ss.set(ss7.get())
				chcat.set(charge1.get()+charge3.get()+charge4.get())
				mucat.set(mult1.get()+mult3.get()+mult4.get()-2)
				xyz1=fileread(species1.get()+'.xyz')
				xyz2=fileread(species3.get()+'.xyz')
				xyz3=fileread(species4.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=(int(xyz2[0]))
				natom3=(int(xyz3[0]))
				natom.set(natom1+natom2+natom3)
				if sampling.get()==3:	
					fr1.set(species1.get()+'_'+species3.get())
					fr2.set(species4.get())
			if i==8:
				ss.set(ss8.get())
				chcat.set(charge1.get()+charge2.get()+charge3.get()+charge4.get())
				mucat.set(mult1.get()+mult2.get()+mult3.get()+mult4.get()-3)
				xyz1=fileread(species1.get()+'.xyz')
				xyz2=fileread(species2.get()+'.xyz')
				xyz3=fileread(species3.get()+'.xyz')
				xyz4=fileread(species4.get()+'.xyz')
				natom1=(int(xyz1[0]))
				natom2=(int(xyz2[0]))
				natom3=(int(xyz3[0]))
				natom4=(int(xyz4[0]))
				natom.set(natom1+natom2+natom3+natom4)
				if sampling.get()==3:	
					fr1.set(species1.get()+'_'+species2.get()+'_'+species3.get())
					fr2.set(species4.get())
	molecule.set(ss.get())
	charge.set(chcat.get())
	mult.set(mucat.get())
        nmode.set(3*natom.get()-6)
	

def proc5():
	ventana.destroy()


#Instancia de la clase Tk
ventana = Tk()
ventana.title('AMK Graphical User Interface')

f0=Frame(ventana,width=1100,height=235,relief=RIDGE,bd=3)
f0.grid(row=0,column=0,columnspan=10,rowspan=7)
f0.grid_propagate(False)
f01=Frame(f0,width=270,height=230,relief=RIDGE,bd=3)
f01.grid(row=0,column=0,columnspan=2,rowspan=7)
f01.grid_propagate(False)
f02=Frame(f0,width=830,height=230,relief=RIDGE,bd=3)
f02.grid(row=0,column=2,columnspan=8,rowspan=7)
f02.grid_propagate(False)

f1=Frame(ventana,width=465,height=105,relief=RIDGE,bd=3)
f1.grid(row=10,column=0,columnspan=4,rowspan=3)
f1.grid_propagate(False)
f2=Frame(ventana,width=635,height=105,relief=RIDGE,bd=3)
f2.grid(row=10,column=5,columnspan=6,rowspan=5)
f2.grid_propagate(False)

f3=Frame(ventana,width=1100,height=180,relief=RIDGE,bd=3)
f3.grid(row=15,column=0,columnspan=10,rowspan=8)
f3.grid_propagate(False)
f31=Frame(f3,width=635,height=155,relief=RIDGE,bd=3)
f31.grid(row=16,column=0,columnspan=5,rowspan=4)
f31.grid_propagate(False)
f311=Frame(f31,width=210,height=150,relief=RIDGE,bd=3)
f311.grid(row=16,column=0,columnspan=2,rowspan=5)
f311.grid_propagate(False)
fcanvas=Frame(f311,width=50,height=50,relief=RIDGE,bd=3)
fcanvas.grid(row=18,column=0,columnspan=1,rowspan=3)
fcanvas.grid_propagate(False)
f312=Frame(f31,width=210,height=150,relief=RIDGE,bd=3)
f312.grid(row=16,column=2,columnspan=2,rowspan=5)
f312.grid_propagate(False)
fcanvas2=Frame(f312,width=50,height=50,relief=RIDGE,bd=3)
fcanvas2.grid(row=18,column=2,columnspan=1,rowspan=3)
fcanvas2.grid_propagate(False)
f313=Frame(f31,width=210,height=150,relief=RIDGE,bd=3)
f313.grid(row=16,column=4,columnspan=2,rowspan=4)
f313.grid_propagate(False)

f32=Frame(f3,width=465,height=155,relief=RIDGE,bd=3)
f32.grid(row=16,column=6,columnspan=4,rowspan=4,sticky=N)
f32.grid_propagate(False)

f4=Frame(ventana,width=1100,height=160,relief=RIDGE,bd=3)
f4.grid(row=33,column=0,columnspan=10,rowspan=7)
f4.grid_propagate(False)
f41=Frame(f4,width=1100,height=55,relief=RIDGE,bd=3)
f41.grid(row=34,column=0,columnspan=10,rowspan=2)
f41.grid_propagate(False)
f421=Frame(f4,width=270,height=80,relief=RIDGE,bd=3)
f421.grid(row=36,column=0,columnspan=2,rowspan=3)
f421.grid_propagate(False)
f422=Frame(f4,width=830,height=80,relief=RIDGE,bd=3)
f422.grid(row=36,column=2,columnspan=8,rowspan=3)
f422.grid_propagate(False)

f5=Frame(ventana,width=1100,height=40,relief=RIDGE,bd=3)
f5.grid(row=40,column=0,columnspan=10,rowspan=1)
f5.grid_propagate(False)

#f5=Frame(f0,width=500,height=500,relief=SUNKEN,bd=2)
#f5.grid(row=18,column=0,columnspan=4,rowspan=2)

#Variables que almacenarn los datos
product=StringVar()
species1 = StringVar()
species2 = StringVar()
species3 = StringVar()
species4 = StringVar()
nvdw=IntVar()
min0=StringVar()
sp1=StringVar()
sp2=StringVar()
sp3=StringVar()
sp4=StringVar()
c01=DoubleVar()
c02=DoubleVar()
c03=DoubleVar()
c04=DoubleVar()
c01m=DoubleVar()
c02m=DoubleVar()
c03m=DoubleVar()
c04m=DoubleVar()
sol1=DoubleVar()
sol2=DoubleVar()
sol3=DoubleVar()
sol4=DoubleVar()
fastmode=IntVar()
refdist=IntVar()
nspecies=IntVar()
atp1=StringVar()
atpd1=DoubleVar()
atp2=StringVar()
atpd2=DoubleVar()
atp3=StringVar()
atpd3=DoubleVar()
charge = IntVar()
charge1 = IntVar()
charge2 = IntVar()
charge3 = IntVar()
charge4 = IntVar()
natom=IntVar()
nmode=IntVar()
mult = IntVar()
mult1 = IntVar()
mult2 = IntVar()
mult3 = IntVar()
mult4 = IntVar()
llcalc = StringVar()
nhlcalc = IntVar()
spcorr=StringVar()
hlcalc1 = StringVar()
title1 = StringVar()
title2 = StringVar()
title3 = StringVar()
title4 = StringVar()
level1 = StringVar()
basis1 = StringVar()
level2 = StringVar()
iop = StringVar()
basis2 = StringVar()
eta = DoubleVar()
tdyn = IntVar()
tkmc = IntVar()
ekmc = IntVar()
edyn = StringVar()
fr1 = StringVar()
fr2 = StringVar()
pivot1 = StringVar()
pivot2 = StringVar()
dbp = DoubleVar()
mind = DoubleVar()
sampling = IntVar()
samplingkmc = IntVar()
ntraj = IntVar()
seed = StringVar()
thmass=DoubleVar()
emaxts = IntVar()
emaxkmc = IntVar()
percentrel=DoubleVar()
emints = IntVar()
freqmin = IntVar()
time = IntVar()
molecule=StringVar()
ss1=StringVar()
ss2=StringVar()
ss3=StringVar()
ss4=StringVar()
ss5=StringVar()
ss6=StringVar()
ss7=StringVar()
ss8=StringVar()
ss=StringVar()
chcat=IntVar()
mucat=IntVar()
job=IntVar()
cs=IntVar()
nmol=IntVar()
timekmc=DoubleVar()
dt=DoubleVar()
#genero.set(1)
sp1.set("")
sp2.set("")
sp3.set("")
sp4.set("")
percentrel.set("0.1")
timekmc.set("1e12")
dt.set("1e11")
nmol.set("1000")
eta.set("")
tkmc.set("298")
ekmc.set("")
tdyn.set("10000")
edyn.set("")
charge.set("0")
atp1.set("")
atpd1.set("")
atpd2.set("")
atpd3.set("")
charge1.set("")
charge2.set("")
charge3.set("")
charge4.set("")
min0.set("ask later")
fr1.set("")
fr2.set("")
dbp.set("")
mind.set("")
pivot1.set("")
pivot2.set("")
llcalc.set("PM7")
hlcalc1.set("DFT")
level1.set("B3LYP")
basis1.set("6-31G*")
mult.set("1")
mult1.set("")
mult2.set("")
mult3.set("")
mult4.set("")
spcorr.set("NO")
seed.set("random")
thmass.set("0.0")
emaxts.set("200")
emaxkmc.set("40")
emints.set("-200")
freqmin.set("200")
time.set("500")
job.set("1")
molecule.set("")
sampling.set("1")
samplingkmc.set("1")
ntraj.set("1")


#title Label
#Label(f0, text='Job Type',font=("Helvetica",14),height=1,anchor=S).grid(row=0, column=0,columnspan=10)


#Radiobutton
Radiobutton(f01,variable=job,text="Single system",font=("Helvetica",11),value=1,command=proc2).grid(row=0,column=0,columnspan=2)
Radiobutton(f02,variable=job,text="Catalysis",font=("Helvetica",11),value=2,command=proc2).grid(row=0,column=2,columnspan=7)
#Button
Button(f02, text='Get subsystems',font=("Helvetica",11),command=proc3).grid(row=0, column=8,columnspan=2)

#title
Label(f01, text='Molecule: ',font=("Helvetica",11),width=12).grid(row=1, column=0,sticky=E)
Entry(f01, textvariable=molecule,font=("Helvetica",11),width=15,bg="salmon").grid(row=1, column=1,sticky=W)

Label(f01, text='Charge: ',font=("Helvetica",11),width=12).grid(row=2, column=0,sticky=E)
Entry(f01, textvariable=charge,font=("Helvetica",11),width=2,bg="salmon").grid(row=2, column=1,sticky=W)
Label(f01, text='Multiplicity: ',font=("Helvetica",11),width=12).grid(row=3, column=0,sticky=E)
Entry(f01, textvariable=mult,font=("Helvetica",11),width=2,bg="salmon").grid(row=3, column=1,sticky=W)

#title

#species
Label(f02, text='Molecule 1: ',font=("Helvetica",11),width=12).grid(row=1, column=2,sticky=E)
Label(f02, text='Molecule 2: ',font=("Helvetica",11),width=12).grid(row=2, column=2,sticky=E)
Label(f02, text='Molecule 3: ',font=("Helvetica",11),width=12).grid(row=3, column=2,sticky=E)
Label(f02, text='Molecule 4: ',font=("Helvetica",11),width=12).grid(row=4, column=2,sticky=E)
Entry(f02, textvariable=species1,font=("Helvetica",11),width=12,bg="pale green").grid(row=1, column=3,sticky=W)
Entry(f02, textvariable=species2,font=("Helvetica",11),width=12,bg="pale green").grid(row=2, column=3,sticky=W)
Entry(f02, textvariable=species3,font=("Helvetica",11),width=12,bg="pale green").grid(row=3, column=3,sticky=W)
Entry(f02, textvariable=species4,font=("Helvetica",11),width=12,bg="pale green").grid(row=4, column=3,sticky=W)
Label(f02, text='Charge: ',font=("Helvetica",11),      width=8).grid(row=1, column=4,sticky=E)
Entry(f02, textvariable=charge1,font=("Helvetica",11), width= 2,bg="pale green").grid(row=1, column=5,sticky=W)
Label(f02, text='Multiplicity: ',font=("Helvetica",11),width=10).grid(row=1, column=6,sticky=E)
Entry(f02, textvariable=mult1,font=("Helvetica",11),   width= 2,bg="pale green").grid(row=1, column=7,sticky=W)
Label(f02, text='Charge: ',font=("Helvetica",11),      width=8).grid(row=2, column=4,sticky=E)
Entry(f02, textvariable=charge2,font=("Helvetica",11), width= 2,bg="pale green").grid(row=2, column=5,sticky=W)
Label(f02, text='Multiplicity: ',font=("Helvetica",11),width=10).grid(row=2, column=6,sticky=E)
Entry(f02, textvariable=mult2,font=("Helvetica",11),   width= 2,bg="pale green").grid(row=2, column=7,sticky=W)
Label(f02, text='Charge: ',font=("Helvetica",11),      width=8).grid(row=3, column=4,sticky=E)
Entry(f02, textvariable=charge3,font=("Helvetica",11), width= 2,bg="pale green").grid(row=3, column=5,sticky=W)
Label(f02, text='Multiplicity: ',font=("Helvetica",11),width=10).grid(row=3, column=6,sticky=E)
Entry(f02, textvariable=mult3,font=("Helvetica",11),   width= 2,bg="pale green").grid(row=3, column=7,sticky=W)
Label(f02, text='Charge: ',font=("Helvetica",11),      width=8).grid(row=4, column=4,sticky=E)
Entry(f02, textvariable=charge4,font=("Helvetica",11), width= 2,bg="pale green").grid(row=4, column=5,sticky=W)
Label(f02, text='Multiplicity: ',font=("Helvetica",11),width=10).grid(row=4, column=6,sticky=E)
Entry(f02, textvariable=mult4,font=("Helvetica",11),   width= 2,bg="pale green").grid(row=4, column=7,sticky=W)
Label(f02, text=u'\u03B7 (Pas): ',font=("Helvetica",11),width=12).grid(row=5, column=2,sticky=E)
Entry(f02, textvariable=eta,font=("Helvetica",11),width=12,bg="pale green").grid(row=5, column=3,sticky=W)

#title Label
Label(f1, text='MOPAC calc: ',font=("Helvetica",11),width=12).grid(row=10, column=1,sticky=E)
OptionMenu(f1, llcalc, "AM1","PM3","PM6","PM7").grid(row=10,column=2,sticky=W)
Label(f1, text='G09 calc: ',font=("Helvetica",11),width=12).grid(row=11, column=0,sticky=E)
OptionMenu(f1, hlcalc1, "HF","MP2","DFT").grid(row=11, column=1,sticky=W)
Entry(f1, textvariable=level1,font=("Helvetica",11),width=12).grid(row=12, column=0,sticky=W)
Entry(f1, textvariable=basis1,font=("Helvetica",11),width=12).grid(row=12, column=1,sticky=W)
Label(f1, text='SP correction?',font=("Helvetica",11),width=12).grid(row=11, column=2)
OptionMenu(f1, spcorr, "NO","HF","MP2","DFT","CCSD(T)").grid(row=11, column=3)
Entry(f1, textvariable=level2,font=("Helvetica",11),width=12).grid(row=12, column=2,sticky=W)
Entry(f1, textvariable=basis2, font=("Helvetica",11),width=12).grid(row=12, column=3,sticky=W)
Label(f1, text='IOP',font=("Helvetica",11),width=12).grid(row=13, column=0)
Entry(f1, textvariable=iop, font=("Helvetica",11),width=12).grid(row=13, column=1,sticky=W)

#title Label
Label(f3, text='Sampling method',font=("Helvetica",11),height=0,anchor=S).grid(row=15, column=0,columnspan=10)

#sampling OptionMenu
Radiobutton(f311,variable=sampling,text="Canonical",font=("Helvetica",11),value=1,command=proc22).grid(row=16,column=0,columnspan=2)
Radiobutton(f312,variable=sampling,text="Microcanonical",font=("Helvetica",11),value=2,command=proc22).grid(row=16,column=2,columnspan=2)
Radiobutton(f32,variable=sampling,text="Association",font=("Helvetica",11),value=3,command=proc22).grid(row=16,column=6,columnspan=4)

Label(f311, text='Temp (K):',font=("Helvetica",11),width=8).grid(row=17, column=0,sticky=E)
Entry(f311, textvariable=tdyn,font=("Helvetica",11),width=8).grid(row=17, column=1,sticky=W)
###############33333333
# create the listbox (note that size is in characters)
canvas=Canvas(fcanvas)
canvas2=Canvas(fcanvas2)
frame=Frame(canvas)
frame2=Frame(canvas2)
listbox1 = Listbox(frame, width=50, height=natom.get(),selectmode='multiple')
listbox1.grid(row=18, column=0,columnspan=1,rowspan=3)

myscrollbar = Scrollbar(fcanvas,command=canvas.yview, orient=VERTICAL)
canvas.configure(yscrollcommand=myscrollbar.set)
myscrollbar.pack(side="right",fill="y")
canvas.pack(side="right")
canvas.create_window((0,0),window=frame,anchor='nw')
frame.bind("<Configure>",myfunction)
Label(f311, text='Excite subset:',font=("Helvetica",11),width=12).grid(row=18, column=1)
button1 = Button(f311, text='Load atoms', command=load_item)
button1.grid(row=19, column=1,sticky=W)
button2 = Button(f311, text='Delete frozen', command=delete_item)
button2.grid(row=20, column=1,sticky=W)
################33
listbox2 = Listbox(frame2, width=50, height=3*natom.get()-6,selectmode='multiple')
listbox2.grid(row=18, column=2,sticky=E)

myscrollbar2 = Scrollbar(fcanvas2,command=canvas2.yview, orient=VERTICAL)
canvas2.configure(yscrollcommand=myscrollbar2.set)
myscrollbar2.pack(side="right",fill="y")
canvas2.pack(side="right")
canvas2.create_window((0,0),window=frame2,anchor='nw')
frame2.bind("<Configure>",myfunction2)
 
Label(f312, text='Excite subset:',font=("Helvetica",11),width=12).grid(row=18, column=3)
button1 = Button(f312, text='Load modes', command=load_item2)
button1.grid(row=19, column=3,sticky=W)
button2 = Button(f312, text='Delete frozen', command=delete_item2)
button2.grid(row=20, column=3,sticky=W)
################33


###############
Label(f312, text='E (kcal): ',font=("Helvetica",11),width=8).grid(row=17, column=2,sticky=E)
Entry(f312, textvariable=edyn,font=("Helvetica",11),width=8).grid(row=17, column=3,sticky=W)
Label(f32, text='Fragment A: ',font=("Helvetica",11),width=16).grid(row=17, column=6,sticky=E)
Entry(f32, textvariable=fr1,font=("Helvetica",11),width=6).grid(row=17, column=7,sticky=W)
Label(f32, text='Fragment B: ',font=("Helvetica",11),width=16).grid(row=17, column=8,sticky=E)
Entry(f32, textvariable=fr2,font=("Helvetica",11),width=6).grid(row=17, column=9,sticky=W)
Label(f32, text='Pivot of A: ',font=("Helvetica",11),width=16).grid(row=18, column=6,sticky=E)
Entry(f32, textvariable=pivot1,font=("Helvetica",11),width=6).grid(row=18, column=7,sticky=W)
Label(f32, text='Pivot of B: ',font=("Helvetica",11),width=16).grid(row=18, column=8,sticky=E)
Entry(f32, textvariable=pivot2,font=("Helvetica",11),width=6).grid(row=18, column=9,sticky=W)
Label(f32, text=u"R betw pivots (\u212B):",font=("Helvetica",11),width=16).grid(row=19, column=6,sticky=E)
Entry(f32, textvariable=dbp,font=("Helvetica",11),width=6).grid(row=19, column=7,sticky=W)


Label(f32, text=u"RMIN betw frgs (\u212B):",font=("Helvetica",11),width=15).grid(row=19, column=8,sticky=E)
Entry(f32, textvariable=mind,font=("Helvetica",11),width=6).grid(row=19, column=9,sticky=W)

#spcorr Entry
Label(f313, text='Simulation details:',font=("Helvetica",11),width=20).grid(row=18, column=4,sticky=E,columnspan=2)
Label(f313, text='Trajectories:',font=("Helvetica",11),width=12).grid(row=19, column=4,sticky=E)
Entry(f313, textvariable=ntraj,font=("Helvetica",11),width=8).grid(row=19, column=5,sticky=W)
Label(f313, text='Time (fs):',font=("Helvetica",11),width=12).grid(row=20, column=4,sticky=E)
Entry(f313, textvariable=time,font=("Helvetica",11),width=8).grid(row=20, column=5,sticky=W)
Label(f313, text='Seed:',font=("Helvetica",11),width=12).grid(row=21, column=4,sticky=E)
Entry(f313, textvariable=seed,font=("Helvetica",11),width=8).grid(row=21, column=5,sticky=W)
Label(f313, text='Mass_exc.(au):',font=("Helvetica",11),width=12).grid(row=22, column=4,sticky=E)
Entry(f313, textvariable=thmass,font=("Helvetica",11),width=8).grid(row=22, column=5,sticky=W)

#title Label
#EMax and Emints and freqmin
Label(f2, text='EMAX_TS (kcal):',font=("Helvetica",11),width=15).grid(row=10, column=5,sticky=E)
Entry(f2, textvariable=emaxts,font=("Helvetica",11),width=7).grid(row=10, column=6,sticky=W)
Label(f2, text='EMIN_TS (kcal):',font=("Helvetica",11),width=15).grid(row=10, column=7,sticky=E)
Entry(f2, textvariable=emints,font=("Helvetica",11),width=7).grid(row=10, column=8,sticky=W)
Label(f2, text='FRMIN_TS (cm-1):',font=("Helvetica",11),width=15).grid(row=10, column=9,sticky=E)
Entry(f2, textvariable=freqmin,font=("Helvetica",11),width=7).grid(row=10, column=10,sticky=W)


Checkbutton(f2,variable=fastmode,text="Fast calculation",font=("Helvetica",11)).grid(row=12,column=5,columnspan=2)
Checkbutton(f2,variable=refdist,text="Modify reference distances",font=("Helvetica",11)).grid(row=13,column=5,columnspan=2)


Label(f2, text='Number add_dist:',font=("Helvetica",11),width=15).grid(row=12, column=7,sticky=E)
Label(f2, text=u'If the above number>=1',font=("Helvetica",11),width=22).grid(row=13, column=7,columnspan=2,sticky=E)
Label(f2, text=u'Atom pairs and dists. (\u212B)',font=("Helvetica",11),width=22).grid(row=14, column=7,columnspan=2,sticky=E)
Entry(f2, textvariable=nvdw,font=("Helvetica",11),width=7).grid(row=12, column=8,sticky=W)

Entry(f2, textvariable=atp1,font=("Helvetica",11),width=7).grid(row=12, column=9,sticky=E)
Entry(f2, textvariable=atpd1,font=("Helvetica",11),width=7).grid(row=12, column=10,sticky=W)
Entry(f2, textvariable=atp2,font=("Helvetica",11),width=7).grid(row=13, column=9,sticky=E)
Entry(f2, textvariable=atpd2,font=("Helvetica",11),width=7).grid(row=13, column=10,sticky=W)
Entry(f2, textvariable=atp3,font=("Helvetica",11),width=7).grid(row=14, column=9,sticky=E)
Entry(f2, textvariable=atpd3,font=("Helvetica",11),width=7).grid(row=14, column=10,sticky=W)


#title Label
Label(f4, text='Kinetic simulations',font=("Helvetica",11),height=0,anchor=S).grid(row=33, column=0,columnspan=10)
Radiobutton(f41,variable=samplingkmc,text="Canonical",font=("Helvetica",11),value=1,command=proc22kmc).grid(row=36,column=0,columnspan=2)
Radiobutton(f41,variable=samplingkmc,text="Microcanonical",font=("Helvetica",11),value=2,command=proc22kmc).grid(row=36,column=2,columnspan=2)
Label(f41, text='Temp (K):',font=("Helvetica",11),width=12).grid(row=37, column=0,sticky=E)
Entry(f41, textvariable=tkmc,font=("Helvetica",11),width=8).grid(row=37, column=1,sticky=W)
Label(f41, text='E (kcal):',font=("Helvetica",11),width=12).grid(row=37, column=2,sticky=E)
Entry(f41, textvariable=ekmc,font=("Helvetica",11),width=8).grid(row=37, column=3,sticky=W)

Label(f421, text='Starting Min:',font=("Helvetica",11),width=12).grid(row=38, column=0,sticky=E)
Entry(f421, textvariable=min0,font=("Helvetica",11),width=10,bg="salmon").grid(row=38, column=1,sticky=W)
Label(f421, text='Relevant Paths (%):',font=("Helvetica",11),width=18).grid(row=39, column=0,sticky=E)
Entry(f421, textvariable=percentrel,font=("Helvetica",11),width=10,bg="salmon").grid(row=39, column=1,sticky=W)

Label(f41, text='Number of molecules:',font=("Helvetica",11),width=20).grid(row=36, column=4,sticky=E)
Entry(f41, textvariable=nmol,font=("Helvetica",11),width=8).grid(row=36, column=5,sticky=W)
Label(f41, text='Time (s/ps):',font=("Helvetica",11),width=12).grid(row=36, column=6,sticky=E)
Entry(f41, textvariable=timekmc,font=("Helvetica",11),width=8).grid(row=36, column=7,sticky=W)
Label(f41, text=u'\u0394t (s/ps):',font=("Helvetica",11),width=12).grid(row=37, column=6,sticky=E)
Entry(f41, textvariable=dt,font=("Helvetica",11),width=8).grid(row=37, column=7,sticky=W)
Label(f41, text='MaxEn (kcal):',font=("Helvetica",11),width=18).grid(row=36, column=8,sticky=E)
Entry(f41, textvariable=emaxkmc,font=("Helvetica",11),width=8).grid(row=36, column=9,sticky=W)
Label(f41, text='Product:',font=("Helvetica",11),width=18).grid(row=37, column=8,sticky=E)
Entry(f41, textvariable=product,font=("Helvetica",11),width=8).grid(row=37, column=9,sticky=W)

Label(f422, text='Species',font=("Helvetica",11),width=11).grid(row=36, column=3,sticky=E)
Label(f422, text=u'C\u2080/P\u2080',font=("Helvetica",11),width=11).grid(row=36, column=4,sticky=E)
Label(f422, text=u'C\u2080max/P\u2080max',font=("Helvetica",11),width=11).grid(row=36, column=5,sticky=E)
Label(f422, text='Solubility',font=("Helvetica",11),width=11).grid(row=36, column=6,sticky=E)
Label(f422, text='Species',font=("Helvetica",11),width=11).grid(row=36, column=7,sticky=E)
Label(f422, text=u'C\u2080/P\u2080',font=("Helvetica",11),width=11).grid(row=36, column=8,sticky=E)
Label(f422, text=u'C\u2080max/P\u2080max',font=("Helvetica",11),width=11).grid(row=36, column=9,sticky=E)
Label(f422, text='Solubility',font=("Helvetica",11),width=11).grid(row=36, column=10,sticky=E)
Entry(f422, textvariable=sp1,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=3,sticky=W)
Entry(f422, textvariable=c01,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=4,sticky=W)
Entry(f422, textvariable=c01m,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=5,sticky=W)
Entry(f422, textvariable=sol1,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=6,sticky=W)
Entry(f422, textvariable=sp2,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=7,sticky=W)
Entry(f422, textvariable=c02,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=8,sticky=W)
Entry(f422, textvariable=c02m,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=9,sticky=W)
Entry(f422, textvariable=sol2,font=("Helvetica",11),width=11,bg="pale green").grid(row=37, column=10,sticky=W)
Entry(f422, textvariable=sp3,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=3,sticky=W)
Entry(f422, textvariable=c03,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=4,sticky=W)
Entry(f422, textvariable=c03m,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=5,sticky=W)
Entry(f422, textvariable=sol3,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=6,sticky=W)
Entry(f422, textvariable=sp4,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=7,sticky=W)
Entry(f422, textvariable=c04,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=8,sticky=W)
Entry(f422, textvariable=c04m,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=9,sticky=W)
Entry(f422, textvariable=sol4,font=("Helvetica",11),width=11,bg="pale green").grid(row=38, column=10,sticky=W)
#Label(ventana, text=u'\u0394t (s/ps):',font=("Helvetica",11),width=12).grid(row=35, column=8,sticky=E)
#Entry(ventana, textvariable=emaxkmc,font=("Helvetica",11),width=8).grid(row=35, column=9,sticky=W)

#title Label
Label(ventana, text='',font=("Helvetica",11,"bold"),height=1,anchor=S).grid(row=42, column=0,columnspan=6)

#boton Button
Button(f5, text='Create input', command=proc1, width=15,height=1,bg='white',font=("Helvetica",11,"bold")).grid(row=36, column=0)
Button(f5, text='Quit', command=proc5,width=15,height=1,bg='white',font=("Helvetica",11,"bold")).grid(row=36, column=1)

#ejecucin de ventana
ventana.mainloop()
