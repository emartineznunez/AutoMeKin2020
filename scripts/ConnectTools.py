import numpy as np
#from Main import dict_ref

atref1,atref2,dref,bondtype=np.genfromtxt("thdist",dtype="|U5",unpack=True)
dict_ref={}
for i in range(len(dref)):
    dict_ref[atref1[i]+atref2[i]]=float(dref[i])
    if atref1[i] != atref2[i]:
       dict_ref[atref2[i]+atref1[i]]=float(dref[i])


# Method to get the center of mass seperation for two fragments.
# Frag1 has the start and end index of the smaller fragment in xyz
def getCOMdist(mol, frag):
    mass1 = 0.0
    mass2 = 0.0
    masses = mol.get_masses()
    COM1 = np.zeros(3)
    COM2 = np.zeros(3)

    # First need total mass of each fragment
    for i in range(0,masses.size):
        if i >= frag[0] and i <= frag[1]:
            mass1 += masses[i]
        else:
            mass2 += masses[i]

    # Then determine center of mass co-ordinates
    for i in range(0,masses.size):
        if i >= frag[0] and i <= frag[1]:
            COM1[0] += masses[i] * mol.get_positions()[i,0]
            COM1[1] += masses[i] * mol.get_positions()[i,1]
            COM1[2] += masses[i] * mol.get_positions()[i,2]
        else:
            COM2[0] += masses[i] * mol.get_positions()[i,0]
            COM2[1] += masses[i] * mol.get_positions()[i,1]
            COM2[2] += masses[i] * mol.get_positions()[i,2]

    COM1 /= mass1
    COM2 /= mass2

    # Finally calculate the distance between COM1 and COM2
    COMdist = np.sqrt( ((COM1[0] - COM2[0]) ** 2) + ((COM1[1] - COM2[1]) ** 2) + ((COM1[2] - COM2[2]) ** 2))
    return COMdist

def getCOMonly(mol):
    mass = 0.0
    COM = np.zeros(3)
    masses = mol.get_masses()
    # First need total mass of each fragment
    for i in range(0,masses.size):
        mass += masses[i]

    # Then determine center of mass co-ordinates
    for i in range(0,masses.size):
        COM[0] += masses[i] * mol.get_positions[i,0]
        COM[1] += masses[i] * mol.get_positions[i,1]
        COM[2] += masses[i] * mol.get_positions[i,2]

    return COM

# Method to return derivative of COM seperation via chain rule
# Needs double CHECKING
def getCOMdel(Mol, frag):
    mass1 = 0.0
    mass2 = 0.0
    masses = Mol.get_masses()
    COM1 = np.zeros(3)
    COM2 = np.zeros(3)
    #First need total mass of each fragment
    for i in range(0,masses.size):
        if i >= frag[0] and i <= frag[1]:
            mass1 += masses[i]
        else :
            mass2 += masses[i]
    #Then determine center of mass co-ordinates
    for i in range(0,masses.size):
        if i >= frag[0] and i <= frag[1]:
            COM1[0] += masses[i] * Mol.get_positions()[i,0]
            COM1[1] += masses[i] * Mol.get_positions()[i,1]
            COM1[2] += masses[i] * Mol.get_positions()[i,2]
        else:
            COM2[0] += masses[i] * Mol.get_positions()[i,0]
            COM2[1] += masses[i] * Mol.get_positions()[i,1]
            COM2[2] += masses[i] * Mol.get_positions()[i,2]

    COM1 /= mass1
    COM2 /= mass2

    # Finally calculate the distance between COM1 and COM2
    COMdist = np.sqrt( ((COM1[0] - COM2[0]) ** 2) + ((COM1[1] - COM2[1]) ** 2) + ((COM1[2] - COM2[2]) ** 2))

    # Now need the derivative component wise
    constraint = np.zeros(Mol.get_positions().shape)
    for i in range(0,masses.size):
        for j in range(0,3):
            constraint[i][j] = 1 / ( 2 * COMdist)
            constraint[i][j] *= 2 * (COM1[j] - COM2[j])
            if i >= frag[0] and i <= frag[1]:
                constraint[i][j] *= -masses[i] / mass1
            else:
                constraint[i][j] *= masses[i] / mass2
    return constraint

# Set up a reference matrix for ideal bond length between any two atoms in the system
# Maps species types onto a grid of stored ideal bond distances stored in the global variables module
def refBonds(mol):
    size =len(mol.get_positions())
    symbols = mol.get_chemical_symbols()
    dRef = np.zeros((size,size))
    for i in range(0 ,size) :
        for j in range(0, size) :
            sp = symbols[i] + symbols[j]
            if i!=j:
               dRef[i][j] = dict_ref[sp]
            else:
               dRef[i][j] = 1
    return dRef

def bondMatrix(dRef,mol):
    size =len(mol.get_positions())
    C = np.zeros((size,size))
    dratio = np.zeros((size,size))
    for i in range(0,size):
        for j in range(0,size):
            C[i][j] = (mol.get_distance(i,j))
            if i != j:
                dratio[i][j] = C[i][j] / dRef[i][j]
            if i == j:
                C[i][j] = 2
            elif C[i][j] < dRef[i][j] * 1.2:
                C[i][j] = 1.0
            else:
                C[i][j] = 0.0
    return C

def getChangedBonds(mol1, mol2):
    r = refBonds(mol1)
    C1 = bondMatrix(r, mol1)
    C2 = bondMatrix(r, mol2)
    indicies = []
    size =len(mol1.get_positions())
    for i in range(1,size):
        for j in range(0,i):
            if C1[i][j] != C2[i][j]:
                indicies.append(i)
                indicies.append(j)
    ind2 = []
    [ind2.append(item) for item in indicies if item not in ind2]
    return ind2

def getChangedBonds2(mol1, mol2):
    r = refBonds(mol1)
    C1 = bondMatrix(r, mol1)
    C2 = bondMatrix(r, mol2)
    indicies = []
    size =len(mol1.get_positions())
    for i in range(1,size):
        for j in range(0,i):
            if C1[i][j] != C2[i][j]:
                indicies.append([i,j])
    ind2 = []
    [ind2.append(item) for item in indicies if item not in ind2]
    return ind2


def getDistMatrix(mol,active):
    if active == "all":
        s1 = len(mol.get_positions())
        s2 = s1*(s1+1)/2
    else:
        s2 = len(active)
    D = np.zeros((s2))
    Dind = []
    if active == "all":
        n = 0
        for i in range(0,s1):
            for j in range(0,(s1 - i)):
                Dist = mol.get_distance(i,j)
                D[n] = Dist
                Dind.append((i,j))
                n += 1
    else:
        for i in range(0,s2):
            D[i] = mol.get_distance(active[i][0],active[i][1])
            Dind.append([active[i][0],active[i][1]])
    return D,Dind

def projectPointOnPath(S,path,type,n,D,reac):
    baseline = S - reac
    if type == 'linear':
        project = np.vdot(baseline,path) / np.linalg.norm(path)
    if type == 'distance':
        project = np.linalg.norm(baseline)
    Sdist = np.vdot(S,n) + D
    return Sdist,project

def genBXDDel(mol,D,Dind,n):
    constraint = np.zeros(mol.get_positions().shape)
    pos = mol.get_positions()
    for i in range(0,len(mol)):
        for j in range(0,len(Dind)):
            if D[j] != 0:
                firstTerm = (1/(2*D[j]))
            else:
                firstTerm=0
            if Dind[j][0] == i:
                constraint[i][0] += firstTerm*2*(pos[i][0]-pos[Dind[j][1]][0])*n[j]
                constraint[i][1] += firstTerm*2*(pos[i][1]-pos[Dind[j][1]][1])*n[j]
                constraint[i][2] += firstTerm*2*(pos[i][2]-pos[Dind[j][1]][2])*n[j]
            if Dind[j][1] == i:
                constraint[i][0] += firstTerm*2*(pos[Dind[j][0]][0]-pos[i][0])*-1*n[j]
                constraint[i][1] += firstTerm*2*(pos[Dind[j][0]][1]-pos[i][1])*-1*n[j]
                constraint[i][2] += firstTerm*2*(pos[Dind[j][0]][2]-pos[i][2])*-1*n[j]
    return constraint
