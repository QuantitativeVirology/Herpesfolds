# -*- coding: utf-8 -*-
"""
input is folder of Alphafold structures to truncate
output is extracted structures as .pdb
"""

# dependencies
import os
import shutil
import numpy as np
import prody
from Bio import PDB
parser = PDB.PDBParser()

THRESHPLDDT = 70
THRESHDOMAIN = 100

dirlistsave = r'C:\lists_save_dir'

dirparent = r'C:\parent_dir'
folder = 'folder_name'
dirpdb = os.path.join(dirparent, folder)

foldersave = 'truncated_proteins'
dirsave = os.path.join(dirparent, foldersave)
if os.path.exists(dirsave)==False:
    os.mkdir(dirsave)
    
## Functions
def Segment(_array, _value):
    _start = 0
    _current = _array[_start]
    _domains = []
    for _i in range(1, len(_array)):
        if not _array[_i] == _current or _i == len(_array) - 1:
            if _current == _value:
                _domains.append([_start, _i])
                _current = _array[_i]
            else:
                _start = _i
                _current = _array[_i]
    return _domains
    
## Main
listpdb = os.listdir(dirpdb)
listFL = []
listcut = []
listremove = []
for i in range(0, len(listpdb)):
    if i%10==0:
        print(f'currently at count: {i} of {len(listpdb)}, file: {listpdb[i]}\n')
    
    # load pdb
    structure = parser.get_structure('listdir[i]', os.path.join(dirpdb, listpdb[i]))
    model = structure[0]
    chain = model['A']

    # extract b-factor/pLDDT
    countresidue = 0
    pLDDT = np.empty(len(chain), dtype='object')
    for residue in chain:
        countatom = 0
        residueLDDT = np.empty(len(residue), dtype='object')
        for atom in residue:
            residueLDDT[countatom] = atom.get_bfactor()
            countatom += 1
        pLDDT[countresidue] = np.mean(residueLDDT)
        countresidue += 1
        
    # threshold pLDDT. identify stretches of low pLDDT. if smaller than THRESHDOMAIN, convert to passed
    pLDDTthresh = pLDDT > THRESHPLDDT
    pLDDTthreshsmooth = np.copy(pLDDTthresh)
    pLDDTdomainslow = Segment(pLDDTthresh, False)
    
    # if disordered region is small, merge it with the structured regions
    for j in range(0, len(pLDDTdomainslow)):
        if pLDDTdomainslow[j][1]-pLDDTdomainslow[j][0] < THRESHDOMAIN:
            for k in range(pLDDTdomainslow[j][0], pLDDTdomainslow[j][1]):
                pLDDTthreshsmooth[k] = True
                
    # identify structured regions
    pLDDTdomains = Segment(pLDDTthreshsmooth, True)
    
    # save cut protein
    if len(pLDDTdomains) == 1: # no cutting required
        listFL.append(listpdb[i][:-4])
        shutil.copy2(os.path.join(dirpdb, listpdb[i]), os.path.join(dirsave, f'{listpdb[i][:-4]}_domain-0.pdb'))
    elif len(pLDDTdomains) == 0: # all disordered
        listremove.append(listpdb[i][:-4])
    else:
        listcut.append([listpdb[i][:-4], str(pLDDTdomains)[1:-1]])
        for j in range(0, len(pLDDTdomains)):
            modeltocut = prody.parsePDB(os.path.join(dirpdb, listpdb[i]), model=1)
            extract = modeltocut.select(f'resnum {pLDDTdomains[j][0]}to{pLDDTdomains[j][1]}')
            prody.writePDB(os.path.join(dirsave, f'{listpdb[i][:-4]}_domain-{j}.pdb'), extract)
            
# save
np.savetxt(os.path.join(dirlistsave, 'PDB_domains_list_FL.tsv'), listFL, fmt="%5s", delimiter="\t")
np.savetxt(os.path.join(dirlistsave, 'PDB_domains_list_cut.tsv'), listcut, fmt="%5s", delimiter="\t")
np.savetxt(os.path.join(dirlistsave, 'PDB_domains_list_removed.tsv'), listremove, fmt="%5s", delimiter="\t")
print('done')