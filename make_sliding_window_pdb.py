# -*- coding: utf-8 -*-
"""   
For every protein, cut it into pieces that are {sizewindow} in length and
{sizestep} spacing along protein
For remaining residues at end, add to last piece

file name convention:
    virus_protein_proteinlength_#pieces_residuestart_residueend
    
input is folder of rank 1 pdb models that passed quality scores
output is sliding window pieces as pdb
"""

# dependencies
import os
import math
import prody
from Bio import SeqIO
import numpy as np

# files and paths
dirparent = 'C:\\parentdir'
folder = 'predictions_passed_quality'  #folder is in above directory
dirpdb = os.path.join(dirparent, folder)
foldersave = 'protein_sliding_window'

# constants
# median herpes protein length is ~300 residues, second peak at 95
SIZEWINDOW = [150] 
SIZESTEP = [20]


## Functions
def padnumber(_value, _numdigits):
    if len(str(_value)) < _numdigits:
        dif = _numdigits - len(str(_value))
        countname=''
        for j in range(0,dif):
            countname = countname +'0'
        countname = countname + str(_value)
    else:
        countname=str(_value)
    return countname

def exportpieces(_filename, _lenprotein, _sizewindow, _sizestep, _pathpdb, _dirsavepieces):
    _allpieces = [[0, 0, 0]]
    _allpieces = np.delete(_allpieces, 0, 0)
    
    numpieces = math.floor((_lenprotein - _sizewindow) / _sizestep) + 1
    _model = prody.parsePDB(_pathpdb, model=1)
    for _i in range(0, numpieces):
        _start = _i*_sizestep + 1
        _end = _start + _sizewindow - 1
        if _i == numpieces - 1: #if the last piece, extend window to the end
            _end = _lenprotein
        
        filename = f'{_filename}_{_lenprotein}_{numpieces}_{padnumber(_start, len(str(_lenprotein)))}_{padnumber(_end, len(str(_lenprotein)))}.pdb'

        _extract = _model.select(f'resnum {_start}to{_end}')
        prody.writePDB(os.path.join(_dirsavepieces, filename), _extract)
        
        _newpiece = [[_filename, _start, _end]]
        _allpieces = np.append(_allpieces, _newpiece, 0)

    return _allpieces #return list of pieces as np 2D array


## Main
# make list of pdb files in dir
respdb = []
for file in os.listdir(dirpdb):
    if file.endswith((".pdb")):
        respdb.append(file)
respdb.sort()

# create save folders and empty array
if os.path.exists(os.path.join(dirparent, foldersave))==False:
    os.mkdir(os.path.join(dirparent, foldersave))
    
listoutput = {}
for i in range(0, len(SIZEWINDOW)):
    listoutput[i] = [[0, 0, 0]]
    listoutput[i] = np.delete(listoutput[i], 0, 0)
    foldersavepieces = f'Protein_pieces_window{padnumber(SIZEWINDOW[i], len(str(SIZEWINDOW[len(SIZEWINDOW)-1])))}_step{SIZESTEP[i]}'
    dirsavepieces = os.path.join(dirparent, foldersave, foldersavepieces)
    if os.path.exists(dirsavepieces)==False:
        os.mkdir(dirsavepieces)

# interate for each protein, save each split
for file in respdb:
    pathpdb = os.path.join(dirpdb, file)
    for record in SeqIO.parse(pathpdb, "pdb-atom"):
        lenprotein = len(record.seq._data)
        for i in range(0, len(SIZEWINDOW)):
            foldersavepieces = f'Protein_pieces_window{padnumber(SIZEWINDOW[i], len(str(SIZEWINDOW[len(SIZEWINDOW)-1])))}_step{SIZESTEP[i]}'
            dirsavepieces = os.path.join(dirparent, foldersave, foldersavepieces)
            listpieces = exportpieces(file[:-4], lenprotein, SIZEWINDOW[i], SIZESTEP[i], pathpdb, dirsavepieces)
            listoutput[i] = np.append(listoutput[i], listpieces, 0)
            print(f'    Exported pieces {file} window {SIZEWINDOW[i]}, step {SIZESTEP[i]}')

# save
for i in range(0, len(SIZEWINDOW)):
    np.savetxt(os.path.join(dirparent, foldersave, f'Protein_pieces_window{padnumber(SIZEWINDOW[i], len(str(SIZEWINDOW[len(SIZEWINDOW)-1])))}_step{SIZESTEP[i]}.tsv'), listoutput[i], fmt="%5s", delimiter="\t")
print('Done')
