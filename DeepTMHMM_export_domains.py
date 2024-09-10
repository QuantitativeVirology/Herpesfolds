# -*- coding: utf-8 -*-
"""
input is the folder outputs from DeepTMHMM
    the outputs for each protein should be in a single folder

this script will copy the .gff3 files from each folder into one folder and rename the file to include the protein name
the topology will be read for all proteins and saved to a single file
output files:
    DeepTMHMM_all.tsv - all proteins with any annotation
    DeepTMHMM_withsignal.tsv - annotations for proteins with a signal peptide
    DeepTMHMM_withTM.tsv - annotations for proteins with a transmembrane domain
    DeepTMHMM_withsignalandTM.tsv - annotations for proteins with both a signal peptide and a transmembrane domain
    DeepTMHMM_withsignalandnoTM.tsv - annotations for proteins with a signal peptide but NOT a transmembrane domain
"""

import os
import shutil
import numpy as np

dirparent = r'parent_dir'
dirfolders = os.path.join(dirparent, 'folder_of_protein_folders')
dirsave = os.path.join(dirparent, 'gff3_files')
if not os.path.exists(dirsave):
    os.mkdir(dirsave)


# Functions
def Parsegff3file(_path, _protein_name):
    _domains = []
    with open(_path, "r") as file:
        _gff_contents = file.read()
    
    _flagsignal = 'signal' in _gff_contents
    _flagTM = 'TMhelix' in _gff_contents
    _gff_contents_split = _gff_contents.split('\t')
    
    for _i in range(0, round((len(_gff_contents_split)-1)/7)): #each domain has 7 '\t'. There is also a header. loop through each domain
        _inew = 1 + 7*_i
        # save protein name, domain type, domain start, domain end
        _new_domain = [_protein_name, _gff_contents_split[_inew], _gff_contents_split[_inew+1], _gff_contents_split[_inew+2]]
        _domains.append(_new_domain)
    return _domains, _flagsignal, _flagTM


# Main
listfolders = os.listdir(dirfolders)
domainsall = []
domainssignal = []
domainsTM = []
domainssecreted = []
domainsboth = []
for i in range(0, len(listfolders)):
    pathcopy = os.path.join(dirfolders, listfolders[i], 'TMRs.gff3')
    pathsave = os.path.join(dirsave, f'{listfolders[i]}_TMRs.gff3')
    if os.path.isfile(pathcopy):
        newdomains, flagsignal, flagTM = Parsegff3file(pathcopy, listfolders[i])
        domainsall.extend(newdomains)
        if flagsignal:
            domainssignal.extend(newdomains)
        if flagTM:
            domainsTM.extend(newdomains)
            
        if flagsignal and flagTM:
            domainsboth.extend(newdomains)
        elif flagsignal and not flagTM:
            domainssecreted.extend(newdomains)
            
        if not os.path.isfile(pathsave):
            shutil.copy2(pathcopy, pathsave)
    else:
        print(f'file not found: {listfolders[i]}')

np.savetxt(os.path.join(dirparent, 'DeepTMHMM_all.tsv'), np.asarray(domainsall), delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dirparent, 'DeepTMHMM_withsignal.tsv'), np.asarray(domainssignal), delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dirparent, 'DeepTMHMM_withTM.tsv'), np.asarray(domainsTM), delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dirparent, 'DeepTMHMM_withsignalandTM.tsv'), np.asarray(domainsboth), delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dirparent, 'DeepTMHMM_withsignalandnoTM.tsv'), np.asarray(domainssecreted), delimiter="\t", fmt="%5s")

print('done')