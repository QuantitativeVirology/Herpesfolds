# -*- coding: utf-8 -*-
"""
input is DeepTMHMM_withsignal.tsv text output from DeepTMHMM_export_domains.py (put it in the parent_dir)
    requires path to folder of fasta files of full length proteins

output is fasta files without signal peptide
    will be in the parent_dir in the created folder 'monomers_fasta-nosignal'
"""

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import os

dirsource = r'monomers_fasta_folder'
dirparent = r'parent_dir'
dirsave = os.path.join(dirparent, 'monomers_fasta-nosignal')
pathinput = os.path.join(dirparent, 'DeepTMHMM_withsignal.tsv')

if not os.path.exists(dirsave):
    os.mkdir(dirsave)

#Main
listinput = np.genfromtxt(pathinput, delimiter='\t', dtype='str')

listproteins = []
for i in range(0, len(listinput)):
    if listinput[i, 1] == 'signal':
        #find protein size
        for j in range(i+1, len(listinput)):
            if listinput[j, 1] == 'signal':
                end = int(listinput[j-1, 3])
                break
        if j == len(listinput)-1:
            end = int(listinput[j, 3])
        
        pathsource = os.path.join(dirsource, f'{listinput[i, 0]}.fasta')
        if os.path.isfile(pathsource):
            for record in SeqIO.parse(pathsource, "fasta"):
                newprotein = [listinput[i, 0], record.id, str(record.seq), str(record.seq)[int(listinput[i, 3]):]]
                listproteins.append(newprotein)
                
                #check that the protein is the right length in case something was switched
                if not end == len(newprotein[2]):
                    print(f'{listinput[i, 0]} not the right length')
                
                record.seq = Seq(newprotein[3])
                SeqIO.write(record, os.path.join(dirsave, f'{listinput[i, 0]}-nosignal.fasta'), "fasta")
        else:
            print(f'missing fasta: {listinput[i, 0]}')
            
np.savetxt(os.path.join(dirparent, 'monomers_fasta-nosignal.tsv'), np.asarray(listproteins), delimiter="\t", fmt="%5s")
print('done')