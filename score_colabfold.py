# -*- coding: utf-8 -*- 
"""
input is ColabFold output folder
exports scores
"""   

# dependencies
import numpy as np
import pandas as pd
import json
import os
from collections import defaultdict

# files and paths
dirparent = 'C:\\parentdir'
folder = 'predictions' #folder is in above directory
path = os.path.join(dirparent, folder)
filepath = os.path.join(dirparent, f'Scores_{folder}.xlsx')

# constants
numoutput = 8
NUMSEED = 1


## Functions

def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record

def read_pdb(pdbfile):
    '''Read a pdb file predicted with AF and rewritten to contain all chains
    '''

    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            #Get CB - CA for GLY
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                    chain_plddt[record['chain']] = [record['B']]

    #Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt


## Main
# lists to store files
listdir = os.listdir(path)
respdb = np.empty(len(listdir), dtype='object')
resjson = np.empty(len(listdir), dtype='object')

# make list of pdb and json files in dir
countpdb = 0
countjson = 0
for file in os.listdir(path):
    if file.endswith((".pdb")):
        respdb[countpdb] = file
        countpdb += 1
    elif file.endswith((".json")) and ("rank" in file):
        resjson[countjson] = file
        countjson += 1
respdb = np.delete(respdb, range(countpdb, len(respdb)), 0)
resjson = np.delete(resjson, range(countjson, len(resjson)), 0)
respdbmissing = respdb.copy()
resjson.sort()

# extract scores
listfailed = []
listpassed = []
scores = np.zeros(((len(resjson)), numoutput), dtype=object)
print("files found:  ", len(resjson))
for i in range(0,len(resjson)):  
  # write file name into scores list
  scores[i,0]=resjson[i][0:-5]

  # Reading metadatafile for a single prediction model
  f = open(os.path.join(path,resjson[i]))
  metadata = json.load(f)
  plddt=np.array(metadata['plddt'])

  # find pdb and json files associated with the current json file
  # files are named: <file>_scores_rank_001_alphafold2_ptm_model_1_seed_000.json
  #                  <file>_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
  delim1 = [k for k,m in enumerate(resjson[i]) if m=='_']
  proteincombnamepdb = resjson[i][0:delim1[len(delim1)-9]+1]+"u"
  proteincombnamejson = resjson[i][0:delim1[len(delim1)-9]+1]+"s"
  proteincombmodel = resjson[i][delim1[len(delim1)-4]:-5]

  # pdb files
  proteincombfilespdb = []
  for j in range(0,len(respdb)):
      if proteincombnamepdb in respdb[j]:
          proteincombfilespdb.append(respdb[j])
  if not len(proteincombfilespdb) == 5*NUMSEED:
      print(f'Problem: {resjson[i]} \n{proteincombnamepdb}')
  proteincombfilespdb.sort()
    
  for j in range(0,len(proteincombfilespdb)):
      if proteincombmodel in proteincombfilespdb[j]:
          i_pdb=proteincombfilespdb[j]     
  i_pdb_pos = np.where(respdbmissing == i_pdb)[0]
  respdbmissing = np.delete(respdbmissing, i_pdb_pos, 0)

  # json files 
  proteincombfilesjson = []
  for j in range(0,len(resjson)):
      if proteincombnamejson in resjson[j]:
          proteincombfilesjson.append(resjson[j])
  proteincombfilesjson.sort()

  # read pdb
  pdbpath=os.path.join(path,i_pdb)
  chain_coords, chain_plddt = read_pdb(pdbpath)

  # Read chain lengths out of pdb file
  chainlength=np.empty(len(chain_coords), dtype="int")
  numchains=0
  numresidues=0
  for j in chain_coords:
      chainlength[numchains]=chain_coords[j].shape[0]
      numresidues=numresidues+chainlength[numchains]
      numchains=numchains+1

  # pTM
  ptm=np.array(metadata['ptm'])
  scores[i,1]=np.round(np.max(ptm),2)
  
  # pLDDT %>70% and %<50%
  counthigh=0
  countlow=0
  for j in range(0,len(plddt)):
      if plddt[j]>70:
          counthigh=counthigh+1
      elif plddt[j]<50:
          countlow=countlow+1
  scores[i,2]=counthigh/len(plddt)
  scores[i,3]=countlow/len(plddt)
  scores[i,4]=1-scores[i,3]-scores[i,4]

  # pLDDT quality scores
  listplddtstats = np.zeros((2, len(proteincombfilesjson)+1), dtype=object)
  for j in range(0, len(proteincombfilesjson)):
      # Reading metadatafile for a single prediction model
      f = open(os.path.join(path,proteincombfilesjson[j]))
      metadata = json.load(f)
      listplddtstats[0,j] = np.mean(np.array(metadata['plddt']))
      listplddtstats[1,j] = np.std(np.array(metadata['plddt']))
  scores[i,5] = np.std(listplddtstats[0])
  scores[i,6] = np.std(listplddtstats[1])
  if scores[i,5] > 3.215218 or scores[i,6] > 1.869734 or scores[i,1] < 0.315: # threshold is mean+2x stdev of fit gaussian
      scores[i,7] = 'fail'
      listfailed.append(scores[i,0])
  else:
      scores[i,7] = 'pass'
      listpassed.append(scores[i,0])

# save files
# model quality lists
np.savetxt(os.path.join(dirparent, f'{folder}_list_colabfold_passed.tsv'), listpassed, delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dirparent, f'{folder}_list_colabfold_failed.tsv'), listfailed, delimiter="\t", fmt="%5s")

# save scores as excel file
scoresDF = pd.DataFrame(scores)
with pd.ExcelWriter(filepath) as writer:
    # the value in [] is the thresholds
    scoresDF.to_excel(writer, header=["file", "pTM [>0.315]", "pLDDT %>70% [>0.31]", "pLDDT %<50% [>0.47]", "pLDDT % 50-70%", \
                                      "stdev mean pLDDT [<3.2]", "stdev stdev pLDDT [<1.9]", "model quality?"])
print('Done')
