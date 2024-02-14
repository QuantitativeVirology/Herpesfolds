# -*- coding: utf-8 -*-
"""
export all signicant reciprocal hits
remove only self hits
        
input is full length protein foldseek output file
output is:
    <file>_foldseek_reciprocal_hits.tsv
    <file>_foldseek_pairs.tsv
    <file>_foldseek_groups_compare_homologs.tsv
"""

# dependencies
import os
import numpy as np

# file paths
dirparent = 'C:\\parentdir'

#input file column order: query[0],target[1],theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue[11],bits
fileinput = 'file.txt'
pathinput = os.path.join(dirparent, fileinput)
listinput = np.genfromtxt(pathinput, delimiter='\t', dtype='str')

pathhomolog = 'C:\\chartdir\\HerpesFolds.tsv' #homolog chart
listinputhomolog = np.genfromtxt(pathhomolog, delimiter="\t", dtype='str')

pathsave = os.path.join(dirparent, f'{fileinput[:-4]}_hits.tsv')
filelistoutput = os.path.join(dirparent, 'foldseek_FL_FL_reciprocal_hits.tsv')

# constants
SIG = 0.001


### Functions
def Deldup(_listsorted, _i):
    if _i < len(_listsorted)-1:
        if _listsorted[_i, 0] == _listsorted[_i+1, 0] and _listsorted[_i, 1]==_listsorted[_i+1, 1]:
            _listsorted = np.delete(_listsorted, _i+1, 0)
            _listsorted = Deldup(_listsorted, _i)
    return _listsorted

def FindGroup(start, _listgroupstart):
    _output = [start]

    _indexstart = np.where(_listgroupstart == start) #y,x
    _indexstartsize = np.shape(_indexstart)[1]

    for _i in range(0, _indexstartsize):
      _listgroupstart[_indexstart[0][_i], _indexstart[1][_i]] = -1

    for _i in range(0, _indexstartsize):
      _partner = _listgroupstart[_indexstart[0][_i], 1-_indexstart[1][_i]]
      if not (_partner == '-1'):
        _next = FindGroup(_partner, _listgroupstart)
        _output = _output + _next
    return _output

def RemoveRedundant(_listoutputpair):
    # remove duplicates
    _listoutputpairsorted = np.asarray(sorted(_listoutputpair, key=lambda x: x[1]))
    _listoutputpairsorted = np.asarray(sorted(_listoutputpairsorted, key=lambda x: x[0]))

    _listoutputpairsorteddel = _listoutputpairsorted
    for _i in range(0, len(_listoutputpairsorteddel)):
        if _i < len(_listoutputpairsorteddel)-1:
            _listoutputpairsorteddel = Deldup(_listoutputpairsorteddel, _i)
            
    # remove reciprical
    _listoutputpairsorteddelrecip = _listoutputpairsorteddel
    for _i in range(0, len(_listoutputpairsorteddelrecip)):
        if _i < len(_listoutputpairsorteddelrecip)-1:
            _indexrecip = np.where(_listoutputpairsorteddelrecip[:, 1] == _listoutputpairsorteddelrecip[_i, 0])[0]
            if not len(_indexrecip) == 0:
                for _j in _indexrecip:
                    if _j < np.shape(_listoutputpairsorteddelrecip)[0]:
                        if _listoutputpairsorteddelrecip[_j, 0] == _listoutputpairsorteddelrecip[_i, 1]:
                            _listoutputpairsorteddelrecip = np.delete(_listoutputpairsorteddelrecip, _j, 0)
    return _listoutputpairsorteddelrecip


## Main
# find pairs that are significant both ways - FL or domain is significant against any FL or domain
print('start finding reciprocal pairs')
listoutput = [[0,0,0,0,0,0,0,0,0,0,0,0,0]]
listoutput = np.delete(listoutput,0,0)
for i in range(0, len(listinput)):
    if i % 10000 == 0:
        print(f'at {i} of {len(listinput)}')
    
    if float(listinput[i, 11]) < SIG and not 'domain' in listinput[i, 0]:
        queryname = listinput[i, 0]
        indexrecip = np.where(listinput[:, 1] == queryname)[0]
        
        if not len(indexrecip) == 0:
            recipsig = False
            for j in indexrecip:
                targetname = listinput[j, 0]
                if targetname == listinput[i,1] and not 'domain' in listinput[j, 0]:
                    if float(listinput[j, 11]) < SIG:
                        recipsig = True
                        break
            if recipsig == True:
                listoutput = np.append(listoutput, [listinput[i,:]], 0)

# make new list of pairs without self matches
print('start making pair list')
listoutputpair = [[0,0,0,0]]
listoutputpair = np.delete(listoutputpair,0,0)
for i in range(0, len(listoutput)):       
    queryname = listoutput[i, 0]
    targetname = listoutput[i, 2]
    #list of protein against protein
    if not queryname == targetname:
        sig_FL = []
        indextemp = [k for k, m in enumerate(listoutput[:, 0]) if queryname in m]
        for j in indextemp:
            if targetname in listoutput[j, 1]:
                sig_FL.append(float(listoutput[j, 11]))
        indextemp2 = [k for k, m in enumerate(listoutput[:, 1]) if queryname in m]
        for j in indextemp2:
            if targetname in listoutput[j, 0]:
                sig_FL.append(float(listoutput[j, 11]))   
        listoutputpair = np.append(listoutputpair, [[queryname[:-4], targetname[:-4], np.mean(sig_FL), np.std(sig_FL)]], 0)
listoutputpairsorteddelrecip = RemoveRedundant(listoutputpair)

# find structural groups
print('start finding structural groups')
listgroupstart = np.copy(listoutputpairsorteddelrecip)
listgroups = []
listgroupssave = []
for i in range(0, len(listgroupstart)):
  if not listgroupstart[i,0] == '-1':
    newgroup = []
    add = FindGroup(listgroupstart[i, 0], listgroupstart)
    if not add == '-1':
      newgroup = newgroup + add
    listgroups.append(newgroup)

# compare groups to homologs
print('start network comparison')
listhomologgroups = np.empty([len(listgroups), 3], dtype = object)
for i in range(0, len(listgroups)):
    listhomologgroups[i, 0] = ''.join([f'{str(element)}, ' for element in listgroups[i]])
    allpos = []
    for j in range(0,len(listgroups[i])):
        index = np.where(listinputhomolog == listgroups[i][j])[0]
        if index.size < 1:
            index = [-1]
        allpos = np.append(allpos, [str(index[0])], axis=0)
    listhomologgroups[i, 1] = ''.join([f'{str(element)}, ' for element in allpos])
    listhomologgroups[i, 2] = len(allpos)
listhomologgroupssorted = np.asarray(sorted(listhomologgroups, key=lambda x: x[2], reverse=True))

# save
np.savetxt(os.path.join(dirparent, 'foldseek_FL_FL_reciprocal_hits.tsv'), listoutput, fmt="%5s", delimiter="\t")
np.savetxt(os.path.join(dirparent, 'foldseek_FL_FL_pairs.tsv'), listoutputpairsorteddelrecip, fmt="%5s", delimiter="\t") #query target mean stdev
np.savetxt(os.path.join(dirparent, 'foldseek_FL_FL_groups_compare_homologs.tsv'), listhomologgroupssorted, fmt="%5s", delimiter="\t")
print("Done")
