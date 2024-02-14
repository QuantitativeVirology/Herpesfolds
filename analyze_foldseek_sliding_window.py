# -*- coding: utf-8 -*-
"""
export all signicant reciprocal hits
remove only EXACT self hits
    if one part of the protein matches another, it is kept
        UNLESS the sliding window overlaps
        
input is sliding window pieces foldseek output file
output is:
    <file>_foldseek_reciprocal_hits.tsv
    <file>_foldseek_reciprocal_hits_noself.tsv
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
print(fileinput)
pathinput = os.path.join(dirparent, fileinput)

pathhomolog = 'C:\\chartdir\\HerpesFolds.tsv' #homolog chart
listinputhomolog = np.genfromtxt(pathhomolog, delimiter="\t", dtype='str')

pathsave = os.path.join(dirparent, f'{fileinput[:-4]}_hits.tsv')
filelistoutput = os.path.join(dirparent, f'{fileinput[:-4]}_foldseek_reciprocal_hits.tsv')

# constants
SIG = 0.001


## Functions
def FindGroup(_start):
    _output = [_start]

    _indexstart = np.where(listgroupstart == _start) #y,x
    _indexstartsize = np.shape(_indexstart)[1]

    for _i in range(0, _indexstartsize):
        listgroupstart[_indexstart[0][_i], _indexstart[1][_i]] = '-1'

    for _i in range(0, _indexstartsize):
        _partner = listgroupstart[_indexstart[0][_i], 1-_indexstart[1][_i]]
        if not (_partner == '-1'):
            _next = FindGroup(_partner)
            _output = _output + _next
    return _output

def SplitSearchIndexInSub(_value, _col, _start, _end):
    #_col: 0 is query and 1 is target in listinputsorted
    #returns position of _value; ASSUMES that _value is unique
    if _end - _start == 0 or _end - _start == 1:
        if _value == listinputpro[_start, _col]:
            _index = _start
        elif _value == listinputpro[_end, _col]:
            _index = _end
        else:
            _index = -1
    else:
        _middle = int(np.round(np.mean([_start, _end])))
        if _value == listinputpro[_middle, _col]:
            _index = _middle
        else:
            if _value > listinputpro[_middle, _col]:
                _index = SplitSearchIndexInSub(_value, _col, _middle + 1, _end)
            else:
                _index = SplitSearchIndexInSub(_value, _col, _start, _middle - 1)  
    return _index

def FindStart(_value, _col, _indexmiddle):  
    _indexstart = -1
    if _indexmiddle == 0:
        _indexstart = 0
    else:
        for _i in range(_indexmiddle, -1, -1):
            if not _value == listinputpro[_i, _col]:
                _indexstart = _i + 1
                break
        if _i == 0:
            _indexstart = 0
    return _indexstart

def FindEnd(_value, _col, _indexmiddle):
    _indexend = -1
    if _indexmiddle == len(listinputpro) - 1:
        _indexend = len(listinputpro) - 1
    else:
        for _i in range(_indexmiddle, len(listinputpro)):
            if not _value == listinputpro[_i, _col]:
                _indexend = _i - 1
                break
        if _i == len(listinputpro) - 1:
            _indexend = len(listinputpro) - 1
    return _indexend

def SplitSearchArrayInWhole(_value, _col, _start, _end):
    #_col: 0 is query and 1 is target in listinputsorted; search array is global variable: listinputsorted
    #returns array of positions
    if _end - _start == 0 or _end - _start == 1:
        if _value == listinputpro[_start, _col]:
            _indexstart = FindStart(_value, _col, _start)
            _indexend = FindEnd(_value, _col, _start)
            _index = [_indexstart, _indexend]
        elif _value == listinputpro[_end, _col]:
            _indexstart = FindStart(_value, _col, _end)
            _indexend = FindEnd(_value, _col, _end)
            _index = [_indexstart, _indexend]
        else:
            _index = -1
    else:
        _middle = int(np.round(np.mean([_start, _end])))
        if _value == listinputpro[_middle, _col]:
            _indexstart = FindStart(_value, _col, _middle)
            _indexend = FindEnd(_value, _col, _middle)
            _index = [_indexstart, _indexend]
        else:
            if _value > listinputpro[_middle, _col]:
                _index = SplitSearchArrayInWhole(_value, _col, _middle, _end)
            else:
                _index = SplitSearchArrayInWhole(_value, _col, _start, _middle)  
    return _index

def RecipSearch(_query, _target):
    #input list is global variable: listinputsorted (MUST be sorted in INCREASING order: a, b, ...)
    #returns True or False
    _indexqueryreciparray = SplitSearchArrayInWhole(_query, 1, 0, len(listinputpro) - 1)
    _recipsig = False
    if not _indexqueryreciparray == -1:
        _indextargetrecippos = SplitSearchIndexInSub(_target, 0, _indexqueryreciparray[0], _indexqueryreciparray[1])
        if not _indextargetrecippos == -1:
            if float(listinputpro[_indextargetrecippos, 11]) < SIG:
                _recipsig = True
    return _recipsig


## Main
# find pairs that are significant both ways
print('\nstart finding reciprocal pairs')
listinput = np.genfromtxt(pathinput, delimiter='\t', dtype='str')

#pre-process input list; take significant hits and sort
print(' Pre-processing list')
listinputpro = np.empty(np.shape(listinput), dtype='object')
count = 0
for i in range(0, len(listinput)):
    if float(listinput[i, 11]) < SIG:
        listinputpro[count] = listinput[i]
        count += 1
listinputpro = np.delete(listinputpro, range(count, np.shape(listinputpro)[0]), 0)
listinputpro = np.asarray(sorted(listinputpro, key=lambda x: x[0]))
listinputpro = np.asarray(sorted(listinputpro, key=lambda x: x[1]))

listoutput = np.empty(np.shape(listinputpro), dtype='object')
count = 0
for i in range(0, len(listinputpro)):
    recipsig = RecipSearch(listinputpro[i, 0], listinputpro[i, 1])
    if recipsig == True:
        temp = np.copy(listinputpro[i,:])
        temp[0] = temp[0][:-4]
        temp[1] = temp[1][:-4]
        temp[2] = temp[2][:-4]
        listoutput[count] = temp
        count += 1
listoutput = np.delete(listoutput, range(count, np.shape(listoutput)[0]), 0)
np.savetxt(filelistoutput, listoutput, fmt="%5s", delimiter="\t")

# make new list of pairs without self matches
print('\nstart making pair list')
listoutputpair = np.empty([np.shape(listoutput)[0], 2], dtype='object')
count = 0
listoutputnoself = np.empty(np.shape(listoutput), dtype='object')
countnoself = 0
for i in range(0, len(listoutput)):
    queryname = listoutput[i, 0]
    querydelim = [k for k,m in enumerate(queryname) if m=='_']
    querynametrunc = f'{queryname[:querydelim[len(querydelim)-4]]}'
    querystart = int(f'{queryname[querydelim[len(querydelim)-2]+1:querydelim[len(querydelim)-1]]}') 
    queryend = int(f'{queryname[querydelim[len(querydelim)-1]+1:]}') 

    targetname = listoutput[i, 1]
    targetdelim = [k for k,m in enumerate(targetname) if m=='_']
    targetnametrunc = f'{targetname[:targetdelim[len(targetdelim)-4]]}'
    targetstart = int(f'{targetname[targetdelim[len(targetdelim)-2]+1:targetdelim[len(targetdelim)-1]]}') 
    targetend = int(f'{targetname[targetdelim[len(targetdelim)-1]+1:]}') 
       
    save = False
    if not querynametrunc == targetnametrunc: # not a self match, keep
        save = True
    elif queryend < targetstart or targetend < querystart: # if it's a self match, keep if ranges don't overlap
        save = True
    if save == True:
        listoutputpair[count] = [querynametrunc, targetnametrunc]
        count += 1
        listoutputnoself[countnoself] = listoutput[i, :]
        countnoself += 1
listoutputpair = np.delete(listoutputpair, range(count, np.shape(listoutputpair)[0]), 0)
listoutputnoself = np.delete(listoutputnoself, range(countnoself, np.shape(listoutputnoself)[0]), 0)

# remove duplicates
print('\n remove duplicates')
listoutputpairsorted = np.asarray(sorted(listoutputpair, key=lambda x: x[1]))
listoutputpairsorted = np.asarray(sorted(listoutputpairsorted, key=lambda x: x[0]))

listoutputpairsorteddel = np.copy(listoutputpairsorted)
count = 0
for i in range(0, len(listoutputpairsorteddel)):
    if i < len(listoutputpairsorteddel)-1:
        if listoutputpairsorteddel[i, 0] == listoutputpairsorteddel[i + 1, 0] and listoutputpairsorteddel[i, 1] == listoutputpairsorteddel[i + 1, 1]:
            listoutputpairsorteddel[i, 0] = '0'
            count += 1
listoutputpairsorteddel = np.asarray(sorted(listoutputpairsorteddel, key=lambda x: x[0]))
if not count ==0:
    listoutputpairsorteddel = np.delete(listoutputpairsorteddel, range(0, count), 0)
        
# remove reciprical
print('\n remove reciprical')
listoutputpairsorteddelrecip = listoutputpairsorteddel
for i in range(0, len(listoutputpairsorteddelrecip)):
    if i < len(listoutputpairsorteddelrecip)-1:
        if not listoutputpairsorteddelrecip[i, 0] == listoutputpairsorteddelrecip[i, 1]: #if it's a self match, it is because it is to a different domain
            indexrecip = np.where(listoutputpairsorteddelrecip[:, 1] == listoutputpairsorteddelrecip[i, 0])[0]
            if not len(indexrecip) == 0:
                for j in indexrecip:
                    if j < np.shape(listoutputpairsorteddelrecip)[0]:
                        if listoutputpairsorteddelrecip[j, 0] == listoutputpairsorteddelrecip[i,1]:
                            listoutputpairsorteddelrecip = np.delete(listoutputpairsorteddelrecip, j, 0)

# find structural similarity groups
print('\nstart finding structural similarity groups')
listgroupstart = np.copy(listoutputpairsorteddelrecip)
listgroups = []
for i in range(0, len(listgroupstart)):
  if not listgroupstart[i,0] == '-1':
    newgroup = []
    add = FindGroup(listgroupstart[i, 0])
    if not add == '-1':
        newgroup = newgroup + add
    listgroups.append(newgroup)

# compare groups to homologs chart
print('\nstart network comparison')
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
np.savetxt(os.path.join(dirparent, f'{fileinput[:-4]}_foldseek_reciprocal_hits_noself.tsv'), listoutputnoself, fmt="%5s", delimiter="\t")
np.savetxt(os.path.join(dirparent, f'{fileinput[:-4]}_foldseek_pairs.tsv'), listoutputpairsorteddelrecip, fmt="%5s", delimiter="\t")
np.savetxt(os.path.join(dirparent, f'{fileinput[:-4]}_foldseek_groups_compare_homologs.tsv'), listhomologgroupssorted, fmt="%5s", delimiter="\t")
print("\nDone")