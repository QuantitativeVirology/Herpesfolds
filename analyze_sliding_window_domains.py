# -*- coding: utf-8 -*-
"""
input is: <file>_foldseek_reciprocal_hits_noself_ranges.tsv
    from from code: Sliding_Window_Analysis
all of the matches for a query protein are combined (if overlap)
output is list of query proteins with the identified domains
"""

# dependencies
import os
import numpy as np
import copy

# files and paths
dirparent = 'C:\\parentdir'
fileinput = 'file_foldseek_reciprocal_hits_noself_ranges.tsv'
pathinput = os.path.join(dirparent, fileinput)
listinput = np.genfromtxt(pathinput, delimiter="\t", dtype='str') #q-Virus_Protein, q-length, q-range from, q-range to, t-Virus_Protein, t-length, t-range from, t-range to

listinputsorted = np.asarray(sorted(listinput[1:, :], key=lambda x: x[4]))
listinputsorted = np.asarray(sorted(listinputsorted, key=lambda x: x[0]))

# constants
OVERLAPBUF = 41 # 1 step is 20


## Functions
def CastInt(_input):
    _output = np.empty(len(_input))
    for _i in range(0, len(_input)):
        _output[_i] = int(_input[_i])
    return _output

def Overlap(_array1, _array2):
    _array1 = CastInt(_array1)
    _array2 = CastInt(_array2)
    if _array1[0] > _array2[1] or _array1[1] < _array2[0]:
        return False
    else:
        return True
    
def Overlapbuffered(_array1, _array2):
    _array1 = CastInt(_array1)
    _array2 = CastInt(_array2)
    if _array1[0] > _array2[1] or _array1[1] < _array2[0]:
        return False
    else:
        _dif1 = (_array2[1] - _array1[0])
        _dif2 = (_array1[1] - _array2[0])
        _mindif = np.min([_dif1, _dif2])
        if _mindif < OVERLAPBUF:
            return False
        else:
            return True
        
        
## Main
listoutput = np.empty([np.shape(listinputsorted)[0], 6], dtype=object) # name, domains, self hit?, duplication?, addition?, summary flag
count = 0
currentprot = listinputsorted[0, 0]
start = 0
for i in range(0, np.shape(listinputsorted)[0]): # don't need to check final entry (there is nothing after it)
    if i % 1000 == 0:
        print(f'at {i} of {len(listinputsorted)}')
    
    endfound = False
    if i < np.shape(listinputsorted)[0]-1:
        if not listinputsorted[i+1, 0] == currentprot: # when start and end of a given protein have been found
            end = i
            endfound = True
    else:
        end = i
        endfound = True
        
    if endfound:        
        # create an array where if there was a match found at that residue the value is 1
        proteinbinary = np.zeros(int(np.max(CastInt(listinputsorted[start:end+1, 3])))+1)
        for j in range(start, end+1):
            isoform = False
            if listinputsorted[j, 0] in listinputsorted[j, 4]:  # if A in B
                isoform = True
            elif listinputsorted[j, 4] in listinputsorted[j, 0]: # if B in A
                isoform = True
            if listinputsorted[j, 0] in listinputsorted[j, 4]: # if they're actually the same
                isoform = False
            
            if not isoform:
                for k in range(int(listinputsorted[j, 2]), int(listinputsorted[j, 3])):
                    proteinbinary[k] = 1
        
        # split this array into domains
        domains = []
        newdomain = {}
        newdomain['start'] = 1
        newdomain['end'] = -1
        for j in range(1, len(proteinbinary)):
            if proteinbinary[j] == 1:
                if proteinbinary[j-1] == 0:
                    newdomain['start'] = j
            else:
                if proteinbinary[j-1] == 1:
                    newdomain['end'] = j
                    domains.append(copy.deepcopy(newdomain))
        if newdomain['end'] == -1:
            newdomain['end'] = len(proteinbinary)-1
            domains.append(copy.deepcopy(newdomain))
        
        domainlist = ''
        for j in range(0, len(domains)):
            domainlist = f'{domainlist},{domains[j]["start"]}-{domains[j]["end"]}'
        
        # check if domain duplication
        # the query protein matches itself or the target protein matches multiple times
        dupself = False
        dupnonself = False
        listoutput[count, 2] = ''
        listoutput[count, 3] = ''
        for j in range(start, end):
            if listinputsorted[j, 0] == listinputsorted[j, 4]:
                listoutput[count, 2] = 'self hit'
                
            if j < end:
                if listinputsorted[j, 4] == listinputsorted[j+1, 4]:
                    if Overlapbuffered(listinputsorted[j, 6:8], listinputsorted[j+1, 6:8]): # overlap means target domain is the same
                        if not Overlap(listinputsorted[j, 2:4], listinputsorted[j+1, 2:4]): # no overlap means query domain is different
                            if not listinputsorted[j, 0] == listinputsorted[j, 4]:
                                listoutput[count, 3] = 'duplication'
        
        # check if domain addition
        # not all found domains in the query are found in any same target protein (e.g. query has 2 domains and target only matches 1)
        hasaddition = False
        if len(domains) > 1:
            templisttarget = listinputsorted[start:end+1, 4]
            templisttargetunique = np.unique(templisttarget)
            for j in range(0, len(templisttargetunique)):
                tempcount = np.count_nonzero(templisttarget == templisttargetunique[j])
                if tempcount < len(domains): # if target has feweer occcurances, then it cannot have all of the domains
                    hasaddition = True
                else:
                    for k in range(0, len(domains)): # for every domain check
                        domainfound = False
                        for m in range(start, end+1): # is it found in the current target protein
                            if templisttargetunique[j] == listinputsorted[m, 4]:
                                if Overlap([domains[k]['start'], domains[k]['end']], listinputsorted[m, 6:8]):
                                    domainfound = True
                        if not domainfound:
                            hasaddition = True
                            break
        if hasaddition:
            listoutput[count, 4] = 'addition'
        else:
            listoutput[count, 4] = '' 
        
        outcome = ''
        if listoutput[count, 2]:
            outcome += f', {listoutput[count, 2]}'
        if listoutput[count, 3]:
            outcome += f', {listoutput[count, 3]}'
        if listoutput[count, 4]:
            outcome += f', {listoutput[count, 4]}'         
        
        listoutput[count, 0] = listinputsorted[i, 0]
        listoutput[count, 1] = domainlist[1:]
        listoutput[count, 5] = outcome[2:]
        count = count+1
        start = i+1
        if i < np.shape(listinputsorted)[0]-1:
            currentprot = listinputsorted[i+1, 0]

# save        
listoutputsave = np.delete(listoutput, range(count, np.shape(listoutput)[0]), 0)
np.savetxt(os.path.join(dirparent, f'{fileinput[:-4]}_merge_domains_annotated.tsv'), listoutputsave, fmt="%5s", delimiter="\t")
print('Done')
