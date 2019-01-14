'''
Created on Sep 28, 2017

@author: vasilis verroios
'''

import time
import sys

from core.hashFunctions import titleANDAuthorsSim, shingling, DBLPprocessV4,\
    titleAuthorsSimANDothers
from random import shuffle, sample
from itertools import combinations
from collections import defaultdict
from lsh.utilities import loadDBLPcomplete
from core.metrics import f1, f1Fast

class titleAuthorsSimpleANDRule:
    def __init__(self, thres1, thres2, baseRule=titleANDAuthorsSim):
        self._thres1 = thres1
        self._thres2 = thres2
        self._baseRule = baseRule
        
    def rule(self, pair):
        return self._baseRule(pair, self._thres1, 
                              self._thres2)

def canopiesTopK(recs, topK, rule1, rule2, rule3):
    rightSet = list(recs)
    shuffle(rightSet)
    
    #build an invertedIndex
    invIndex = {_getId(rec): rec 
                for rec in recs}
    
    #count the pairwise comparisons
    pairwiseCompsStage1 = 0
    
    #start the timer now
    timerStage1 = time.time()
    canopies = list()
    while len(rightSet) > 0:
        currentCanopy = list([rightSet.pop(0)])
        i=0
        while i < len(rightSet):
            rec = rightSet[i]
            pairwiseCompsStage1 += 1
            if rule2((currentCanopy[0], rec)):
                currentCanopy.append(rightSet.pop(i))
                #do not increment i
                continue
            elif rule1((currentCanopy[0], rec)):
                currentCanopy.append(rec)
            i += 1
        #append the new canopy to the list of all canopies
        canopies.append(currentCanopy)
    #end of stage 1
    timerStage1 = time.time() - timerStage1
    
    print 'stage 1 completed: ' + str(timerStage1)
    
    print 'pairwise comparisons stage 1: ' + str(pairwiseCompsStage1)
    sys.stdout.flush()
    
    #start stage 2 
    timerStage2 = time.time()
    parentPointers = {recId:-1
                      for recId in invIndex.keys()}
    
    pairwiseCompsStage2 = 0
    for canopy in canopies:
        for pair in combinations(canopy, 2):
            pairwiseCompsStage2 += 1
            id1 = _getId(pair[0]); id2 = _getId(pair[1])
            root1 = _getRoot(id1, parentPointers) 
            root2 = _getRoot(id2, parentPointers) 
            if root1 != root2:
                #check if they should merge
                if rule3(pair):
                    parentPointers[min(root1, root2)] = max(root1, root2)
    #end of stage 2
    timerStage2 = time.time() - timerStage2
    
    print 'stage 2 completed: ' + str(timerStage2)
    
    print 'pairwise comparisons stage 2: ' + str(pairwiseCompsStage2)
    sys.stdout.flush()
    
    entityDict = defaultdict(list)
    
    #stage3 - build the entities
    timerStage3 = time.time()
    for recId in invIndex.keys():
        entityDict[_getRoot(recId, parentPointers)
                   ].append(recId)
    entities = [[invIndex[recId]
                 for recId in entityDict[entityId]]
                for entityId in entityDict.keys()]
    entities = sorted(entities, key=lambda x: len(x), reverse=True)
    #end of stage 3
    timerStage3 = time.time() - timerStage3
    
    print 'stage 3 completed: ' + str(timerStage3)
    
    #let's also compute an lower bound for stage 1 time,
    #i.e., how low we could get on stage 1 with a heavily optimized implementation
    
    recsLeft = [sample(recs, 1)[0]
                for _ in range(pairwiseCompsStage1)]
    recsRight = [sample(recs, 1)[0]
                 for _ in range(pairwiseCompsStage1)]
    timerLowerBound = time.time()
    for i in range(pairwiseCompsStage1):
        rule2((recsLeft[i], recsRight[i]))
    timerLowerBound = time.time() - timerLowerBound
    
    timer = timerStage1+timerStage2+timerStage3
    
    return (entities[:topK], timer, timerStage1,
            timerStage2, timerStage3, timerLowerBound)
    
def _getId(rec):
    return rec['original']['other']['internalId']

def _getRoot(recId, parentPointers):
    while parentPointers[recId] != -1:
        recId = parentPointers[recId]
    return recId
                
                
def main(argv=None):
    fields = ['title', 'authors']
    K=20
    allRecords, _ = loadDBLPcomplete(filename=argv[1])
    allRecords = sample(allRecords, 1000)
    records = shingling(allRecords, fields, 
                        processingFunction=DBLPprocessV4)
    nR1 = titleAuthorsSimpleANDRule(0.3, 0.3)
    nR2 = titleAuthorsSimpleANDRule(0.6, 0.6)
    nR3 = titleAuthorsSimpleANDRule(0.5, 0.5)
    entities, timer, timerStage1,\
    timerStage2, \
    timerStage3, timerLB = canopiesTopK(records,
                                        K, 
                                        rule1=nR1.rule, 
                                        rule2=nR2.rule, 
                                        rule3=nR3.rule)
    for k in range(1,K+1):
        print 'k: ' + str(k)
        print 'stage 1 time: ' + str(timerStage1)
        print 'stage 2 time: ' + str(timerStage2)
        print 'stage 3 time: ' + str(timerStage3)
        print 'overall time: ' + str(timer)
        print 'Lower Bound time: ' + str(timerLB)
        print 'canopies (p,r,f1): ' + str(f1Fast(allRecords, entities[:k], k))
    
if __name__ == '__main__':
    sys.exit(main(sys.argv)) 
                
            
         
    
