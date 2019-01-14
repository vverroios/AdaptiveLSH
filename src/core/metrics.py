'''
Created on Jul 15, 2016

@author: vasilis verroios
'''

import sys

from collections import defaultdict
from itertools import combinations

def _getGoldClusters(gold):
    #let's create the gold clusters for all records
    goldClusters = defaultdict(list)
    for record in gold:
        goldClusters[record['entity']].\
        append(int(record['other']['internalId']))
        
    return goldClusters

def _getOutClusters(output):
    #let's create the gold clusters for all records
    goldClusters = defaultdict(list)
    for record in output:
        goldClusters[record['original']['entity']].\
        append(int(record['original']['other']['internalId']))
        
    return goldClusters

def _flatten(outputK):
    output = list()
    for entityRecs in outputK:
        output.extend(entityRecs)
        
    return output

def _getAllOutRecords(output):
    listOfIds = list()
    for record in output:
        listOfIds.append(int(record['original']['other']['internalId']))
        
        
    entities = set([record['original']['entity']
                    for record in output])
    print str([(entity,
                sum([1 if record['original']['entity'] == entity
                     else 0
                     for record in output]))
                for entity in entities])
        
    return listOfIds

def _getTargetClusters(target):
    clusters = dict()
    i = 0
    for cluster in target:
        clusters[i] = [int(rec['original']['other']['internalId'])
                       for rec in cluster]
        i += 1
    return clusters

def reduction(dataset, output):
    return len(_flatten(output))/float(len(dataset))

def recoverF1(gold, output, k, target=None, maximize=2):
    '''
    compute the recovered output by reconstructing all entities 
    (assume an oracle comparison function between all records in the 
    reduced set R and the remaining records, i.e., F-R, where F is the 
    full set of records - see Verroios and Molina top-K Entity Resolution) 
    return the precision, recall, and f1 when comparing 
    the recovered output to the 
    set of all records in the gold standard top-K 
    and the percentage of the top-k entities covered
    
    Arguments:
    gold -- a list of all records, where each record is a dict 
    (under ['entity'] we find the entity the record refers to)
    (under ['other']['internalId'] we find the internal id)
    output -- the list of the output records, where each record 
    is a dict (under ['original']['entity'] we find the entity 
    the record refers to) (under ['original']['other']['internalId'] 
    we find the internal id)
    k -- parameter k in the top-K
    target -- in case we compare with the output of another (target) 
    algorithm
    maximize -- this is used in case there are equalities in the 
    gold (target) clusters (e.g., the k-th and the (k+1)-th clusters
    have the same size). 
    maximize: precision for 0, recall for 1, f1 for 2
    
    Output:
    (percentageOfEntitiesCovered, precision, recall, F1)
    '''
    
    if target is None:
        goldClusters = _getGoldClusters(gold)
    else:
        goldClusters = _getTargetClusters(target)
        
    allClusters = sorted(goldClusters.iteritems(), 
                         key = lambda c: len(c[1]), 
                         reverse = True)
    
    #let's get the top-k entities
    #we will get more than k if there are many entities 
    #with the same number of records at the k-th position
    topkClusters = allClusters[:k]
    
    sizeOfKth = len(topkClusters[-1][1])
    choose = len([c for c in topkClusters
                  if len(c[1]) == sizeOfKth])
    allKth = [c for c in allClusters 
              if len(c[1]) == sizeOfKth]
    topkClusters = allClusters[:k-choose+len(allKth)]
    
    #topK entities
    topEntities = set([c[0] for c in topkClusters])
    
    #output entities
    outEntities = set([rec['original']['entity'] 
                       for entity in output
                       for rec in entity])
    
    #we assume that we recover the full entities from the output
    recovered = [cluster 
                 for cluster in allClusters
                 if cluster[0] in outEntities][:k]
                 
    #let's get the set of records referring to the top-K gold entities
    topKrecordList, _ = _getTopKrecords(goldClusters, k)
    #now let's get the set of records in the output
    outputRecordSet = set([recId 
                           for entity in recovered
                           for recId in entity[1]])
    
    best = None
    
    for topKrecords in topKrecordList:
        topKrecordSet= set(topKrecords)
        #now let's compute the precision
        precision = len(outputRecordSet.intersection(topKrecordSet)) / float(len(outputRecordSet))
        recall = len(outputRecordSet.intersection(topKrecordSet)) / float(len(topKrecordSet))
        if (precision+recall) == 0.0:
            current = (0.0, 0.0, 0.0)
        else:
            f1 = 2*precision*recall / (precision + recall)
            current = (precision, recall, f1)
        
        if best is None or best[maximize] < current[maximize]:
            best = current
        
    percentageEntitiesCovered = len(outEntities.intersection(topEntities)) / float(len(topEntities))
    
    return (percentageEntitiesCovered, 
            best[0], best[1], best[2])

def actualPercentage(gold, k):
    goldClusters = _getGoldClusters(gold)
    #let's get the set of records referring to the top-K gold entities
    topKrecords = _getTopKrecords(goldClusters, k)[0][0]
    allRecords = [recId
                  for cluster in goldClusters.iteritems()
                  for recId in cluster[1]]
    
    return float(len(topKrecords))/len(allRecords)

def meanPrecision(gold, output, k, target=None, recover=False):
    if target is None:
        goldClusters = _getGoldClusters(gold)
    else:
        goldClusters = _getTargetClusters(target)
    
    #let's get the list of possible topK clusters
    _, topKclusterList = _getTopKrecords(goldClusters, k)
    
    #now let's find the true clusters in the output
    outputClusters = _getOutClusters(_flatten(output))
    
    cc = compareClusters(goldClusters)
    if recover:
        outputClustersRecords = sorted(outputClusters.iteritems(),
                                       cmp=lambda x,y: cc.cmp2(x, y)
                                       )[:k]
    else:
        outputClustersRecords = sorted(outputClusters.iteritems(),
                                       cmp=lambda x,y: cc.cmp(x, y)
                                       )[:k]
                                   
    bestMP = None
    bestMR = None
    for topKclusters in topKclusterList:
        if recover:
            mp, mr = _meanAveragePrecisionBulk(topKclusters, outputClustersRecords, 
                                               goldClusters=goldClusters)
        else:
            mp, mr = _meanAveragePrecisionBulk(topKclusters, outputClustersRecords)
        if bestMP is None or bestMP < mp:
            bestMP = mp
            bestMR = mr
    
    return bestMP, bestMR

class compareClusters:
    def __init__(self, gold):
        self._gold = gold
        
    def cmp(self, x, y):
        if len(x[1]) > len(y[1]):
            return -1
        elif len(x[1]) < len(y[1]):
            return 1
        else:
            if len(self._gold[x[0]]) > len(self._gold[y[0]]):
                return -1
            else:
                return 1
            
    def cmp2(self, x, y):
        if len(self._gold[x[0]]) > len(self._gold[y[0]]):
            return -1
        else:
            return 1

def _meanAveragePrecisionEntities(topKclusters, outputClustersRecords):
    sizesInTopK = set([len(c[1])
                       for c in topKclusters])
    checkAt = sorted([max([i
                           for i in range(len(topKclusters))
                           if len(topKclusters[i][1]) == s])
                      for s in sizesInTopK])
    mp = list()
    previousCheck = -1
    for check in checkAt:
        goldEntities = set([topKclusters[i][0]
                            for i in range(check+1)])
        topIEntities = set([outputClustersRecords[i][0]
                           for i in range(check+1)])
        precision = len(topIEntities.intersection(goldEntities))/float(len(topIEntities))
        for _ in range(check-previousCheck):
            mp.append(precision)
        previousCheck = check
        
    print 'precisions: ' + str(mp)
    
    return sum(mp)/len(mp)
    

def _meanAveragePrecision(topKclusters, outputClustersRecords, 
                          goldClusters=None):
    sizesInTopK = set([len(c[1])
                       for c in topKclusters])
    checkAt = sorted([max([i
                           for i in range(len(topKclusters))
                           if len(topKclusters[i][1]) == s])
                      for s in sizesInTopK])
    
    previousCheck = -1
    precisions = list()
    recalls = list()
    for check in checkAt:
        goldRecords = set([record
                           for i in range(previousCheck+1, check+1)
                           for record in topKclusters[i][1]])
        if goldClusters is not None:
            topIEntities = set([outputClustersRecords[i][0]
                                for i in range(previousCheck+1, check+1)])
            topIRecords = set([record
                               for entity in topIEntities
                               for record in goldClusters[entity]])
        else:
            topIRecords = set([record
                               for i in range(previousCheck+1, check+1)
                               for record in outputClustersRecords[i][1]])
            
        precision = len(topIRecords.intersection(goldRecords))/float(len(topIRecords))
        recall = len(topIRecords.intersection(goldRecords))/float(len(goldRecords))
        for _ in range(check-previousCheck):
            precisions.append(precision)
            recalls.append(recall)
        previousCheck = check
        
    APs = [sum(precisions[:i])/i
           for i in range(1,len(precisions)+1)]
    ARs = [sum(recalls[:i])/i
           for i in range(1,len(recalls)+1)]
        
    print 'precisions: ' + str(APs)
    print 'recalls: ' + str(ARs)
    
    return sum(APs)/len(APs), sum(ARs)/len(ARs) 

def _meanAveragePrecisionBulk(topKclusters, outputClustersRecords,
                              goldClusters=None):
    sizesInTopK = set([len(c[1])
                       for c in topKclusters])
    checkAt = sorted([max([i
                           for i in range(len(topKclusters))
                           if len(topKclusters[i][1]) == s])
                      for s in sizesInTopK])
    
    mp = list()
    mr = list()
    previousCheck = -1
    for check in checkAt:
        goldRecords = set([record
                           for i in range(check+1)
                           for record in topKclusters[i][1]])
        if goldClusters is not None:
            topIEntities = set([outputClustersRecords[i][0]
                                for i in range(check+1)])
            topIRecords = set([record
                               for entity in topIEntities
                               for record in goldClusters[entity]])
        else:
            topIRecords = set([record
                               for i in range(check+1)
                               for record in outputClustersRecords[i][1]])
        precision = len(topIRecords.intersection(goldRecords))/float(len(topIRecords))
        recall = len(topIRecords.intersection(goldRecords))/float(len(goldRecords))
        for _ in range(check-previousCheck):
            mp.append(precision)
            mr.append(recall)
        previousCheck = check
        
    print 'precisions: ' + str(mp)
    print 'recalls: ' + str(mp)
    
    return sum(mp)/len(mp), sum(mr)/len(mr) 

def f1(gold, output, k, target=None, maximize=2):
    '''
    return the precision, recall, and f1 when comparing an
    the set of all records the algorithm outputs to the 
    set of all records in the gold standard top-K
    
    Arguments:
    gold -- a list of all records, where each record is a dict 
    (under ['entity'] we find the entity the record refers to)
    (under ['other']['internalId'] we find the internal id)
    output -- the list of the output records, where each record 
    is a dict (under ['original']['entity'] we find the entity 
    the record refers to) (under ['original']['other']['internalId'] 
    we find the internal id)
    k -- parameter k in the top-K
    target -- in case we compare with the output of another (target) 
    algorithm
    maximize -- this is used in case there are equalities in the 
    gold (target) clusters (e.g., the k-th and the (k+1)-th clusters
    have the same size). 
    maximize: precision for 0, recall for 1, f1 for 2
    '''
    
    if target is None:
        goldClusters = _getGoldClusters(gold)
    else:
        goldClusters = _getTargetClusters(target)
    
    #let's get the set of records referring to the top-K gold entities
    topKrecordList, _ = _getTopKrecords(goldClusters, k)
    print 'cluster sizes: ' + str(sorted([len(x) for x in output], reverse=True))
    #now let's get the set of records in the output
    outputRecords = _getAllOutRecords(_flatten(output))
    
    
    outputRecordSet = set(outputRecords) 
    
    best = None
    
    for topKrecords in topKrecordList:
        topKrecordSet= set(topKrecords)
        #now let's compute the precision
        precision = len(outputRecordSet.intersection(topKrecordSet)) / float(len(outputRecordSet))
        recall = len(outputRecordSet.intersection(topKrecordSet)) / float(len(topKrecordSet))
        if (precision+recall) == 0.0:
            current = (0.0, 0.0, 0.0)
        else:
            f1 = 2*precision*recall / (precision + recall)
            current = (precision, recall, f1)
        
        if best is None or best[maximize] < current[maximize]:
            best = current  
        
    return best

def f1Fast(gold, output, k, target=None, log=True):
    '''
    "Efficient implementation"
    return the precision, recall, and f1 when comparing an
    the set of all records the algorithm outputs to the 
    set of all records in the gold standard top-K
    
    Arguments:
    gold -- a list of all records, where each record is a dict 
    (under ['entity'] we find the entity the record refers to)
    (under ['other']['internalId'] we find the internal id)
    output -- the list of the output records, where each record 
    is a dict (under ['original']['entity'] we find the entity 
    the record refers to) (under ['original']['other']['internalId'] 
    we find the internal id)
    k -- parameter k in the top-K
    target -- in case we compare with the output of another (target) 
    algorithm
    maximize -- this is used in case there are equalities in the 
    gold (target) clusters (e.g., the k-th and the (k+1)-th clusters
    have the same size). 
    '''
    
    if log:
        print 'starting F1 fast'
        sys.stdout.flush()
    
    if target is None:
        goldClusters = _getGoldClusters(gold)
    else:
        goldClusters = _getTargetClusters(target)
    
    #let's get the set of records referring to the top-K gold entities
    print 'cluster sizes: ' + str(sorted([len(x) for x in output], reverse=True))
    #now let's get the set of records in the output
    outputRecords = _getAllOutRecords(_flatten(output))
    
    
    outputRecordSet = set(outputRecords) 
    
    topKrecords, _ = _getTopKrecordsFast(goldClusters, k, outputRecordSet)
    
    topKrecordSet= set(topKrecords)
    precision = len(outputRecordSet.intersection(topKrecordSet)) / float(len(outputRecordSet))
    recall = len(outputRecordSet.intersection(topKrecordSet)) / float(len(topKrecordSet))
    if (precision+recall) == 0.0:
        return (0.0, 0.0, 0.0)
    else:
        f1 = 2*precision*recall / (precision + recall)
        
    return (precision, recall, f1)    

def _getTopKrecords(goldClusters, k):
    allClusters = sorted(goldClusters.iteritems(), 
                         key = lambda c: len(c[1]), 
                         reverse = True)
    
    topkClusters = allClusters[:k]

    sizeOfKth = len(topkClusters[-1][1])
    
    choose = len([c for c in topkClusters
                  if len(c[1]) == sizeOfKth])
    
    allKth = [c for c in allClusters 
              if len(c[1]) == sizeOfKth]
    
    fixed = topkClusters[:(k-choose)]

    topKrecordsList = list()
    topKclustersList = list()
    
    for chosen in combinations(allKth, choose):
        topKrecords = list()
        for cluster in fixed+list(chosen):
            topKrecords.extend(cluster[1])
        topKrecordsList.append(topKrecords)
        topKclustersList.append(fixed+list(chosen))
        
    return topKrecordsList, topKclustersList

def _getTopKrecordsFast(goldClusters, k, outputRecSet, log=True):
    
    if log:
        print 'starting getTopK fast'
        sys.stdout.flush()
    
    allClusters = sorted(goldClusters.iteritems(), 
                         key = lambda c: len(c[1]), 
                         reverse = True)
    
    if log:
        print 'get topK gold-fast: sorted all clusters'
        sys.stdout.flush()
    
    topkClusters = allClusters[:k]

    sizeOfKth = len(topkClusters[-1][1])

    if log:
        print 'get topK gold-fast: size of Kth - ' + str(sizeOfKth)
        sys.stdout.flush()
    
    choose = len([c for c in topkClusters
                  if len(c[1]) == sizeOfKth])
    
    if log:
        print 'get topK gold-fast: choose - ' + str(choose)
        sys.stdout.flush()
    
    allKth = [c for c in allClusters 
              if len(c[1]) == sizeOfKth]
    
    if log:
        print 'get topK gold-fast: found allKth'
        sys.stdout.flush()

    allKth = sorted(allKth, 
                    key=lambda x: len(outputRecSet.intersection(x[1])), 
                    reverse=True)
    
    if log:
        print 'get topK gold-fast: sorted allKth based on intersection with output'
        sys.stdout.flush()
    
    fixed = topkClusters[:(k-choose)]

    topKrecords = [rec 
                   for cluster in (fixed+allKth[:choose])
                   for rec in cluster[1]]
    topKclusters = fixed+allKth[:choose]
    
    if log:
        print 'get topK gold-fast: ready to return'
        sys.stdout.flush()
    
    return topKrecords, topKclusters
            
    
    
