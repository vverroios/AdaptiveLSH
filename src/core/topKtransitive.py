'''
Created on Jul 12, 2016

@author: vasilis verroios
'''

import time
import sys
import pickle
import os

from collections import defaultdict
from itertools import combinations
from math import log

from hashFunctions import minHash, shingling, jaccardSim, \
titleAuthorsSimANDothers, advancedProcessField, spotsigs, \
titleAuthorsSimANDvenue, titleANDAuthorsSim, titleANDAuthorsSim2, \
DBLPprocess, DBLPprocessV1, titleAuthorsSimANDothers2, \
DBLPprocessV2, sketches
from lsh.utilities import loadCora, loadSpotSigs, \
loadPreProSpotSigs, preprocessSpotSigs, loadDBLP, \
computeANDHashSchemeGivenBudget, computeHashSchemeGivenBudget, loadDBLPcomplete
from cost_model.model import pairCost
from metrics import f1
from random import randint, sample
from pip._vendor.html5lib.constants import entities
from numpy.ma.core import mean
from core.hashFunctions import advancedProcessFieldDBLP, CORAOTHERS,\
    CORATITLEAUTHORS, DBLPTITLE, DBLPAUTHORS, DBLPprocessV3, DBLPprocessV4
from core.traditionalBlocking import traditionalBlocking
from core.metrics import f1Fast

MODES = ['adaLSH', 'adaLSHnoP', 'Pairs', 'LSH', 'LSHnoP']

class allPairs:
    '''
    goes over all pairs of records and 
    applies transitive closure on the pairs 
    with a similarity greater than the thres
    '''
    
    def __init__(self, records, light=True, skipMerged=False):
        """
        initialization: just stores the records
        
        Arguments:
        records -- a list with all the records.
        light -- True: just compute the pairwise similarities
        without performing the merging (for performance evaluation)
        skipMerged -- True: before checking if there is an edge 
        between two records, check if the records are already merged
        """
        self.records = records
        self._light = light
        self._skipMerged = skipMerged
    
    def topK(self, k, thres, fields, erRule=None, focus=[]):
        """
        finds the top-K(in terms of number of records-per-entity) entities 
        based on a similarity threshold; this threshold is used for transitive closure
        
        Arguments:
        k -- the number of top entities to return
        fields -- the fields used in the similarity, 
        we take the average similarity over all those fields
        thres -- the similarity threshold on which records 
        are assumed to be of the same entity
        """
        #initialize clusters
        self._clusters = {int(record['original']['other']['internalId']):
                          set([int(record['original']['other']['internalId'])])
                          for record in self.records}
        
        self._recordIdMap = {int(record['original']['other']['internalId']):
                             record
                             for record in self.records}
        
        #go over all pairs
        self.timer = time.time()
        
        pairsCount = (len(self.records)*(len(self.records)-1))/2
        count = 0
        print '|',
        for pair in combinations(self.records, 2):
            count += 1
            if count % (pairsCount/15) == 0:
                print '-',
                sys.stdout.flush()
                
            if self._skipMerged:
                ids = [int(record['original']['other']['internalId'])
                       for record in pair]
                if ids[1] in self._clusters[ids[0]]:
                    continue
            
            if erRule is None:
                if jaccardSim(pair, fields) >= thres:
                    self._merge(pair)
                else:
                    for ent in focus:
                        if (pair[0]['original']['entity'] == ent 
                            and pair[1]['original']['entity'] == ent):
                            print 'did not merge: ' + str(pair)
                            break
            else:
                mustMatch = False
                for ent in focus:
                        if (pair[0]['original']['entity'] == ent 
                            and pair[1]['original']['entity'] == ent):
                            mustMatch = True
                            break
                if erRule(pair, mustMatch=mustMatch):
                    self._merge(pair)
                    if ((pair[0]['original']['entity'] in focus or
                        pair[1]['original']['entity'] in focus) 
                        and
                        not mustMatch): 
                        print 'should not merge: ' + str(pair)
                    
        print '|'
        
        if not self._light:
            self._buildEntities(k)
        
        self.timer = time.time() - self.timer
        
        
        if self._light:
            return 
        
        return self._entities
        
    def _buildEntities(self, k):
        #find the top-k clusters
        self._entities = buildEntities(k, 
                                       self._clusters, 
                                       self._recordIdMap)
                
    def _merge(self, pair):
        if self._light: 
            return
        
        ids = [int(record['original']['other']['internalId'])
               for record in pair]
        #note that ids[0] is always less than ids[1] 
        unionSet = self._clusters[ids[0]].union(self._clusters[ids[1]])
        for recId in unionSet:
            self._clusters[recId] = unionSet

def buildEntities(k, clusters, idMap):
    topClusters = list()
    while len(topClusters) < k:
        topCluster = max(clusters.iteritems(), 
                         key = lambda c: len(c[1]))
        for recId in topCluster[1]:
            clusters[recId] = list()
        topClusters.append(topCluster)
        
    return [[idMap[recId]
             for recId in cluster[1]]
            for cluster in topClusters]
    

class adaLSH:
    """
    runs LSH in an adaptive way so that more hash evaluations are 
    applied on strong (for being a topK entity) candidates and less 
    hash evaluations are applied on weak candidates 
    """
    
    EPSUPPERS = 0.005
    EPSSGAP = 0.3
    EPSILON = 0.1
    PAIRSRECSAMPLE = 50#4#100
    
    def __init__(self, records, numberOfSchemes=None, 
                 topBuckets=False, mode='adaLSH', 
                 schemeMode='expo', transitivity=True):
        """
        initialization: just stores the records and all necessary information 
        to apply the algorithm
        
        Arguments:
        ---------------------
        records -- a list with all the records
        ---------------------
        numberOfSchemes -- (Not used in the current implementation) 
        None if we want to run the adaptive LSH.
        ---------------------
        topBuckets -- (Not used in the current implementation)
        True: run the topBuckets algorithm
        1 if we want to run the classic(single-Scheme) LSH.
        ---------------------
        mode -- one of ['adaLSH', 'adaLSHnoP', 'Pairs', 'LSH', 'LSHnoP']
        ---------------------
        schemeMode -- 'expo' OR 'sX' for sequential with a step of X
        ---------------------
        transitivity -- False in case we just want to run LSH without transitivity
        """
        
        if mode not in MODES:
            raise ValueError('Unknown mode: ' + str(mode))
        
        #pairsOn -- True: run the pairs function based on a cost model
        if 'noP' in mode: 
            self._pairsOn = False
        else:
            self._pairsOn = True
        
        self._trans = transitivity
        self._mode = mode
        if 'CN' in schemeMode:
            #e.g., for a costNoise of 0.5, the pairs cost computed is half the real one
            self._costNoise = float(schemeMode[-3:])
            self._schemesMode = schemeMode[:-5]
        else:
            self._costNoise = 1.0
            self._schemesMode = schemeMode
        self.records = records
        self._numberOfSchemes = numberOfSchemes
        self._topBuckets = topBuckets
        self._ids = 0
        self._globalBuckets = [dict() for _ in range(int(log(len(records),2))+1)]
        
        #we cache optimization results, in order not to run the program again and again
        self._optFilename = './optimization/optResults.pickle'
        if os.path.isfile(self._optFilename): 
            with open(self._optFilename, 'rb') as handle:
                self._optRes = pickle.load(handle)
        else:
            self._optRes = dict()
            with open(self._optFilename, 'wb') as handle:
                pickle.dump(self._optRes, handle)
        
    def _linearCost(self, trials=10000):#main exps run with trials=100
        """
        returns the avg time needed to apply a single hash function
        
        Arguments:
        trials -- the number of trials to compute the average
        """
        tScheme = self._schemes[-1]
        if len(tScheme) == 3:
            x,y = tScheme[:2]
        else:
            x,y,z = tScheme[:3]; x+=z
            
        #select uniformly at random the functions and records
        probes = [(randint(0,x-1), 
                   randint(0,y-1), 
                   randint(0,len(self.records)-1))
                  for _ in range(trials)]
        
        t = time.time()
        for i in range(trials):
            self._f.h(probes[i][0], 
                      probes[i][1], 
                      self.records[probes[i][2]])
        
        self._linear = (time.time() - t) / trials
    
    def topK(self, k, thres, fields=None, erRule=None, 
             tracesOn=False, cosine=False):
        """
        finds the top-K(in terms of number of records-per-entity) entities 
        based on a similarity threshold; this threshold is used for transitive closure
        
        Arguments:
        k -- the number of top entities to return
        thres -- the jaccard similarity to decide if two records refer
        fields -- the fields where similarity is computed on
        erRule -- the ER rule used to decide if two records refer 
        to the same entity
        tracesOn -- True/False: indicates if we should keep track metrics 
        regarding the implementation
        cosine -- False/True: indicates if cosine distance is used 
        (default - False)
        """
        
        self._fields = self.records[0]['data'].keys() if erRule is None else None
        self._erRule = erRule
        self._thres = thres[0] if len(thres) == 1 else None
        self._tracesOn = tracesOn
        self._cos = cosine
        
        
        #first find the schemes for this similarity threshold
        self._buildSchemes(thres)
        tScheme = self._schemes[-1]
        if not self._cos:
            self._f = minHash(self.records, 
                              tScheme[0], 
                              tScheme[1],
                              None if len(tScheme) == 3 
                              else tScheme[2],
                              fields)
        else:
            self._f = sketches(self.records, 
                               tScheme[0], 
                               tScheme[1])
        
        #create the all tables
        y = tScheme[1]
        self._splitPoint = tScheme[0]   
        
        print 'top-scheme: ' + str(self._schemes[-1])
        print 'all schemes: ' + str(self._schemes)
        
        self._hashtables = [defaultdict(list) for _ in range(y)]
        
        #let's compute the cost model parameters
        self._linearCost()
        self._pairCost = pairCost(sample(self.records,
                                         adaLSH.PAIRSRECSAMPLE), 
                                  self._fields, erRule)
        print 'linear cost: ' + str(self._linear)
        print 'pairwise cost: ' + str(self._pairCost)
        
        #if mode is Pairs then build a tree to run the P function on
        #Note that the time to build this tree is not taken into account
        if self._mode == 'Pairs':
            recordsTree = self._buildGlobalTree()
        
        self.timer = time.time()
        self.hashEvaluations = 0
        
        #start with the minimum scheme for all records
        self._phase = 0
        if self._tracesOn:
            self._stepsToRoot = list()
            self.stepsToRootAvg = list()
            self.stepsToRootMin = list()
            self.stepsToRootMax = list()
                
        if self._topBuckets:
            self._schemes = [self._schemes[-1]]
            self._applyScheme()
            self.timer = time.time() - self.timer
            return self._topBucketEntities(k)
        else:
            print '---------------------'
            print 'Round ' + str(self._phase)
            t1 = time.time()
            if self._mode == 'Pairs':
                self._applyPairs(recordsTree)
            else:    
                self._applyScheme()
            print 'processed in: ' + str(time.time() - t1)
            
        if not self._trans:
            #if we just want to run LSH
            return self._hashtables
        
        entities = list()
        #while less than k entities have been approved 
        while len(entities) < k:
            if self._tracesOn:
                self.stepsToRootAvg.append(mean(self._stepsToRoot))
                self.stepsToRootMin.append(min(self._stepsToRoot))
                self.stepsToRootMax.append(max(self._stepsToRoot))
                self._stepsToRoot = list()
                
            self._phase += 1
            #get the top global bucket
            topRecords = self._getTopBucket()
            print '---------------------'
            print 'Round ' + str(self._phase)
            print 'top-size: ' + str(topRecords['size'])
            print 'level: ' + str(topRecords['first']['level'])

            self._delTop(topRecords)
            #if it is under the top scheme 
            #(i.e., the likelihood of FP or FN is very low)
            #and since there is NO global bucket (entity)
            #that is larger (even the ones under lower level 
            #schemes), we can assume this entity is a topK
            #entity
            if (topRecords['first']['level'] == self._topScheme 
                or self._mode == 'LSHnoP'):
                entities.append(self._children(topRecords))
                print 'added in final entities'
            else:
                t1 = time.time()
                #let's apply a higher level scheme on the topBucket
                self._applyScheme(topRecords=topRecords)
                print 'processed in: ' + str(time.time() - t1)
            
                
        self.timer = time.time() - self.timer
        
        #let's serialize the new (if any) optimization results 
        with open(self._optFilename, 'wb') as handle:
            pickle.dump(self._optRes, handle)
        
        return entities
    
    def _buildGlobalTree(self):
        previousRecord = self.records[0]
        previousRecord['level'] = 0
        self._ids += 1
        root = {'first': previousRecord, 
                'last': self.records[-1], 
                'size': len(self.records), 
                'id': self._ids, 
                'parent': None}
        for record in self.records[1:]:
            previousRecord['right'] = record
            record['level'] = 0
            previousRecord = record
        record['right'] = None
            
        return root
    
    def _buildSchemesBasedOnSshape(self, thres):
        self._schemes = list()
        x=2;y=1
        #while the S shape is not sharp enough
        while True:
            while self._pr(x,y,thres) < 1.0-adaLSH.EPSUPPERS:
                y += 1
            self._schemes.append((x,y))
            if self._pr(x,y,thres-adaLSH.EPSSGAP) < adaLSH.EPSUPPERS:
                break
            x += 1 
            
        if self._numberOfSchemes is not None:
            if self._numberOfSchemes == 1:
                self._schemes = [self._schemes[-1]]
            else:
                schemes = [self._schemes[i]
                           for i in range(0,
                                          len(self._schemes),
                                          len(self._schemes)/self._numberOfSchemes)]
                schemes[-1] = self._schemes[-1]
                self._schemes = schemes
        
        if self._mode == 'LSH':    
            #take the P function into account as well
            self._topScheme = len(self._schemes)
        else:
            self._topScheme = len(self._schemes)-1  
          
        
    def _buildSchemes(self, thresholds):
        if 'expo' in self._schemesMode:
            if self._schemesMode == 'expo':
                til = 10
            else:
                til = int(self._schemesMode[4:])
            budgets = [20*(2**i) for i in range(til)]
        elif self._schemesMode[0] == 's':
            maxBudget = 3840#5120
            perStep = int(self._schemesMode[1:])
            if (('LSH' in self._mode and 'ada' not in self._mode)
                or 'Pairs' in self._mode):
                howmany = 1
            else:
                howmany = maxBudget/perStep
            budgets = [perStep*i for i in range(1, howmany + 1)]
            
        if 'LSH' in self._mode and 'ada' not in self._mode:
            budgets = [budgets[0]]
            
        self._schemes = list()
        epsilon = adaLSH.EPSILON
        exact = True
        
        if len(thresholds) == 2:
            #AND between two groups of fields
            res = (0,0)
            for B in budgets:
                optCase = (thresholds[0], B, epsilon, exact, self._cos, 
                           thresholds[1], res[0], res[1])
                if optCase in self._optRes:
                    res = self._optRes[optCase]
                else:
                    res = computeANDHashSchemeGivenBudget(B, epsilon, thresholds[0], 
                                                          thresholds[1], exact=exact, 
                                                          wConstraint=res[0],
                                                          uConstraint=res[1])
                    self._optRes[optCase] = res
                self._schemes.append([int(res[0]), int(res[2]), 
                                      int(res[1]), int(res[3]), 
                                      int(res[4])])
        elif len(thresholds) == 1:
            for B in budgets:
                optCase = (thresholds[0], B, epsilon, exact, self._cos)
                if optCase in self._optRes:
                    res = self._optRes[optCase]
                else:
                    res = computeHashSchemeGivenBudget(B, epsilon, 
                                                       thresholds[0], 
                                                       exact=exact,
                                                       cosine=self._cos)
                    self._optRes[optCase] = res
                self._schemes.append(map(int,res))
        else:
            raise ValueError('cannot support more than two thresholds \
                              in the current implementation')
            
        if self._mode == 'LSH':    
            #take the P function into account as well
            self._topScheme = len(self._schemes)
        else:
            self._topScheme = len(self._schemes)-1
    
    def _pr(self, x, y, s):
        return 1.0 - (1.0-s**x)**y 
    
    def _getTopBucket(self):
        for i in range(int(log(len(self.records),2)), -1, -1):
            if len(self._globalBuckets[i]) > 0:
                return max(self._globalBuckets[i].iteritems(), 
                           key=lambda b: b[1]['size'])[1]
    
    def _children(self, topRecords):
        recs = list()
        record = topRecords['first']
        while record is not None:
            recs.append(record)
            record = record['right']
        return recs
    
    def _applyScheme(self, topRecords=None):
        if topRecords is None:
            #apply the level 0 scheme to all records
            for record in self.records:
                record['level'] = 0
                record['parent'] = None
                record['right'] = None
                self._applyHashing(record, None, 
                                   self._schemes[0])
                for i in range(self._schemes[0][1]):
                    self._addtoTable(record,i)
        else:
            record = topRecords['first']
            nextLevel = record['level'] + 1
            if self._pairsOn:
                if ((self._mode == 'LSH' and self._phase > 0)
                    or self._pairsCostLower(topRecords['size'], 
                                            self._schemes[nextLevel], 
                                            self._schemes[nextLevel-1])):
                    self._applyPairs(topRecords)
                    return
                
            while record is not None:
                #reset the global node tree for this node
                record['level'] = nextLevel
                record['parent'] = None
                nextRecord = record['right'] 
                record['right'] = None 
                self._applyHashing(record, 
                                   self._schemes[nextLevel-1], 
                                   self._schemes[nextLevel],)
                for i in range(self._schemes[nextLevel][1]):
                    self._addtoTable(record,i)
                record = nextRecord
                    
    def _applyPairs(self, recordsRoot):
        record = recordsRoot['first']
        records = list()
        print 'applying P at level:' + str(record['level'])
        while record is not None:
            records.append(record)
            #reset the global node tree for this node
            record['level'] = self._topScheme
            record['parent'] = None
            nextRecord = record['right'] 
            record['right'] = None 
            record = nextRecord
        
        #for all pairs
        print '|',
        allP = (len(records)*(len(records)-1))/2
        j = 0
        for pair in combinations(records,2):
            j+=1
            if j%(allP/15)==0:
                print '-',
            #check if they are under the same tree
            root1 = self._root(pair[0]); root2 = self._root(pair[1]) 
            if ('id' in root1 and 
                'id' in root2 and 
                root1['id'] == root2['id']):
                continue
            
            #if not check if there is an edge
            match = False
            if self._erRule is None:
                if jaccardSim(pair, self._fields) > self._thres:
                    match = True
            else:
                if self._erRule(pair):
                    match = True
                    
            if not match:
                continue
            
            for i in [0,1]:
                record = pair[i]
                if record['parent'] is None:
                    self._ids += 1
                    node = {'first': record, 'last': record, 
                            'size': 1, 'id': self._ids, 
                            'parent': None}
                    record['parent'] = node
                    record['right'] = None
                    #it is also a top node
                    self._addTop(node)
                    if i==0:
                        root1 = node
                    else:
                        root2 = node
            
            #merge them 
            self._ids += 1
            newRoot = {'first': root1['first'], 
                       'last': root2['last'], 
                       'size': root1['size']+root2['size'], 
                       'id': self._ids, 
                       'parent': None} 
            root1['parent'] = newRoot; root2['parent'] = newRoot
            root1['last']['right'] = root2['first']
            #update the top global buckets
            self._delTop(root1); self._delTop(root2)
            self._addTop(newRoot)
        print '|'
            
    def _pairsCostLower(self, recs, 
                        schemeNext, 
                        schemePrevious):
        hashFunctionCostNext = (schemeNext[0] + 
                                (0 if len(schemeNext) == 3
                                 else schemeNext[2]))*schemeNext[1] 
        hashFunctionCostPrevious = (schemePrevious[0] + 
                                    (0 if len(schemePrevious) == 3
                                     else schemePrevious[2]))*schemePrevious[1]
        
        if (self._costNoise*(self._pairCost*(recs*(recs-1))/2) 
            < 
            (self._linear*
             (hashFunctionCostNext-hashFunctionCostPrevious)*
             recs)):
            return True
        else:
            return False
                
        
    def _topBucketEntities(self, k):
        
        clusters = [set([int(record['original']['other']['internalId'])])
                    for record in self.records]
        
        recordIdMap = {int(record['original']['other']['internalId']):
                       record
                       for record in self.records}
        
        while len([c for c in clusters if len(c) > 1]) < k:
            #get the largest bucket
            bucket = max([(b,i) 
                          for i in range(len(self._hashtables))
                          for b in self._hashtables[i].iteritems()],
                         key=lambda x: len(x[0][1]))
            #remove this bucket
            self._hashtables[bucket[1]][bucket[0][0]] = list()
            
            #record ids
            ids = set([int(record['original']['other']['internalId'])
                       for record in bucket[0][1]])
            
            #reconstruct clusters
            newClusters = list()
            merged = set()
            for cluster in clusters:
                if len(cluster.intersection(ids)) == 0:
                    newClusters.append(cluster)
                else:
                    merged.update(cluster)
            newClusters.append(merged)
            clusters = newClusters
            
        clusters = [c for c in clusters if len(c) > 1]
        return [[recordIdMap[recId] 
                 for recId in cluster]
                for cluster in clusters]
    
    def _addtoTable(self, record, i):
        #find the bucket
        bucket = self._hashtables[i][record['hashes'][i]]
        bucket.append(record)
        
        #if we just want to run LSH
        if not self._trans:
            return
        
        #if this is the first record to be added in this bucket
        #and this record does not have a parent yet
        if len(bucket) == 1 and record['parent'] is None:
            self._ids += 1
            node = {'first': record, 'last': record, 
                    'size': 1, 'id': self._ids, 
                    'parent': None}
            record['parent'] = node
            record['right'] = None
            #it is also a top node
            self._addTop(node)
        #in case this is the first record to be added in this bucket
        #but this record already has a parent, we don't have to do something
        
        #in case there are more records already added in this bucket
        if len(bucket) > 1:
            #the last node in the bucket has a better chance of 
            #having a straight pointer to the root
            currentRoot = self._root(bucket[-2])
            #if this record does not have a parent yet
            if record['parent'] is None:
                #give this record a parent
                record['parent'] = currentRoot
                previousLast = currentRoot['last'] 
                currentRoot['last'] = record
                previousLast['right'] = record
                #make sure that the node is in the right heap
                self._delTop(currentRoot)
                currentRoot['size'] += 1
                self._addTop(currentRoot)
            else:    
                #if the roots are different
                recordRoot = self._root(record) 
                if currentRoot['id'] != recordRoot['id']:
                    #merge them 
                    self._ids += 1
                    newRoot = {'first': currentRoot['first'], 
                               'last': recordRoot['last'], 
                               'size': currentRoot['size']+recordRoot['size'], 
                               'id': self._ids, 
                               'parent': None} 
                    currentRoot['parent'] = newRoot; recordRoot['parent'] = newRoot
                    currentRoot['last']['right'] = recordRoot['first']
                    #update the top global buckets
                    self._delTop(currentRoot); self._delTop(recordRoot)
                    self._addTop(newRoot)
                    
    def _root(self, node):
        i = 0
        while node['parent'] is not None:
            node = node['parent']
            i += 1
            
        if self._tracesOn:
            self._stepsToRoot.append(i)
                
        return node
    
    def _delTop(self, node):
        key = self._key(node)
        del self._globalBuckets[key][node['id']]
        
    def _addTop(self, node):
        key = self._key(node)
        self._globalBuckets[key][node['id']] = node
        
    def _key(self, node):
        return int(log(node['size'],2))
    
    def _applyHashing(self, r, previous, current):
        #first phase
        if previous is None:
            self.hashEvaluations += (current[0]*current[1] 
                                     if len(current) == 3
                                     else (current[0]+current[2])*current[1])
            
            r['hashes'] = [tuple([self._f.h(self._findIndex(i,current), 
                                            j, r)
                                  for i in range(current[0]+
                                                 (0 if len(current) == 3
                                                  else current[2]))])
                           +tuple(['#|*phase*|#N'+str(self._phase)])
                           for j in range(current[1])]
        #subsequent phases
        else:
            if len(current) == 3:
                self.hashEvaluations += (current[0]*current[1] - previous[0]*previous[1])
                #extend the hash for the pre-existing hashes
                r['hashes'] = [r['hashes'][j][:-1]
                               +tuple([self._f.h(i+previous[0], j, r)
                                       for i in range(current[0]-previous[0])])
                               +tuple(['#|*phase*|#N'+str(self._phase)])
                               for j in range(len(r['hashes']))]
                #add more hashes for the new tables
                r['hashes'].extend([tuple([self._f.h(i, j+previous[1], r)
                                           for i in range(current[0])])
                                    +tuple(['#|*phase*|#N'+str(self._phase)])
                                    for j in range(current[1]-previous[1])])
            else:
                self.hashEvaluations += ((current[0]+current[2])*current[1] 
                                         - (previous[0]+previous[2])*previous[1])
                #extend the hash for the pre-existing hash tables
                r['hashes'] = [r['hashes'][j][:-1]
                               +tuple([self._f.h(self._findIndex(i,current,previous), 
                                                 j, r)
                                       for i in range(current[0]+current[2]-
                                                      previous[0]-previous[2])])
                               +tuple(['#|*phase*|#N'+str(self._phase)])
                               for j in range(len(r['hashes']))]
                #add more hashes for the new tables
                r['hashes'].extend([tuple([self._f.h(self._findIndex(i,current), 
                                                     j+previous[1], r)
                                           for i in range(current[0]+current[2])])
                                    +tuple(['#|*phase*|#N'+str(self._phase)])
                                    for j in range(current[1]-previous[1])])

    def _findIndex(self, i, current, previous=None):
        if len(current) == 3:
            if previous is None:
                return i
            else:
                return i+previous[0]
        else:
            if previous is None:
                return i-current[0]+self._splitPoint if i >= current[0] else i
            else:
                return (i-(current[0]-previous[0])+previous[2]+self._splitPoint 
                        if i>= (current[0]-previous[0])
                        else i+previous[0])

def testDBLP():
    allRecords, _, focus = loadDBLP() 
    
    print 'recs test: ' + str(len(allRecords)) 
    fields = ['title', 'authors']#, 'venue']
    K = 20
    
    records = shingling(allRecords, fields)
    for thres in [0.5]:#[0.4]:#, 0.3]:
        algo2 = allPairs(records, light=False)
        entities = algo2.topK(4*K, thres, fields, focus=focus, erRule=titleANDAuthorsSim)
        entities = sorted(entities, key=lambda x: len(x), reverse=True)
        for k in range(1,K+1):#, 2, 4, 8, 16]:
            print 'k - thres: ' + str(k) + ' - ' + str(thres)
            print 'overall allPairs time: ' + str(algo2.timer)
            print 'allPairs (p,r,f1): ' + str(f1(allRecords, entities[:4*k], k))
    
    return

def testDBLP4(tryCoraRules=True, 
              tryAll=False, 
              tryAND=False):
    allRecords, focus = loadDBLPcomplete(filename='./datasets/DBLP/complete/DBLP-72KP3.csv')
    print 'top entities: ' + str(focus)
    
    print 'recs test: ' + str(len(allRecords)) 
    K = 20
    
    fields = (['title', 'authors', 'other']
              if tryCoraRules 
              else (['title', 'authors', 'venue', 'year']
                    if tryAll
                    else ['title', 'authors']))#, 'venue'])
    if tryCoraRules or tryAll:
        records = shingling(allRecords, fields, 
                            processingFunction=(advancedProcessFieldDBLP 
                                                if tryCoraRules
                                                else DBLPprocess))
    else:
        records = shingling(allRecords, fields, 
                            processingFunction=DBLPprocessV4)
        
    if tryCoraRules:
        for rec in records:
            rec['data']['author'] = rec['data']['authors']
            del(rec['data']['authors'])
    elif tryAll or not tryAND:
        for rec in records:
            rec['data'] = {'all': 
                           set([shingle 
                                for field in fields
                                for shingle in rec['data'][field]])}
    for thres in range(1):
        for mode in [#('LSH', 's1280'), 
                     ('adaLSH', 'expo7')]:
    #     for thres in [0.5]:#[0.4]:#, 0.3]:
            algoADA = adaLSH(records, mode=mode[0], schemeMode=mode[1])
            if tryCoraRules:
                entitiesADA = algoADA.topK(4*K, [CORATITLEAUTHORS, CORAOTHERS], 
                                           [['title', 'author'], ['other']], 
                                           erRule=titleAuthorsSimANDothers)
            elif tryAll:
                entitiesADA = algoADA.topK(4*K, [thres], None)
            else:
                if tryAND:
                    entitiesADA = algoADA.topK(4*K, [DBLPTITLE, DBLPAUTHORS], 
                                               [['title'], ['authors']], 
                                               erRule=titleANDAuthorsSim)
                else:
                    entitiesADA = algoADA.topK(4*K, [thres], None)
            entitiesADA = sorted(entitiesADA, key=lambda x: len(x), reverse=True)
        
            for k in range(1,K+1):#, 2, 4, 8, 16]:
                print 'k - thres: ' + str(k) + ' - ' + str(thres)
                print 'overall adaLSH time: ' + str(algoADA.timer)
                print 'adaLSH (p,r,f1): ' + str(f1Fast(allRecords, entitiesADA[:k], k))
    
    return

def testDBLP2():
    allRecords, _, focus = loadDBLP()
    
    print 'recs test: ' + str(len(allRecords)) 
    fields = ['title', 'authors', 'venue', 'year']
    K = 20
    
    records = shingling(allRecords, fields, processingFunction=DBLPprocess)
    for rec in records:
        rec['data'] = {'all': 
                       set([shingle 
                            for field in fields
                            for shingle in rec['data'][field]])}
        
    for thres in [0.5]:#[0.4]:#, 0.3]:
        algo2 = allPairs(records, light=False)
        entities = algo2.topK(4*K, thres, ['all'])#, focus=focus, erRule=titleANDAuthorsSim)
        entities = sorted(entities, key=lambda x: len(x), reverse=True)
        for k in range(1,K+1):#, 2, 4, 8, 16]:
            print 'k - thres: ' + str(k) + ' - ' + str(thres)
            print 'overall allPairs time: ' + str(algo2.timer)
            print 'allPairs (p,r,f1): ' + str(f1(allRecords, entities[:k], k))
    
    return

def testDBLP3():
    allRecords, _, focus = loadDBLP()
    _glimpseTopEntities(allRecords, focus)
    
    print 'recs test: ' + str(len(allRecords)) 
    fields = ['title', 'authors', 'year']
    K = 20
    
    records = shingling(allRecords, fields, processingFunction=DBLPprocessV2)
        
    for thres in [None]:#[0.4]:#, 0.3]:
        algo2 = allPairs(records, light=False)
        entities = algo2.topK(4*K, thres, fields, 
                              erRule=titleAuthorsSimANDothers2, focus=focus[:1])
        entities = sorted(entities, key=lambda x: len(x), reverse=True)
        for k in range(1,K+1):#, 2, 4, 8, 16]:
            print 'k - thres: ' + str(k) + ' - ' + str(thres)
            print 'overall allPairs time: ' + str(algo2.timer)
            print 'allPairs (p,r,f1): ' + str(f1(allRecords, entities[:4*k], k))
    
    return

def _glimpseTopEntities(allRecords, focus):
    i=1
    for entity in focus:
        j=0
        for rec in allRecords:
            if rec['entity'] == entity:
                j += 1
                print ('top' + str(i) + 
                       '(' + str(j) + '): ' + 
                       str(rec))
                print '*** processed: ' + str({field: DBLPprocess(field, rec[field], rec)
                                               for field in ['title', 'authors', 
                                                             'venue', 'year']})
        i += 1

def testSpotSigs():
    (allRecords, records) = loadPreProSpotSigs()
    K=50
    
    for thres in [0.1,0.2,0.3,0.4]:
        algo2 = allPairs(records, light=False)
        entities = algo2.topK(K, thres, ['data'], focus=[])
        entities = sorted(entities, key=lambda x: len(x), reverse=True)
        for k in [5]:#range(1,K+1):#, 2, 4, 8, 16]:
            print 'k - thres: ' + str(k) + ' - ' + str(thres)
            print 'overall allPairs time: ' + str(algo2.timer)
            for j in range(15):
                print ('allPairs (p,r,f1): ' + 
                       str(f1(allRecords, entities[:(k+j)], k)))
            
    return

def checkEmpty(records):
    for rec in records:
        if len(rec['data']['data']) == 0:
            print rec

def testCora():
    allRecords = loadCora(venueConcat=False)#[:100]
    fields = ['title', 'author', 'other']
    thres = 0.7
    K=20
     
    records = shingling(allRecords, fields, processingFunction=advancedProcessField)
     
    algo2 = allPairs(records, light=False, skipMerged=True)
    targetEntities = algo2.topK(2*K, thres, fields, erRule=titleAuthorsSimANDothers, focus=[])
    print 'overall time: ' + str(algo2.timer)
    algo2 = adaLSH(records,
                   mode='Pairs', 
                   schemeMode='s40')
    entities = algo2.topK(K, [0.7, 0.2], 
                        [['title', 'author'], ['other']], 
                        erRule=titleAuthorsSimANDothers)
    entities = sorted(entities, key=lambda x: len(x), reverse=True)
    for k in range(1,K+1):#, 2, 4, 8, 16]:
            print 'k - thres: ' + str(k) + ' - ' + str(thres)
            print 'overall time: ' + str(algo2.timer)
            print '(p,r,f1): ' + str(f1(allRecords, entities[:k], k, 
                                                 target=targetEntities))
    

def main(argv=None):
    
    testDBLP4()
#     testCora()
#     testSpotSigs()
    return
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))

    