'''
Created on Jun 20, 2016

@author: vasilis verroios
'''

import sys
import numpy as np
import scipy.integrate as integrate
import pickle

from collections import defaultdict
from os import listdir
from os.path import isfile, join, isdir
from HTMLParser import HTMLParser
from random import random, sample
from parsing.utils import parse_file, MISSING
from core.hashFunctions import shingling, spotsigs
from numpy.ma.core import floor, ceil

def loadCora(venueConcat=True):
    return parse_file('./datasets/Cora/CoraData.xml', 
                      includeId=True, 
                      concat=venueConcat)
    
class MLStripper(HTMLParser):
    def __init__(self):
        self.reset()
        self.fed = []
    def handle_data(self, d):
        self.fed.append(d)
    def get_data(self):
        return ''.join(self.fed)

def strip_tags(html):
    s = MLStripper()
    s.feed(html)
    return s.get_data()

def preprocessSpotSigs():
    allRecords = loadSpotSigs()
    records = shingling(allRecords, ['data'], spotsigs)
    with open('./datasets/spotsigs.pickle', 'wb') as handle:
        pickle.dump((allRecords, records), handle)

def loadPreProSpotSigs():
    with open('./datasets/spotsigs.pickle', 'rb') as handle:
        (allRecords, records) = pickle.load(handle)
    return (allRecords, records)

def printSpotSigsFiles():
    folder = './datasets/GoldSetOfDuplicates/'
    onlydirs = [f for f in listdir(folder) if isdir(join(folder,f))]
    for d in onlydirs:
        directory = folder+d
        for frec in [f for f in listdir(directory) 
                    if isfile(join(directory,f))]:
            if frec[0] == '.':
                continue
            print d+'/'+frec

def loadSpotSigs(gooseProcessed=True, folder=None, additionalFields=[]):
    if folder is None:
        if gooseProcessed:
            folder = './datasets/SpotSigs/'
        else:
            folder = './datasets/GoldSetOfDuplicates/'
        
    onlydirs = [f for f in listdir(folder) if isdir(join(folder,f))]
    records = list()
    currentId = 0
    for d in onlydirs:
        directory = folder+d
        for frec in [f for f in listdir(directory) 
                    if isfile(join(directory,f))]:
            if frec[0] == '.':
                continue
            with open(directory+'/'+frec) as f:
                currentId += 1
                records.append({'entity':d,
                                'other':{'internalId':currentId, 
                                         'fileName':frec},
                                'data': f.read().replace('\n', ' ') 
                                if gooseProcessed else  
                                strip_tags(unicode(f.read().replace('\n', ' '), 
                                                   errors='ignore'))})
                
    for record in records:
        for field in additionalFields:
            record[field] = record['data']
        
    return records

def computeANDHashSchemeGivenBudget(B, epsilon, dThr1, 
                                    dThr2, exact=False, 
                                    wConstraint=0, 
                                    uConstraint=0):
    #TODO: you also need a zConstraint cause u may end up with 
    #less tables from one sequence position to the next. However, 
    #this case is complicated because you wil have to consider 
    #solutions that do not utilize the full budget, in order not to 
    #add useless tables, before the budget becomes enough.
    
    #build (w,u,z)-scheme,
    #for (w,u) in [1,budget], z = floor(budget/w)
    best = None
    for s in range(1,B+1):
        for allocation in range(1,10):
            z = floor(B/float(s))
            sp = int(B - s*z)
            w,u = _allocate(s,allocation/10.0)
            if w < wConstraint or u < uConstraint:
                continue
            wp,up = _allocate(sp,allocation/10.0) if not exact else (0.0, 0.0)
            #if constraint is satisfied
            if _probANDSameBucket(w,u,z,dThr1,
                                  dThr2,wp=wp,
                                  up=up) >= 1.0-epsilon:
                obj = _areaAND(w,u,z,wp,up)
                if best is None or best[0] > obj:
                    best = [obj, w, u, z, wp, up]
    #return the best result found
    return best[1:]

def _allocate(amount,allocation):
    a1 = amount*allocation; a2 = amount*(1.0-allocation)
    if (a1-int(a1)) > (a2-int(a2)):
        a1 = ceil(a1); a2 = floor(a2)
    elif (a1-int(a1)) < (a2-int(a2)):
        a2 = ceil(a2); a1 = floor(a1)
    elif (a1-int(a1)) == 0.5:
        a1 = ceil(a1); a2 = floor(a2)  
        
    return (a1,a2)

def computeHashSchemeGivenBudget(B, epsilon, dThr, 
                                 exact=False, cosine=False):
    #build (w,z)-scheme,
    #for w in [1,budget], z = floor(budget/w), wp = budget - w*z
    best = None
    for w in range(1,B+1):
        z = floor(B/float(w))
        wp = int(B - w*z) if not exact else 0.0
        #if constraint is satisfied
        if _probSameBucket(w,z,dThr,wp=wp,
                           cosine=cosine) >= 1.0-epsilon:
            obj = _area(w,z,wp)
            if best is None or best[0] > obj:
                best = [obj, w, z, wp]
    #return the best result found
    return best[1:]
                
def _area(w,z,wp):
    return integrate.quad(lambda x: 1.0 - ((1.0-x**w)**z)*
                          (1.0 if wp==0.0 else 1.0-x**wp), 
                          0, 1.0)
    
def _areaAND(w,u,z,wp,up):
    return integrate.dblquad(lambda x1, x2: 
                             1.0 - ((1.0-(x1**w)*(x2**u))**z)*
                             (1.0 
                              if wp==0.0 and up==0.0 
                              else 1.0-(x1**wp)*(x2**up)), 
                             0, 1.0, 
                             lambda x: 0, 
                             lambda x: 1)

def _probANDSameBucket(w,u,z,dThr1,dThr2,
                       wp=0.0,up=0.0,px=None):
    if px is None:
        if wp == 0.0 and up == 0.0:
            pAdd = 1.0
        else:
            pAdd = 1.0 - (dThr1**wp)*(dThr2**up)
        return  1.0 - ((1.0-(dThr1**w)*(dThr2**u))**z)*pAdd

def _probSameBucket(w,z,dThr,wp=0.0,px=None, 
                    cosine=False):
    #in case of cosine distance the threshold is an angle in degrees
    if cosine:
        dThr = 1.0 - dThr/180.0
    
    if px is None:
        if wp == 0.0:
            pAdd = 1.0
        else:
            pAdd = 1.0 - dThr**wp
        return  1.0 - ((1.0-dThr**w)**z)*pAdd

def loadDBLPFiles(folder = './datasets/DBLP/'):
    acmRecs = readCSV(folder+'ACM.csv')
    dblp2Recs = readCSV(folder+'DBLP2.csv')#with acm
    dblpAcmMap = readCSV(folder+'DBLP-ACM_perfectMapping.csv')
    scholarRecs = readCSV(folder+'Scholar.csv')
    dblp1Recs = readCSV(folder+'DBLP1.csv')#with scholar
    dblpScholarMap = readCSV(folder+'DBLP-Scholar_perfectMapping.csv')
    
    dblpMap = _normalizeIds(acmRecs, dblp1Recs, dblp2Recs, 
                           scholarRecs, dblpAcmMap, 
                           dblpScholarMap)
    
    entities, stamp = _joinDBLPsets(acmRecs, dblp1Recs, dblp2Recs, 
                                    scholarRecs, dblpAcmMap, 
                                    dblpScholarMap, dblpMap)
    
    return (acmRecs, dblp2Recs, dblpAcmMap, scholarRecs, 
            dblp1Recs, dblpScholarMap, dblpMap, entities, 
            stamp) 

def loadDBLPcomplete(filename='./datasets/DBLP/complete/DBLP-72KP3.csv'):
    allRecs = readCSV(filename)
    for rec in allRecs:
        rec['other'] = {'internalId': rec['id']}
        rec['id'] = rec['originalID']
        del rec['originalID']
    
    gold = defaultdict(list)
    for rec in allRecs:
        gold[rec['entity']].append(rec['other']['internalId'])
    entities = [(entity, len(l))
                for entity, l in gold.iteritems()]
    
    topK = 40
    topOnes = [entity 
               for entity, _ in 
               sorted(entities, key=lambda x:x[1], 
                      reverse=True)][:topK]
    topOnesF = [(entity,size) 
                for entity, size in 
                sorted(entities, key=lambda x:x[1], 
                       reverse=True)][:topK]
    for i in range(len(topOnesF)):
        print 'entity top-'+str(i)+': ' +str(topOnesF[i])
        
    return allRecs, topOnes

def loadDBLP(folder = './datasets/DBLP/', 
             extraSmpl=0.0, fields=None):
    
    acmRecs, dblp2Recs, dblpAcmMap, scholarRecs, \
    dblp1Recs, dblpScholarMap, _, entities, \
    _ = loadDBLPFiles(folder=folder)
    
    extraRecs = list()
    bigOnes = [entity 
               for entity, size in entities
               if size > 7]
    topK = 20
    topOnes = [entity 
               for entity, size in 
               sorted(entities, key=lambda x:x[1], 
                      reverse=True)][:topK]
    topOnesF = [(entity,size) 
                for entity, size in 
                sorted(entities, key=lambda x:x[1], 
                       reverse=True)][:topK]
    for i in range(len(topOnesF)):
        print 'entity top-'+str(i)+': ' +str(topOnesF[i])
    
    return ([rec 
             for dset in [acmRecs, dblp1Recs, 
                          dblp2Recs, scholarRecs, 
                          extraRecs]
             for rec in dset
             if rec['entity'] in bigOnes
             or random() < extraSmpl], 
            (acmRecs, dblp1Recs, dblp2Recs, 
             dblpAcmMap, scholarRecs, dblpScholarMap),
            topOnes)
 
def _normalizeIds(acmRecs, dblp1Recs, dblp2Recs, 
                  scholarRecs, dblpAcmMap, dblpScholarMap):
    id1s = [rec['id']
            for rec in dblp1Recs]
    newMap = [{'idDBLP1':'DBLP1'+id1, 
               'idDBLP2':'DBLP2'+id1,}
              for id1 in id1s]
    
    tags = ['ACM', 'DBLP1', 'DBLP2', 'SCHOLAR']
    i = 0
    for dset in [acmRecs, dblp1Recs, 
                 dblp2Recs, scholarRecs]:
        for rec in dset:
            rec['id'] = tags[i] + rec['id']
        i += 1
    
    for t in dblpAcmMap:
        t['idDBLP'] = tags[2] + t['idDBLP']
        t['idACM'] = tags[0] + t['idACM']
    for t in dblpScholarMap:
        t['idDBLP'] = tags[1] + t['idDBLP']
        t['idScholar'] = tags[3] + t['idScholar']
        
    return newMap

def _joinDBLPsets(acmRecs, dblp1Recs, dblp2Recs, 
                  scholarRecs, dblpAcmMap, 
                  dblpScholarMap, dblpMap):
    stamp = 0
    for dset in [acmRecs, dblp1Recs, 
                 dblp2Recs, scholarRecs]:
        for rec in dset:
            stamp += 1
            rec['other'] = {'internalId': stamp}
    
    #for each record create a node
    nodes = {rec['id']: {'root': None, 'data': rec, 
                         'id': rec['id'], 'right': None, 
                         'first': None, 'last': None} 
             for dset in [acmRecs, dblp1Recs, 
                          dblp2Recs, scholarRecs]
             for rec in dset}
    
    #let's now process the mappings
    trees = list()
    for t in dblpAcmMap:
        _mergeTrees(t['idDBLP'], t['idACM'], nodes, trees)
    for t in dblpScholarMap:
        _mergeTrees(t['idDBLP'], t['idScholar'], nodes, trees)
    for t in dblpMap:
        _mergeTrees(t['idDBLP1'], t['idDBLP2'], nodes, trees)
        
    #find all the trees without roots
    nonSingles = 0
    entities = list()
    for tree in [t for t in trees
                 if t['root'] is None]:
        node = tree['first']
        allRecs = 0
        while node is not None:
            nodes[node['id']]['data']['entity'] = tree['id']
            node = node['right']
            allRecs += 1
            if allRecs > 100000:
                print node['id']
                print node['right']['id']
        entities.append((tree['id'], allRecs))
        nonSingles += allRecs
    
    #for all the nodes without an entity assigned to them
    singles = [n for n in nodes.values()
               if 'entity' not in n['data']] 
    
    print 'allRecs: ' + str(len(nodes))
    print 'singles: ' + str(len(singles))
    print 'nonSingles: ' + str(nonSingles)
    
    for node in singles:
        node['data']['entity'] = node['id']
        
    return entities, (stamp+1)
        
def _mergeTrees(id1, id2, nodes, trees):
    #find the roots of the two nodes
    root1 = _findRoot(nodes[id1]); root2 = _findRoot(nodes[id2])
    if root1['id'] == root2['id']:
        #already merged
        return
    
    first1 = root1['first'] if root1['data'] is None else root1
    first2 = root2['first'] if root2['data'] is None else root2
    last1 = root1['last'] if root1['data'] is None else root1
    last2 = root2['last'] if root2['data'] is None else root2
     
    #merge the two trees
    newRoot = {'root': None, 'data': None, 
               'id': root1['id'], 'right': None, 
               'first': first1, 'last': last2}
    #update the old root pointers
    root1['root'] = newRoot; root2['root'] = newRoot
    #update the node pointers
    last1['right'] = first2
    
    trees.append(newRoot)

def _findRoot(node):
    while node['root'] is not None:
        node = node['root']
        
    return node
 
def _someSimpleStats(dblpAcmMap, dblpScholarMap):
    id1s = [t['idDBLP'] for t in dblpAcmMap]
    id2s = [t['idACM'] for t in dblpAcmMap]
    id3s = [t['idDBLP'] for t in dblpScholarMap]
    id4s = [t['idScholar'] for t in dblpScholarMap]
    print 'id1s: ' + str(len(id1s)) + ', set: ' + str(len(set(id1s)))
    print 'id2s: ' + str(len(id2s)) + ', set: ' + str(len(set(id2s)))
    print 'id3s: ' + str(len(id3s)) + ', set: ' + str(len(set(id3s)))
    print 'id4s: ' + str(len(id4s)) + ', set: ' + str(len(set(id4s)))
    
    dupIds = [(idDup, len([1
                           for recId in id4s
                           if idDup == recId])) 
              for idDup in set(id4s)]
    print 'dupIds: ' + str('\n'.join(map(str,[t for t in dupIds
                                              if t[1] > 1])))
    
def _compareDBLPs(dblp1Recs, dblp2Recs):
    fields = ['id', 'title', 'authors', 
              'venue', 'year']
    id2Dic = {rec['id']: rec 
              for rec in dblp2Recs}
    id1Dic = {rec['id']: rec 
              for rec in dblp1Recs}
    matches = [rec['id'] 
               for rec in dblp1Recs
               if rec['id'] in id2Dic]
    absMatches = [recId
                  for recId in matches
                  if sum([1 
                          if id1Dic[recId][field] == id2Dic[recId][field]
                          else 0
                          for field in fields]) == len(fields)]
    
    print 'matches: ' + str(len(matches))
    print 'absMatches: ' + str(len(absMatches))

def readCSV(filename, delim='|'):
    headers = None
    recs = list()
    with open(filename, 'r') as f:
        for line in f:
            if headers is None:
                headers = [field.strip() 
                           for field in line.split(delim)]
            else:
                terms = line.split(delim)
                recs.append({headers[i]: (terms[i].strip()
                                          if len(terms[i].strip()) > 0 
                                          else MISSING)
                             for i in range(len(headers))})
    return recs
            

def allSameSign(vectors, hyperplane):
    """ Returns true if all vectors are on the same side of the hyperplane.  """
    first = np.sign(np.dot(vectors[0], hyperplane))
    for i in range(1,len(vectors)):
        if np.sign(np.dot(vectors[i], hyperplane)) != first:
            return False
    return True

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def jaccardStringSim(string1, string2, shingleSize):
    #convert to lowercase and replace 
    #all whitespaces with a single space
    string1 = ' '.join(string1.lower().split())
    string2 = ' '.join(string2.lower().split())
    
    #perform the shingling
    set1 = set([string1[i:i+shingleSize]
                for i in range(len(string1)-shingleSize+1)])
    set2 = set([string2[i:i+shingleSize]
                for i in range(len(string2)-shingleSize+1)])
    
    #jaccard sim of the two sets
    return float(len(set1.intersection(set2)))/len(set1.union(set2))

def shingleSet(records, field, shingleSize):
    text = [' '.join(record[field].lower().split())
            for record in records]
    allShingles = [set([t[i:i+shingleSize] 
                        for i in range(len(t)-shingleSize+1)])
                   for t in text]
    return set([s 
                for shingles in allShingles
                for s in shingles])
    
def applyPermute(record, permutation, allShingles, 
                 field, shingleSize):
    #TODO: this is a dummy way to apply min-hashing
    text = ' '.join(record[field].lower().split())
    recShingles = set([text[i:i+shingleSize] 
                       for i in range(len(text)-shingleSize+1)])
    for i in permutation:
        if allShingles[i] in recShingles:
            return allShingles[i]
    
def main(argv=None):
    printSpotSigsFiles()
    return
    
    print computeHashSchemeGivenBudget(200, 0.011, 0.4, exact=True)#0.916666666)
    B=100
    res = [0,0]
    while B <= 12800:
        res = computeANDHashSchemeGivenBudget(B, 0.02, 0.8, 0.7, 
                                              exact=True, 
                                              wConstraint=res[0],
                                              uConstraint=res[1])
        print res
        B *= 2
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
    
