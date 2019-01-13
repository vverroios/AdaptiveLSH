'''
Created on Mar 15, 2016

@author: vasilis verroios
'''

from itertools import combinations
from random import sample as smpl
from collections import defaultdict

from lsh.utilities import jaccardStringSim

def precisionRecallEdges(gold, outcome):
    alltopKMatches = [sorted(edge)
                      for cluster in gold
                      for edge in combinations(cluster[1], 2)]
    allOutMatches = [sorted(edge)
                      for cluster in outcome
                      for edge in combinations(cluster, 2)]
    
    return (percentage(allOutMatches, alltopKMatches), 
            percentage(alltopKMatches, allOutMatches))
   
def fOne(p,r):
    if p==0.0 and r==0.0:
        return 0.0
    return 2*p*r/float(p+r)
     
def sample(gold, prior, size):
    entities = defaultdict(list)
    records = list()
    
    #let's pick the records first
    allRecords = [record 
                  for cluster in gold.values()
                  for record in cluster]
    
    while len(records) < size:
        #pick record
        record = smpl(allRecords,1)[0] 
        allRecords.remove(record)
        #find entity
        entity = [e[0] 
                  for e in gold.iteritems()
                  if record in e[1]][0]
        entities[entity].append(record)
        records.append(record)
        
    #now let's create the edges
    edges = {tuple(sorted(edge)): prior[tuple(sorted(edge))] 
             for edge in combinations(records, 2)}
    
    return entities, edges
    
def percentage(target, other):
    edges = [edge 
             for edge in target 
             if edge in other]
    return len(edges) / float(len(target))

def readGold(filename):
    gold = dict()
    with open(filename) as f:
        for line in f:
            items = line.split(',')
            gold[items[0]] = set(map(int,items[1:]))
            
    return gold

def readEdges(filename):
    prior = dict()
    with open(filename) as f:
        for line in f:
            items = line.split(',')
            records = map(int, items[:2])
            prior[(min(records), max(records))] = float(items[2])
            
    return prior

def precisionRecall(tK, topKGold):
    taken = set()
    overall_precision = 0.0
    
    for cluster in tK:
        #match it with one that is not taken yet
        max_p = 0
        assigned = None
        for goldC in topKGold:
            if goldC[0] in taken:
                continue
            #compute the precision
            p = len(cluster.intersection(goldC[1])) / float(len(cluster))
            if p > max_p:
                max_p = p
                assigned = goldC[0]
                
        overall_precision += max_p
        if assigned is not None:
            taken.add(assigned)
            
    return overall_precision, len(taken)/float(len(topKGold))
    
    

def epsilonRes(prior, answers, k, eps):
    edgesSorted = [p for p in prior.iteritems() 
                   if p[0] not in answers]
    edgesSorted = sorted(edgesSorted, 
                         key=lambda p: 1.0-p[1])
    
    best_eps = 1.0
    for edge in edgesSorted:
        #let's check the lower bound
        #if the epsilon of the edge is lower than the current epsilon
        #AND the probability of NO multiplied by the current epsilon 
        #is larger than the best we found so far
        # we can stop 
        if edge[1] < eps and (1.0 - edge[1])*eps >= best_eps:
            break
        
        #else we need to evaluate the current edge
        #epsilon in case of a YES
        pY = edge[1]
        prior[edge[0]] = 1.0
        tK, epsY, _ = computeEpsilon(prior, k)
        #epsilon in case of a NO
        prior[edge[0]] = 0.0
        tK, epsN, _ = computeEpsilon(prior, k)
        #restore the right prior for this edge
        prior[edge[0]] = pY
        
        overall_eps = pY*epsY + (1.0-pY)*epsN
        if overall_eps < best_eps:
            best_eps = overall_eps
            best_edge = edge 
        
    return best_edge
        
    

def maxFirst(prior, answers, k, eps):
    return max([p for p in prior.iteritems() 
                if p[0] not in answers], 
               key = lambda p: p[1])
    
def topK(gold, k):
    return sorted(gold.iteritems(), 
                  key=lambda c: len(c[1]),
                  reverse=True)[:k]
                  
def topKClusters(clusters, k):
    return sorted(clusters, 
                  key=lambda c: len(c),
                  reverse=True)[:k]
                  
def getAnswer(pair, gold):
    for cluster in gold.values():
        if pair[0] in cluster:
            if pair[1] in cluster:
                return 1
            else: 
                return 0
        if pair[1] in cluster:
            return 0
        
def onlyHuman(answers, k):
    yesAns = [a[0] 
              for a in answers.iteritems() 
              if a[1] == 1]
    allNodes = set()
    for ans in yesAns:
        allNodes.update(set(ans))
    
    clusters = [set([n]) for n in allNodes]
    
    for answer in yesAns:
        for cluster in clusters:
            if answer[0] in cluster:
                c1 = cluster
            if answer[1] in cluster:
                c2 = cluster
                
        clusters.remove(c1)
        clusters.remove(c2)
        clusters.append(c1.union(c2))
        
    return sorted(clusters, 
                  key=lambda c: len(c), 
                  reverse=True)[:k]
            
def clus(clusters,components):
    return clusters

def compo(clusters,components):
    return [findCompo(cluster, components) 
            for cluster in clusters]
    
def findCompo(cluster, components):
    for component in components:
        if len(component.intersection(cluster)) > 0:
            return component
        
def computeEpsilon(prior, k, whatToReturn=clus, 
                   logging=True, gold=None):
    #first compute all clusterings
    sortedPYedges = sorted(prior.iteritems(), 
                           key=lambda p: 1.0-p[1])
    
    allNodes = set()
    for ans in sortedPYedges:
        allNodes.update(set(ans[0]))
    
    clusterings = dict()
    clusters = [set([n]) for n in allNodes]
    
    if gold is not None:
        tK = topK(gold, k)
        topKgold = {tK[i][0]:{'rank': i, 'size': len(tK[i][1])} 
                    for i in range(len(tK))}
    
    for ans in sortedPYedges:
        answer = ans[0]
        for cluster in clusters:
            if answer[0] in cluster:
                c1 = cluster
            if answer[1] in cluster:
                c2 = cluster
        
        if c1 != c2: 
            if logging:
                #the two clusters merge
                print (isMatch(ans, gold) + ': ' + str(ans) + 
                       ', merges clusters of sizes: ' +
                       str(len(c1))+','+str(len(c2))+' - ' +
                       entity(ans[0][0],gold,topKgold)+
                       ','+entity(ans[0][1],gold,topKgold))
                       
            clusters.remove(c1)
            clusters.remove(c2)
            clusters.append(c1.union(c2))
        else:
            if logging:
                #the two clusters are already merged
                print (isMatch(ans, gold) + ': ' + str(ans) + 
                       ', inside cluster of size: ' +
                       str(len(c1))+' - ' +
                       entity(ans[0][0],gold,topKgold)+
                       ','+entity(ans[0][1],gold,topKgold))
        
        clusterings[ans[1]] = list(clusters)
    
    # let's sort the edges
    sortedEdges = sorted(prior.iteritems(), 
                         key=lambda p: 1.0-p[1] 
                         if p[1]>0.5 else p[1])
    
    # we will start enabling edges until we have found the top-K
    l=0
    r=len(sortedEdges)
    
    while True:
        #probe the one in the middle
        ans=sortedEdges[(l+r)/2]
        
        epsilon = 1.0-ans[1] if ans[1]>0.5 else ans[1]
        if epsilon == 0:
            continue
        
        #let's see if for this epsilon it is guaranteed that topK is found
        topKClusters, isGuarant = isGuaranteed(sorted(clusterings.iteritems(),
                                                      key=lambda c: c[0],
                                                      reverse=True), 
                                               epsilon, k, 
                                               whatToReturn=whatToReturn)
        if isGuarant:
            r = (l+r)/2
        else:
            l = (l+r)/2
        
        if l == r-1:
            #TODO be more precise with the epsilon value you return
            return topKClusters, epsilon, clusterings

def entity(node, gold, topKgold):
    if gold is None:
        return 'unknown'
    else:
        for entities in gold.iteritems():
            if node in entities[1]:
                return topKstats(entities[0],topKgold) 
        
def topKstats(entity,topKgold):
    if entity in topKgold:
        return (entity+'(rank:'+ str(topKgold[entity]['rank']) + 
                ', size:' + str(topKgold[entity]['size'])+')')
    else:
        return entity+'(non top-K)' 

def isMatch(ans, gold):
    if gold is None:
        return 'unknown'
    else:
        if getAnswer(ans[0], gold) == 1:
            return '*M'
        else:
            return 'N-M'

def isGuaranteed(clusterings, epsilon, k, whatToReturn=clus):
    #get the components and clusters for this epsilon
    components = list([c for c in clusterings 
                       if c[0] <= epsilon][0][1])
    
    cs = [c for c in clusterings 
          if c[0] <= 1.0-epsilon][0]
    clusters = sorted(cs[1],
                      key=lambda c: len(c),
                      reverse=True)[:k]
                  
    #for the current state to be guaranteed
    # we want 
    #1) all clusters to be in different components
    #2) the k-th cluster size to be larger than the largest "reduced" component
    
    if not allDifferent(clusters, components):
        return whatToReturn(clusters, components), False
    
    compoSizes = [reducedSize(c, clusters) for c in components]
    
    if min([len(c) for c in clusters]) >= max(compoSizes):
        return whatToReturn(clusters, components), True
    else:
        return whatToReturn(clusters, components), False
    
def reducedSize(component, clusters):
    for cluster in clusters:
        if list(cluster)[0] in component:
            return len(component) - len(cluster)
        
    return len(component)

def allDifferent(clusters, components):
    for component in components:
        contains = 0
        for cluster in clusters:
            if len(component.intersection(cluster)) > 0:
                contains += 1
            if contains > 1:
                return False
            
    return True
    
def createEdges(records, shingleSize=5):
    prior = {(min(int(records[i]['other']['internalId']),
                  int(records[j]['other']['internalId'])),
              max(int(records[i]['other']['internalId']),
                  int(records[j]['other']['internalId']))):
             jaccardStringSim(records[i]['title'], 
                              records[j]['title'], 
                              shingleSize)
             for i in range(len(records))
             for j in range(i+1, len(records))}
    return prior
    
def createGold(records):
    allEntities = set([record['entity']
                       for record in records])
    gold = {entity: set([int(record['other']['internalId'])
                         for record in records
                         if record['entity'] == entity])
            for entity in allEntities}
    return gold
            
    
