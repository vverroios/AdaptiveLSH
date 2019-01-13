'''
Created on Jun 13, 2016

@author: bill
'''

import sys

from utilities import topK, readEdges, readGold, \
computeEpsilon, compo, precisionRecallEdges, sample, \
topKClusters, fOne, createGold, createEdges

from lsh.utilities import loadCora

RUNS = 1

def main(argv=None):
    #e.g., ./datasets/Cora/ 1 100
    print testForSizeAndK(argv[1], int(argv[2]), int(argv[3])) 
    #,loadCora())

def testForSizeAndK(dataset, k, size, records=None):
    
    if records is None:
        #get the whole graph and entities
        gold = readGold(dataset+'gold.txt')
        #gold = readGold(dataset+'goldTest.txt')
        #gold = readGold(dataset+'goldTiny.txt')
        #prior = readEdges(dataset+'graphTiny.txt')
        prior = readEdges(dataset+'graph.txt')
    else:
        gold = createGold(records)
        prior = createEdges(records)
    
    precision = list()
    recall = list()
    f1 = list()
    epsilons = list()
    
    bestF1 = list()
    bestEpsilon = list()
    
    print '|',
    for _ in range(RUNS):
        print '-',
        #sample from the whole dataset
        if size > 0:
            goldS, priorS = sample(gold, prior, size)
        else:
            goldS = gold; priorS = prior

        #let's get the gold for this sample
        topKGold = topK(goldS, k)
        #let's compute the topK components for this sample
        topKComponents, epsilon, allClusterings = computeEpsilon(priorS, k, 
                                                                 whatToReturn=compo, 
                                                                 gold=goldS)
        epsilons.append(epsilon)
        
        p,r = precisionRecallEdges(topKGold, topKComponents)
        precision.append(p)
        recall.append(r)
        f1.append(fOne(p,r))
        
        #let's see what is the best we can get over all clusterings
        pr = [(precisionRecallEdges(topKGold, topKClusters(clustering[1], k)), 
               clustering[0], len(clustering[1]))  
              for clustering in allClusterings.iteritems()]
        f1s = [(fOne(x[0][0],x[0][1]), x[1])  
               for x in pr]
        
        #TODO: remove when done
        print '*********'
        print '*********'
        print '*********'
        for prerec in sorted(pr, 
                             key=lambda x: x[1]):
            print prerec
        
        bestF = max(f1s,key=lambda x:x[0])
        bestF1.append(bestF[0])
        bestEpsilon.append(bestF[1])
        
    print '|'
        
    return (sum(precision)/len(precision), 
            sum(recall)/len(recall), 
            sum(f1)/len(f1), 
            sum(bestF1)/len(bestF1))

if __name__ == '__main__':
    sys.exit(main(sys.argv))