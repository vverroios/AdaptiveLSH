'''
Created on Jun 30, 2016

@author: bill
'''

import sys

from random import shuffle
from collections import defaultdict

from lsh.utilities import shingleSet, applyPermute, loadCora

def andorScheme(x,y,records,field='title',
                shingleSize=5):
    #first we generate the set of all shingles
    shingles = list(shingleSet(records, field, 
                               shingleSize))
    
    #we produce x*y permutations
    permutations = [range(len(shingles))
                    for _ in range(x*y)]
    for permutation in permutations:
        shuffle(permutation)
        
    #we will produce y hashtables
    hashtables = list()
    for i in range(y):
        hashtable = defaultdict(list)
        for record in records:
            recHash = tuple([applyPermute(record, 
                                          permutations[j], 
                                          shingles, 
                                          field, 
                                          shingleSize)
                             for j in range(i*x,(i+1)*x)])
            hashtable[recHash].append(record)
        hashtables.append(hashtable)
        
    return hashtables

def runAndOr(allRecords, xRange, 
             yUpper, howMany, thres):
    shuffle(allRecords)
    
    topKEntities(allRecords[:howMany], 5, thres=thres)
    
    for x in xRange:
        print '|||||||||||||||||||'
        print '|||||||||||||||||||'
        print '|||||||X=' + str(x) + '|||||||'
        print '|||||||||||||||||||'
        print '|||||||||||||||||||'
        
        hashTables = andorScheme(x, yUpper, 
                                 allRecords[:howMany])
        for y in range(1, yUpper+1):
            print '|||||||||||||||||||'
            print '|||||||||||||||||||'
            print '|||||||y=' + str(y) + '|||||||'
            print '|||||||||||||||||||'
            print '|||||||||||||||||||'
            bucketStats(hashTables[:y], thres=thres)
    
def bucketStats(hashTables, thres=4):
    #sort the buckets from all hashtables
    allBuckets = [(bucket[0], bucket[1], i) 
                  for i in range(len(hashTables))
                  for bucket in hashTables[i].iteritems()]
    topBuckets = [bucket 
                  for bucket in allBuckets
                  if len(bucket[1]) >= thres]
    topBuckets = sorted(topBuckets, 
                        key=lambda b: len(b[1]), 
                        reverse=True)
    for bucket in topBuckets:
        entitiesInBucket(bucket, allBuckets)
        
def topKEntities(records, k, thres=None):
    entities = set([record['entity'] 
                    for record in records])
    perEntity = {entity: len([record 
                              for record in records
                              if record['entity'] == entity])
                 for entity in entities}
    if thres is not None:
        topEntities = [entity 
                       for entity in perEntity.iteritems()
                       if entity[1] >= thres]
        topEntities = sorted(topEntities, key=lambda x: x[1], 
                             reverse=True)
    else:
        topEntities = [entity 
                       for entity in perEntity.iteritems()]
        topEntities = sorted(topEntities, key=lambda x: x[1], 
                             reverse=True)
        topEntities = topEntities[:k]
        
    for entity in topEntities:
        print 'entity: ' + str(entity[0]) + ', ' + str(entity[1]) + ' records'
    
def entitiesInBucket(bucket, allBuckets):
    entities = set([record['entity']
                    for record in bucket[1]])
    print '------'
    print 'BucketHash: ' + str(bucket[0]) + ', size: ' + str(len(bucket[1])) 
    for entity in entities:
        #find the records for this entity
        records = [record
                   for record in bucket[1]
                   if record['entity'] == entity]
        print 'entity: ' + entity + ', ' + str(len(records)) + ' records'
        #find records of this entity in other buckets
        for otherBucket in allBuckets:
            
            #if different hashtable skip
            if otherBucket[2] != bucket[2]:
                continue
            
            #if same hashtable, same bucket skip
            if otherBucket[0] == bucket[0]:
                continue
            
            records = [record
                       for record in otherBucket[1]
                       if record['entity'] == entity]
            if len(records) == 0:
                continue
            
            print 'in bucket of size: ' + str(len(otherBucket[1])) + ', ' + str(len(records)) + ' records'

def main(argv=None):
    allRecords = loadCora()
    runAndOr(allRecords,
             xRange=range(1,17,15), 
             yUpper=1,
             howMany=2000, 
             thres=20)
        
if __name__ == '__main__':
    sys.exit(main(sys.argv))
     
