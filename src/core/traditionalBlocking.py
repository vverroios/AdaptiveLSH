'''
Created on Sep 27, 2017

@author: vasilis verroios
'''

from collections import defaultdict

def traditionalBlocking(recs, fields):
    #build an inverted index
    recIdToRed = {rec['original']['other']['internalId'] : rec
                  for rec in recs}
    tables = {field: defaultdict(set)
              for field in fields}
    print '|',
    count = 0
    for rec in recs:
        count += 1
        if count%(len(recs)/15) == 0:
            print '-',
        for field in fields:
            if field in rec['data']:
                for token in rec['data'][field]:
                    tables[field][token].add(rec['original']['other']['internalId'])
    print '|'
    
    print 'largest ten buckets per field:'
    for field in fields:
        print 'field ' + str(field) + ' buckets'
        print 'sizes: ' + str(sorted([(k, len(s))
                                      for k,s in tables[field].iteritems()],
                                     key=lambda x: x[1], 
                                     reverse=True)[:10])              
                    
