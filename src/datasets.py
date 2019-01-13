'''
Created on Mar 16, 2016

@author: vasilis verroios
'''

import sys

from utilities import readGold, readEdges

def main(argv=None):
    createGraph(argv[1], argv[2], argv[3])
    
def createGraph(filename, filename2, filename3):
    gold = readGold(filename)
    allItems = set([c for entity in gold.values()
                    for c in entity])
    
    graph = readEdges(filename2)
    graph = [p for p in graph.iteritems()
             if p[0][0] in allItems and p[0][1] in allItems]
    
    writeEdges(graph, filename3)
    
def writeEdges(graph, filename):
    with open(filename, 'w') as f:
        for p in graph:
            f.write(str(p[0][0]) + ',' + str(p[0][1]) + ',' + str(p[1]) + '\n')

if __name__ == '__main__':
    sys.exit(main(sys.argv))

