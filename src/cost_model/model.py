'''
Created on Oct 13, 2016

@author: bill
'''

import time

from core.hashFunctions import jaccardSim
from itertools import combinations

def pairCost(records, fields, erRule=None):
    trials = 0
    t = time.time()
    for pair in combinations(records, 2):
        trials += 1
        if erRule is None:
            jaccardSim(pair, fields)
        else:
            erRule(pair)
            
    return (time.time()-t)/trials
