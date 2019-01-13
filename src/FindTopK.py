'''
Created on Mar 15, 2016

@author: bill
'''
import sys

from utilities import *

def main(argv=None):
    #first read the ground truth and prior prob on edges
    gold = readGold(argv[1])
    prior = readEdges(argv[2])
    k = int(argv[3])
    method = argv[4]
    
    pickPair = maxFirst
    if method == 'EpsilonResolve':
        pickPair = epsilonRes
    
    topKGold = topK(gold, k)
    
    #run for B questions
    answers = dict()
    B = 100
    while B > 0:
        B -= 1
        #find the current topK
        tK, epsilon, _ = computeEpsilon(prior, k)
        #what is the accuracy now
        p, r = precisionRecall(tK, topKGold)
        
        print str(B) + ' - p: ' + str(p) + ', r: ' + str(r) + '\n'
        print 'epsilon: ' + str(epsilon)
        
        #what if we use only the human verified edges
        tK = onlyHuman(answers, k)
        p, r = precisionRecall(tK, topKGold)
        print 'Only human-verified, p: ' + str(p) + ', r: ' + str(r) + '\n'
        
        #pick the next pair to ask a question
        pair = pickPair(prior, answers, k, epsilon)
        
        print 'pair: ' + str(pair[0]) + ', prob: ' + str(pair[1]) + '\n'
        
        #get answer
        answer = getAnswer(pair[0], gold)
        print 'answer: ' + str(answer) + '\n'
        
        answers[pair[0]] = answer
        prior[pair[0]] = float(answer)
        

if __name__ == '__main__':
    sys.exit(main(sys.argv))
