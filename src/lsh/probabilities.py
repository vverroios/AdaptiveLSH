'''
Created on Jun 20, 2016

@author: vasilis verroios
'''

import sys

import numpy as np

from math import cos, sin, pi, degrees, radians, factorial

from utilities import angle_between, allSameSign, loadCora
from random import random

RND_VECS = 10000

def cliqueRndJaccard(simThres, q):
    sig = 2*simThres/(simThres+1.0)
    unionTerms = [sig**i for i in range(q)]
    unionTerms = [(-1)**i*
                  unionTerms[i]*
                  factorial(q)/(factorial(q-i-1)*factorial(i+1))
                  for i in range(len(unionTerms))]
    return sig**(q-1) / sum(unionTerms)
    

def generateClique(psi, q, D):
    #first let's find how many 
    #rotations per dimension we need
    rotations = np.ones(D)
    where = 1
    while reduce(lambda x,y:x*y, rotations) < q:
        rotations[where] += 1
        where = (where+1)%D
        if where == 0:
            where += 1
    #just in case of q=2
    if q == 2:
        rotations[2] = 2.0
            
    vectors = list()
    where = np.ones(D)
    while len(vectors) < q:
        vectors.append(newVector(where, rotations, psi))
        advanceWhere(where, rotations)
        
    return vectors
        
def newVector(where, rotations, psi):
    return np.array([element(where, rotations, psi, i)
                     for i in range(len(where))])
    
def element(where, rotations, psi, i):
    if i == 0:
        return cos(psi/2)
    if rotations[i] == 1.0:
        return 0.0
    factors = [sin(angle(where, rotations, j)) 
               for j in range(1,i)]
    if lastElement(rotations, i):
        return sin(psi/2)*reduce(lambda x,y:x*y, factors, 1.0)
    else:
        return (sin(psi/2)*reduce(lambda x,y:x*y, factors, 1.0)
                *cos(angle(where, rotations, i)))
                
def lastElement(rotations, i):
    if i == len(rotations)-1:
        return True
    if rotations[i+1] == 1.0:
        return True
    return False

def angle(where, rotations, j):
    return (where[j]*2*pi) / rotations[j]

def advanceWhere(where, rotations):
    for i in range(1,len(where)):
        if where[i] == rotations[i]:
            where[i] = 1
        else:
            where[i] += 1
            return 
        
def pss(q, D):
    probs = list()
    for psiD in range(180):
        vectors = generateClique(radians(psiD), q, D)
        probs.append(probAllSameSign(vectors, D))
        
    return range(180), probs
        
def probAllSameSign(vectors, D):
    allSame = 0.0
    for _ in range(RND_VECS):
        if allSameSign(vectors, rnd_vec(D)):
            allSame += 1
            
    return allSame / RND_VECS

def rnd_vec(D):
    angles = [np.arccos(np.random.uniform(-1,1))
              if i == 0
              else
              np.random.uniform(0,np.pi*2)
              for i in range(D-1)]
    return [reduce(lambda x,y: x*y, 
                   [sin(angles[j])
                    for j in range(i)], 
                   1.0)
            if i == D-1
            else
            cos(angles[i])*
            reduce(lambda x,y: x*y, 
                   [sin(angles[j])
                    for j in range(i)], 
                   1.0)
            for i in range(D)]
    
def random_three_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return (x,y,z)
        
def main(argv=None):
    p = pss(8*4, 10)
    for i in range(len(p[0])):
        print ('Angle: ' + str(p[0][i]) +
               ', prob: ' + str(p[1][i]) + 
               ', (180-a)/180: ' + str((180.0-p[0][i])/180))
        
if __name__ == '__main__':
    sys.exit(main(sys.argv))
        
    
        
