'''
Created on Jul 13, 2016

@author: bill
'''

import re
import string

from random import shuffle, random
from numpy import array, dot, arccos
from math import degrees
from scipy.spatial.distance import cosine

from parsing.utils import MISSING 

CORATITLEAUTHORS = 0.7 #Note that we run Cora with 0.7 in the paper exps
CORAOTHERS = 0.2 #Note that we run Cora with 0.2 in the paper exps
SPOTSIGS = 0.4
DBLPTITLE = 0.5#0.7#0.4
DBLPAUTHORS = 0.5#0.3#0.2
DBLP2THRES = 0.5#0.5
#threshold for cosine distance in angle degrees
#the lower the angle the more likely the two 
#vectors are to refer to the same entity
IMAGESTHRES = 2#5#3#5#10

EMPTY = '##*#||EMPTY||#*##'

class sketches:
    """
    Sketches for cosine distance. Current implementation supports 
    only the case where each record consists of a single field, 
    i.e., a vector
    """
    
    def __init__(self, records, x, y):
        """
        Creates a set of x*y sketching functions, that can be invoked 
        as self.h(i,j,r); for the (i,j)-th function on record r.
        
        Arguments:
        records -- the list of records. Each record is a dict. 
        Here we are interested in the 'data' key, which contains 
        a dict: for each field name there is a vector 
        (e.g., RGB histogram for images)
        x -- the functions will be used on a x-AND_y-OR scheme
        y -- the functions will be used on a x-AND_y-OR scheme 
        """
        
        field = records[0]['data'].keys()[0]
        dim = len(records[0]['data'][field])
        self._functions = dict()
        
        for j in range(y):
            for i in range(x):
                self._functions[(i,j)] = {'rndV': self._rndVector(dim),
                                          'field': field}
                
    def _rndVector(self, dim):
        return array([1 if random() < 0.5 else -1 
                      for _ in range(dim)])
        
    def h(self, i, j, r):
        """
        applies function i,j on record r
        
        Arguments:
        i -- the i-th hash on the j-th table
        j -- the j-th table
        r -- the record
        """
        f = self._functions[(i,j)]
        rndV = f['rndV']
        field = f['field']
        recV = r['data'][field]
        
        return 1 if dot(recV, rndV) > 0.0 else -1

class minHash:
    """
    Min Hashing functions for a set of records, where each record 
    consists of fields and each field is a set of elements
    """
    
    def __init__(self, records, x, y, 
                 x2=None, andFields=None, 
                 fast=True):
        #Note original experiments we performed with fast=False
        #moreover, there are cases where the fast implementation 
        #is slower, after all. 
        """
        Creates a set of x*y min hash functions, that can be invoked 
        as self.h(i,j,r); for the (i,j)-th function on record r.
        
        Arguments:
        records -- the list of records. Each record is a dict. 
        Here we are interested in the 'data' key, which contains 
        a dict: for each field name there is a set of elements 
        (e.g., title shingles for that records)
        we use a uniform aggregation on the similarity of different fields. 
        For example, for fields title and authors 
        the overall similarity is 0.5*title_sim + 0.5*authors_sim
        x -- the functions will be used on a x-AND_y-OR scheme
        y -- the functions will be used on a x-AND_y-OR scheme 
        x2 -- in case of a (x+x2)-AND_y-OR scheme
        andFields -- fields used in case of a (x+x2)-AND_y-OR scheme
        """
        
        fields = andFields if andFields is not None else [records[0]['data'].keys()]
        self._fast=fast
        self._preprocessRecords(records, [field 
                                          for fieldList in fields
                                          for field in fieldList])
        self._functions = dict()
        
        
        if x2 is not None:
            whichField2 = 0
            for j in range(y):
                for i in range(x+x2):
                    whichField = 0 if i < x else 1
                    whichField2 = 0  if i == x else whichField2
                    self._functions[(i,j)] = {'perm': self._pickPermutation(fields[whichField][whichField2]),
                                              'field': fields[whichField][whichField2]}
                    
                    if self._fast:
                        self._functions[(i,j)]['invIndex'] = self._createInvIndex(self._functions[(i,j)]['perm'], 
                                                                                  self._functions[(i,j)]['field'])
                    whichField2 = (whichField2+1)%len(fields[whichField])
            return
                    
        whichField = 0
        for j in range(y):
            for i in range(x):
                self._functions[(i,j)] = {'perm': self._pickPermutation(fields[0][whichField]),
                                          'field': fields[0][whichField]}
                
                if self._fast:
                    self._functions[(i,j)]['invIndex'] = self._createInvIndex(self._functions[(i,j)]['perm'], 
                                                                              self._functions[(i,j)]['field'])
                #we are using a uniform aggregation of similarities across all fields 
                whichField = (whichField+1)%len(fields[0])
                
    def _preprocessRecords(self, records, fields):
        self._globals = {field: list(set([element 
                                   for record in records
                                   for element in record['data'][field]])) 
                         for field in fields}
        if self._fast:
            tmpIndex = {field: {self._globals[field][i] : i
                                for i in range(len(self._globals[field]))}
                        for field in fields}
            for record in records:
                record['coding'] = {field: [tmpIndex[field][element]
                                            for element in record['data'][field]]
                                    for field in fields}
        for field, value in self._globals.iteritems():
            print ('total elements for field: ' + str(field) 
                   + ' is: ' + str(len(value)))
        
    def _pickPermutation(self, field):
        permutation = range(len(self._globals[field]))
        shuffle(permutation)
        return permutation
    
    #old fast implementation
    def _createInvIndexOld(self, permutation, field):
        return {self._globals[field][permutation[i]]: i
                for i in range(len(permutation))}
    
    def _createInvIndex(self, permutation, field):
        tmpIndex = {permutation[i]: i
                    for i in range(len(permutation))}
        return [tmpIndex[i]
                for i in range(len(permutation))]
    
    def h(self, i, j, r):
        """
        applies function i,j on record r
        
        Arguments:
        i -- the i-th hash on the j-th table
        j -- the j-th table
        r -- the record
        """
        f = self._functions[(i,j)]
        perm = f['perm']
        field = f['field']
        allElements = self._globals[field]
        if self._fast:
            recElements = r['coding'][field]
        else:
            recElements = r['data'][field] 
        
        if len(recElements) == 0:
            return EMPTY+str(int(1000*random()))
        
        if self._fast:
            invIndex = f['invIndex']
            scores = [invIndex[rElement]
                      for rElement in recElements]
            whichOne = scores.index(min(scores))
            elemToReturn = allElements[recElements[whichOne]]
            return elemToReturn
        else:
            for k in perm:
                if allElements[k] in recElements:
                    return allElements[k]
            
            
def shingling(records, fields, processingFunction=None):
    """
    applies shingling on a set of records
    
    Arguments: 
    records -- each record is a dictionary: field name / text value
    fields -- the fields to be processed
    processingFunction -- the function that processes each field
    """
    
    if processingFunction is None:
        processingFunction = _processField
        
    return [{'original': record,
             'data': {field: processingFunction(field, record[field], record)
                      for field in fields}
             }
            for record in records]
    
def allArticle(field, data, fullrecord=None):
    if field == 'spotsigs':
        return spotsigs(field, data, fullrecord)
    elif field.startswith('shingles'):
        return shingles(field, data, fullrecord, 
                        ssize=int(field[len('shingles'):len('shingles')+1]))
    elif field.startswith('exact'):
        return exactLength(data, wsize=int(field[len('exact'):]))
    raise ValueError('no processing function for field ' + str(field))

def exactLength(data, wsize=3):
    return set([term.lower() 
                for term in re.findall(r"[\w']+", data)
                if len(term) == wsize])

def shingles(field, data, fullrecord=None, ssize=5):
    if ssize == 0:
        return random()
    terms = ' '.join([term.lower() 
                      for term in re.findall(r"[\w']+", data)])
    return set([terms[i:i+ssize]
                for i in range(len(terms)+1-ssize)])
    
def spotsigs(field, data, fullrecord=None):
    #we are using the following stop words
    sw = ['a', 'there', 'was', 'said', 'the', 'is']
    allStopWords = ['a', 'about', 'above', 'after', 
                    'again', 'against', 'all', 'am', 
                    'an', 'and', 'any', 'are', "aren't", 
                    'as', 'at', 'be', 'because', 'been', 
                    'before', 'being', 'below', 'between', 
                    'both', 'but', 'by', "can't", 'cannot', 
                    'could', "couldn't", 'did', "didn't", 
                    'do', 'does', "doesn't", 'doing', 
                    "don't", 'down', 'during', 'each', 'few', 
                    'for', 'from', 'further', 'had', "hadn't", 
                    'has', "hasn't", 'have', "haven't", 'having', 
                    'he', "he'd", "he'll", "he's", 'her', 'here', 
                    "here's", 'hers', 'herself', 'him', 'himself', 
                    'his', 'how', "how's", 'i', "i'd", "i'll", "i'm", 
                    "i've", 'if', 'in', 'into', 'is', "isn't", 'it', 
                    "it's", 'its', 'itself', "let's", 'me', 'more', 
                    'most', "mustn't", 'my', 'myself', 'no', 'nor', 
                    'not', 'of', 'off', 'on', 'once', 'only', 'or', 
                    'other', 'ought', 'our', 'ours', 'ourselves', 
                    'out', 'over', 'own', 'same', "shan't", 'she', 
                    "she'd", "she'll", "she's", 'should', "shouldn't", 
                    'so', 'some', 'such', 'than', 'that', "that's", 
                    'the', 'their', 'theirs', 'them', 'themselves', 
                    'then', 'there', "there's", 'these', 'they', 
                    "they'd", "they'll", "they're", "they've", 
                    'this', 'those', 'through', 'to', 'too', 'under', 
                    'until', 'up', 'very', 'was', "wasn't", 'we', 
                    "we'd", "we'll", "we're", "we've", 'were', 
                    "weren't", 'what', "what's", 'when', "when's", 
                    'where', "where's", 'which', 'while', 'who', 
                    "who's", 'whom', 'why', "why's", 'with', "won't", 
                    'would', "wouldn't", 'you', "you'd", "you'll", 
                    "you're", "you've", 'your', 'yours', 'yourself', 
                    'yourselves']
    distance = 2
    chain_len = 3
    delim = ':'
    
    terms = [term.lower() 
             for term in re.findall(r"[\w']+", data)]
    spotsignatures = list()
    
    for word in sw:
        for i in range(len(terms)):
            if terms[i] == word:
                #we start a signature
                sigs = list()
                at = distance
                while len(sigs) < chain_len:
                    if i+at >= len(terms):
                        break
                    if terms[i+at] not in allStopWords:
                        sigs.append(terms[i+at])
                        at += distance
                    else:
                        at += 1
                        
                if len(sigs) > 0:
                    spotsignatures.append(word+delim+(delim.join(sigs)))
                    
    return set(spotsignatures)
    
def advancedProcessField(field, data, fullrecord):
    if field == 'other':
        #get the venue
        venue = fullrecord['published']
        #get the pages 
        pages = None
        if 'pages' in fullrecord['other']:
            pages = fullrecord['other']['pages']
        volume = None
        if 'volume' in fullrecord['other']:
            volume = fullrecord['other']['volume']
        #is there is no info really
        if (venue == MISSING and 
            pages is None and 
            volume is None):
            return set([venue])
        
        termSet = set()
        #process venue 
        if venue != MISSING:
            termSet.update(_processVenue(venue))
        
        #process pages
        if pages is not None:
            non_decimal = re.compile(r'[^\d.]+')
            termSet.update([number 
                            for number in [non_decimal.sub('', term)
                                           for term in 
                                           re.split('\s|[,.-]', pages)]
                            if len(number) >= 1])
        #process the volume
        if volume is not None:
            non_decimal = re.compile(r'[a-z]+')
            termSet.update([number 
                            for number in [non_decimal.sub('', term)
                                           for term in 
                                           re.split('\s|[,.-]', volume)]
                            if len(number) >= 1])
        
        #if there is anything to return
        if len(termSet) == 0:
            return set([MISSING])
        else:
            return termSet
            
        
    return _processField(field, data)

def advancedProcessFieldDBLP(field, data, fullrecord):
    if field == 'other':
        termsSet = _processVenue(fullrecord['venue'])
        if 'year' in fullrecord:
            years = [int(year) 
                     for year in re.findall(r'[12]\d{3}', 
                                            fullrecord['year'])]
            termsSet.update(set([year + i 
                                 for year in years
                                 for i in [-2, -1, 0, 1, 2]]))
        
        #if there is anything to return
        if len(termsSet) == 0:
            return set([MISSING])
        else:
            return termsSet
            
    return _processField(field, data)

def _processVenue(venue):
    full_terms = [term.lower() for term in 
                  re.split('\s|[,.]', venue)
                  if len(term) >= 1
                  and term.lower() not in ['in', 'proceedings', 
                                           'and', 'on', 
                                           'of', 'the']]
    
    if len(full_terms) == 0:
        return set([MISSING])
    
    acro = ''.join([term[0] 
                    for term in full_terms])
    shingleSize = 2
    termSet = set([acro[i:i+shingleSize] 
                   for i in range(len(acro)-shingleSize+1)])
    termSet.update(set(full_terms))
    return termSet

def DBLPprocessV4(field, data, fullrecord=None):
    if field == 'title':
        tokens = [token.lower()
                  for token in re.findall(r"[\w']+", data)
                  if len(token) > 4 
                  and token.isalnum()]
        return set(tokens)
    if field == 'venue':
        allUpperTerms = [term 
                         for term in 
                         re.findall(r"[\w']+", data)
                         if term.isupper() and 
                         len(term) > 2 and 
                         (term != 'ACM' and term != 'IEEE')]
        ssize = 2
        shingles = [term[i:i+ssize] 
                    for term in allUpperTerms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 3
        shingles.extend([term[i:i+ssize] 
                         for term in allUpperTerms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        shingles.extend([term[:i]+
                         ('' if i == (len(term)-1) 
                          else term[(i+1):]) 
                         for term in allUpperTerms
                         for i in range(len(term))])
        return set(shingles)
    if field == 'authors':
        return set([lastName.lower()
                    for lastName in 
                    re.split('\s|[,.]', data)
                    if len(lastName) > 1
                    and lastName != 'and'])
    if field == 'year':
        years = [int(year) 
                 for year in
                 re.findall(r'[12]\d{3}', data)]
        return set([year + i
                    for year in years
                    for i in [-2, -1, 0, 1, 2]])

def DBLPprocessV3(field, data, fullrecord=None):
    if field == 'title':
        tokens = [token.lower()
                  for token in re.findall(r"[\w']+", data)
                  if len(token) > 4 
                  and token.isalnum()]
        return set(tokens)
    if field == 'venue':
        allUpperTerms = [term 
                         for term in 
                         re.findall(r"[\w']+", data)
                         if term.isupper() and 
                         len(term) > 2 and 
                         (term != 'ACM' and term != 'IEEE')]
        ssize = 2
        shingles = [term[i:i+ssize] 
                    for term in allUpperTerms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 3
        shingles.extend([term[i:i+ssize] 
                         for term in allUpperTerms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        shingles.extend([term[:i]+
                         ('' if i == (len(term)-1) 
                          else term[(i+1):]) 
                         for term in allUpperTerms
                         for i in range(len(term))])
        return set(shingles)
    if field == 'authors':
        terms = [term.lower() 
                 for term in 
                 re.findall(r"[\w']+", data)
                 if (not term.isupper()) 
                 and term != 'and'
                 and len(term) > 1]
        ssize = 3
        shingles = [term[i:i+ssize] 
                    for term in terms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 5
        shingles.extend([term[i:i+ssize] 
                         for term in terms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        return set(shingles)
    if field == 'year':
        years = [int(year) 
                 for year in
                 re.findall(r'[12]\d{3}', data)]
        return set([year + i
                    for year in years
                    for i in [-2, -1, 0, 1, 2]])

def DBLPprocess(field, data, fullrecord=None):
    if field == 'title':
        text = ''.join(re.findall(r"[\w']+", data))
        ssize = 5
        shingles = [text.lower()[i:i+ssize] 
                    for i in range(len(text)-ssize+1)]
        return set(shingles)
    if field == 'venue':
        allUpperTerms = [term 
                         for term in 
                         re.findall(r"[\w']+", data)
                         if term.isupper() and 
                         len(term) > 2 and 
                         (term != 'ACM' and term != 'IEEE')]
        ssize = 2
        shingles = [term[i:i+ssize] 
                    for term in allUpperTerms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 3
        shingles.extend([term[i:i+ssize] 
                         for term in allUpperTerms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        shingles.extend([term[:i]+
                         ('' if i == (len(term)-1) 
                          else term[(i+1):]) 
                         for term in allUpperTerms
                         for i in range(len(term))])
        return set(shingles)
    if field == 'authors':
        terms = [term.lower() 
                 for term in 
                 re.findall(r"[\w']+", data)
                 if (not term.isupper()) 
                 and term != 'and'
                 and len(term) > 1]
        ssize = 3
        shingles = [term[i:i+ssize] 
                    for term in terms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 5
        shingles.extend([term[i:i+ssize] 
                         for term in terms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        return set(shingles)
    if field == 'year':
        years = [int(year) 
                 for year in
                 re.findall(r'[12]\d{3}', data)]
        return set([year + i
                    for year in years
                    for i in [-2, -1, 0, 1, 2]])
        
def DBLPprocessV2(field, data, fullrecord=None):
    if field == 'title':
        text = ''.join(re.findall(r"[\w']+", data))
        ssize = 5
        shingles = [text.lower()[i:i+ssize] 
                    for i in range(len(text)-ssize+1)]
        if len(shingles) == 0:
            return set([MISSING])
        return set(shingles)
    if field == 'venue':
        allUpperTerms = [term 
                         for term in 
                         re.findall(r"[\w']+", data)
                         if term.isupper() and 
                         len(term) > 2 and 
                         (term != 'ACM' and term != 'IEEE')]
        ssize = 2
        shingles = [term[i:i+ssize] 
                    for term in allUpperTerms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 3
        shingles.extend([term[i:i+ssize] 
                         for term in allUpperTerms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        shingles.extend([term[:i]+
                         ('' if i == (len(term)-1) 
                          else term[(i+1):]) 
                         for term in allUpperTerms
                         for i in range(len(term))])
        if len(shingles) == 0:
            return set([MISSING])
        return set(shingles)
    if field == 'authors':
        terms = [term.lower() 
                 for term in 
                 re.findall(r"[\w']+", data)
                 if (not term.isupper()) 
                 and term != 'and'
                 and len(term) > 1]
        ssize = 3
        shingles = [term[i:i+ssize] 
                    for term in terms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 5
        shingles.extend([term[i:i+ssize] 
                         for term in terms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        if len(shingles) == 0:
            return set([MISSING])
        return set(shingles)
    if field == 'year':
        years = [int(year) 
                 for year in
                 re.findall(r'[12]\d{3}', data)]
        if len(years) == 0:
            return set([MISSING])
        return set([year + i
                    for year in years
                    for i in [-1, 0, 1]])
        
def DBLPprocessV1(field, data, fullrecord=None):
    if field == 'title':
        text = ''.join(re.findall(r"[\w']+", data))
        ssize = 5
        shingles = [text.lower()[i:i+ssize] 
                    for i in range(len(text)-ssize+1)
                    if len(text) >= ssize]
        if len(shingles) == 0:
            return set([MISSING])
        return set(shingles)
    if field == 'authors':
        terms = [term.lower() 
                 for term in 
                 re.findall(r"[\w']+", data)
                 if (not term.isupper()) 
                 and term != 'and'
                 and len(term) > 1]
        ssize = 5
        shingles = [term[i:i+ssize] 
                    for term in terms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        ssize = 3
        shingles.extend([term[i:i+ssize] 
                         for term in terms
                         if len(term) >= ssize
                         for i in range(len(term)-ssize+1)])
        if len(shingles) == 0:
            return set([MISSING])
        return set(shingles)
    if field == 'other':
        data = fullrecord['venue']
        allUpperTerms = [term 
                         for term in 
                         re.findall(r"[\w']+", data)
                         if term.isupper() and 
                         len(term) > 2 and 
                         (term != 'ACM' and term != 'IEEE')]
        ssize = 3
        shingles = [term[i:i+ssize] 
                    for term in allUpperTerms
                    if len(term) >= ssize
                    for i in range(len(term)-ssize+1)]
        venue = set(shingles)
        
        data = fullrecord['year']
    
        years = [int(year) 
                 for year in
                 re.findall(r'[12]\d{3}', data)]
        years = set([year + i
                     for year in years
                     for i in [-1, 0, 1]])
        
        if len(venue) == 0 and len(years) == 0:
            return set([MISSING])
        else:
            return venue.union(years)
        
def _processField(field, data, fullrecord=None):
    
    if field == 'html':
        return
    
    if field == 'venue':
        return _processVenue(data)
    
    if field == 'author' or field == 'authors':
        return set([lastName.lower()
                    for lastName in 
                    re.split('\s|[,.]', data)
                    if len(lastName) > 1
                    and lastName != 'and'])
    
    text = ' '.join(data.lower().split())
    text = text.replace('-',' ').replace(':',' ').translate(string.maketrans("",""), string.punctuation)
    if field == 'title':
        shingleSize = 5
    elif field == 'published':
        shingleSize = 4
    return set([text[i:i+shingleSize] 
                for i in range(len(text)-shingleSize+1)])
    
def jaccardSim(pair, fields):
    if (len(pair[0]['data'][fields[0]]) == 0 
        and len(pair[1]['data'][fields[0]]) == 0): 
            return 0.0
    
    sims = [len(pair[0]['data'][field].intersection(pair[1]['data'][field])) 
            / float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
            for field in fields]
    
    return sum(sims) / len(sims)

def average(pair, fields):
    return (pair[0]['data'][fields[0]]+pair[1]['data'][fields[0]])/2.0 

def titleAuthorsSimANDothers(pair, thres1=CORATITLEAUTHORS, 
                             thres2=CORAOTHERS, mustMatch=None):
    sims = [len(pair[0]['data'][field].intersection(pair[1]['data'][field])) 
            / float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
            for field in ['title', 'author', 'other']]
    
    if ((sims[0]+sims[1])/2.0 >= thres1) and (sims[2] >= thres2):
        return True
    else:
        return False

def imagesCosineRecs(rec1, rec2):
    field = 'hist'
    return 1.0 - cosine(rec1['data'][field], rec2['data'][field])    

def imagesCosine(pair, thres=IMAGESTHRES, mustMatch=None):
    field = 'hist'
    cosdis = cosine(pair[0]['data'][field], 
                    pair[1]['data'][field])
    if degrees(arccos(1.0-cosdis)) <= thres:
        return True
    else:
        return False
    
def titleAuthorsSimANDothers2(pair, thres1=0.5, 
                              thres2=0.5, mustMatch=None):
    try:
        sims = [len(pair[0]['data'][field].intersection(pair[1]['data'][field])) 
                / float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
                for field in ['title', 'authors', 'year']]
    except:
        print pair
        raise
    
    if ((0.7*sims[0]+0.3*sims[1])/1.0 >= thres1):
        return True
    else:
        return False
    
def titleAuthorsSimANDvenue(pair, thres1=0.7, thres2=0.2):
    sims = [len(pair[0]['data'][field].intersection(pair[1]['data'][field])) 
            / float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
            for field in ['title', 'author', 'venue']]
    
    if ((sims[0]+sims[1])/2.0 >= thres1) and (sims[2] >= thres2):
        return True
    else:
        return False
    
def titleANDAuthorsSim(pair, thres1=DBLPTITLE, thres2=DBLPAUTHORS, mustMatch=False):
    denoms = [float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
              for field in ['title', 'authors']]
    
    if denoms[0] == 0.0 or denoms[1] == 0.0:
        if mustMatch:
            print 'not merged denoms: ' + str(denoms)
        return False
    
    sims = [len(pair[0]['data'][field].intersection(pair[1]['data'][field])) 
            / float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
            for field in ['title', 'authors']]
    
    if (sims[0] >= thres1) and (sims[1] >= thres2):
        return True
    else:
        if mustMatch:
            print 'not merged sims: ' + str(sims)
        return False
    
def titleANDAuthorsSim2(pair, mustMatch=False):
    sims = [len(pair[0]['data'][field].intersection(pair[1]['data'][field])) 
            / float(len(pair[0]['data'][field].union(pair[1]['data'][field])))
            for field in ['title', 'authors']]
    
    if (0.7*sims[0] + 0.3*sims[1]) >= 0.4:
        return True
    else:
        if mustMatch:
            print 'not merged sims: ' + str(sims)
        return False 
    
    
    
                
            
        
        
        
        
