'''
Created on Oct 19, 2016

@author: bill
'''
from goose import Goose
import os

def removeNonAlphaNumeric(s):
    return ''.join(ch for ch in s if ch.isalnum())

files = list()
with open('/Users/bill/python-goose/SpotSigsFileList.txt', 'r') as f:
    for line in f:
        files.append(line.strip())

directory = 'file:///Users/bill/Documents/fromAthens07Sep2012-selectionsOnly/fromOffice/workHomeAug2012/workspace/ERTopK/datasets/GoldSetOfDuplicates/'

datasetName = 'SpotSigs/'

for fileName in files[2126:2127]:
    try:
        g = Goose()
        print 'extracting article: ' + fileName
        article = g.extract(url=directory+fileName)
        print 'article extracted: ' + fileName
        terms = [removeNonAlphaNumeric(term.lower()) 
                 for term in 
                 article.cleaned_text.replace('\n', ' ').split()
                 if len(removeNonAlphaNumeric(term.lower())) > 0]
        if not os.path.exists(datasetName+str(fileName.split('/')[0])):
            os.makedirs(datasetName+str(fileName.split('/')[0]))
        with open(datasetName+fileName, 'w') as f:
            f.write(' '.join(terms))
        print 'Successfully processed: ' + fileName
    except Exception,e: 
        print 'Could not process file: ' + fileName
        print str(e)
        


######################
##OTHER THINGS I TRIED
######################

def check(test_str):
    import re
    #http://docs.python.org/library/re.html
    #re.search returns None if no position in the string matches the pattern
    #pattern to search for any character other then . a-z 0-9
    pattern = r'[^\.\,\\a-z0-9]'
    if re.search(pattern, test_str):
        #Character other then . a-z 0-9 was found
        return True
    else:
        #No character other then . a-z 0-9 was found
        return False

processedFiles = dict()
for fileName in files:
    try:
        g = Goose()
        article = g.extract(url=directory+fileName)
        processedFiles[fileName] = article.cleaned_text
        print 'Successfully processed: ' + fileName
    except:
        print 'Could not process file: ' + fileName

#unicodes involved
# codes = set()
terms = list()
for fileName in ['Smoking/original.htm',
                 'Smoking/s1.htm',
                 'Smoking/s2.htm',
                 'Smoking/s3.htm',
                 'Smoking/s4.htm',
                 'Smoking/s5.htm',
                 'Smoking/s6.htm',
                 'Smoking/s7.htm',
                 'Smoking/s8.htm']:
    try:
        article = g.extract(url=directory+fileName)
        terms.append([removeNonAlphaNumeric(term.lower()) 
                     for term in 
                     article.cleaned_text.replace('\n', ' ').split()
                     if len(removeNonAlphaNumeric(term.lower())) > 0])
#         codes.update(set([code.lower() 
#                           for code in 
#                           article.cleaned_text.replace('\n', ' ').split() 
#                           if check(code.lower())]))
    except:
        print 'could not process file: ' + fileName
    