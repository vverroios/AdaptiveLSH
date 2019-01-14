'''
Created on Jan 13, 2019

@author: vasilis verroios
'''
import xml.etree.ElementTree as ET

STOP_WORDS = []
MISSING = '####'

def parse_file(filename, includeId=False, concat=True):
    #list of records
    record_list = list()
    tree = ET.parse(filename)
    root = tree.getroot()
    foundEntity = False
    foundTitle = False
    foundAuthor = False
    foundPublished = False
    #get all the records
    i = 0
    for record in root.findall('record'):
        other_fields = dict()
        for attribute in record.findall('attribute'):
            #get the name attribute
            attrName = attribute.find('name')
            if attrName.text == 'title':
                title = attribute.find('value').text
                foundTitle = True
            elif attrName.text == 'entity':
                entity = attribute.find('value').text
                foundEntity = True 
            elif attrName.text == 'author':
                author = attribute.find('value').text
                foundAuthor = True
            elif attrName.text == 'booktitle':
                published = attribute.find('value').text
                foundPublished = True
            elif attrName.text == 'tech':
                published = attribute.find('value').text
                foundPublished = True
            elif attrName.text == 'type':
                published = attribute.find('value').text
                foundPublished = True
            elif attrName.text == 'journal':
                published = attribute.find('value').text
                foundPublished = True
            elif attrName.text == 'internal id' and includeId:
                other_fields['internalId'] = attribute.find('value').text
            else:
                other_fields[attrName.text] = attribute.find('value').text
        #filter the published
        if foundPublished and concat:
            published = remove_punc(remove_stop(published))
        #fill in the fields of this record
        if not foundTitle or len(remove_punc(title)) < 5:
            title = MISSING
        if not foundAuthor or len(remove_punc(author)) < 5:
            author = MISSING
        if not foundPublished or len(remove_punc(published)) < 5:
            published = MISSING
        if not foundEntity:
            print "Corrupted Entity: " + str(i)
            exit -1
        
        record_list.append({'entity': entity, 'title': title, 
                            'author': author, 'recid': i,
                            'published': published, 
                            'other': other_fields})
        i += 1
        foundEntity = False
        foundTitle = False
        foundAuthor = False
        foundPublished = False
        
    print "Successfully processed " + str(i) + " records"
    
    return record_list
   
def remove_punc(s):
    return ''.join(e for e in s if e.isalnum())

def remove_stop(published):
    for word in STOP_WORDS:
        published = published.replace(word, '')
    return published
