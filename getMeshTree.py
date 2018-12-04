import re
import json
 
terms = {}
numbers = {}
 
meshFile = 'd2017.bin'
with open(meshFile, mode='rb') as file:
    mesh = file.readlines()
 
outputFile = open('mesh.json', 'w')
 
for line in mesh:
    meshTerm = re.search(b'MH = (.+)$', line)
    if meshTerm:
        term = meshTerm.group(1)
    meshNumber = re.search(b'MN = (.+)$', line)
    if meshNumber:
        number = meshNumber.group(1)
        numbers[number.decode('utf-8')] = term.decode('utf-8')
        if term in terms:
            terms[term] = terms[term] + ' ' + number.decode('utf-8')
        else:
            terms[term] = number.decode('utf-8')
 
meshData = {}


meshNumberList = []
meshTermList = terms.keys()
for term in meshTermList:
    item_list = terms[term].split(' ')
    for phrase in item_list:
      # Create term in dict if not already created
      if term not in meshData:
        meshData.setdefault(term, [])
      meshData[term].append(phrase) 

json.dump(meshData, outputFile)
