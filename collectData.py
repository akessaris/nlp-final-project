# you need to install Biopython:
# pip install biopython

'''
Baseline code for querying PubMed using Biopython was adapted from the following source
https://gist.github.com/bonzanini/5a4c39e4c02502a8451d
'''

from Bio import Entrez
import random
from random import shuffle
import math
import json

# Uses good practices for retrieving large datasets
# https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large
def batched_search(query, num_files, batch_amount):
    num_processed = 0
    results = []
    while num_processed < num_files:
        handle = Entrez.esearch(db='pubmed', 
                                sort='relevance', 
                                retmax=batch_amount,
                                retstart = num_processed,
                                retmode='xml', 
                                term=query)
        results.append(Entrez.read(handle))
        num_processed += batch_amount
    return results

# Gets the data for all ids in id_list, pre-batch beforehand!
# id_list = array of ids ** standard return type of Entrez.read **
def fetch_details(id_list):
    # To pass the id_list first we must stringify it
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

def read_ids(ids):
    handle = Entrez.efetch(db='pubmed',
                       retmode='xml',
                       id=ids)
    return Entrez.read(handle)


if __name__ == '__main__':
    # Config variables
    num_files = 2
    batch_amount = 1
    Entrez.email = 'your.email@example.com'

    #Initialize file to write article number, mesh codes, and abstract to test/training/dev files
    test = open("data.test", "w")    
    train = open("data.train", "w")
    dev = open("data.dev", "w")
    files = [train, test, dev]


    # Store outputs
    outputs = []

    batchedResults = batched_search('skin cancer', num_files, batch_amount)
    batchedPapers = []
    for results in batchedResults:
        id_list = results['IdList']
        batchedPapers.append(fetch_details(id_list))

    # Store dictionary of abstracts
    abstract_dict = {}

    # Keep array of all mesh codes from corresponding articles
    meshList = []

    # Get JSON of MeSH Tree #'s
    with open("mesh.json", "r") as f:
      meshTreeDict = json.load(f)

    for papers in batchedPapers:
        for i, paper in enumerate(papers['PubmedArticle']):
            # Get PubMed ID of paper and article object
            pmid = int(str(paper['MedlineCitation']['PMID']))
            article = paper['MedlineCitation']['Article']

            #print article['Abstract']['AbstractText']

            # If article doesn't have an abstract, don't process it
            if 'Abstract' not in article:
                continue

            # Get all available text from the article
            text = article['Abstract']['AbstractText']
            textList = []
            for i in text:
              print (i)
              if i not in textList:
                textList.append(i)

            # Append all text to abstract
            abstract = ""
            for i in textList:
              abstract += i

            abstract_dict[pmid] = abstract

            # Store dictionary of mesh codes for paper
            meshDict = {}

            # Get mesh codes for paper (if paper contains them)
            if (u'MeshHeadingList' in paper['MedlineCitation'].keys()):
                meshHeadingList = paper['MedlineCitation'][u'MeshHeadingList'] 
            # Otherwise, move on to next paper
            else:
                continue 
            for head in enumerate(meshHeadingList):
                # Get MeSH unique ID and term
                meshUniqueID =  head[1][u'DescriptorName'].__getattribute__("attributes")[u'UI']
                meshTerm = head[1][u'DescriptorName'].__str__()
                  
                # Add mesh term and id to meshDict
                meshDict[meshUniqueID] = meshTerm

            #Append mesh dict to meshList
            meshList.append({pmid: meshDict})

            # Convert MeSH Unique IDs to Tree Number

            # Add to mesh codes to output string
            meshString = ""
            for code in meshDict.keys():
                # THIS IS WHERE YOU QUERY THE JSON FILE
                #if meshDict[code] in meshTreeDict:
                #  print meshDict[code]
                #  print meshTreeDict[meshDict[code]] 
                meshString += code + "\t"

            outputs.append(meshString + "\n" + abstract)
      
    # Shuffle the list of outputs before writing to files
    shuffle(outputs)

    # Split list of outputs in 80-10-10 ration for train, dev, test among files
    splitN = [int(math.ceil(len(outputs)*0.8)),int(math.ceil(len(outputs)*0.9)),len(outputs)]
    splitOutputs = []
    for i in range(3):
        start = 0
        if i > 0:
            start = splitN[i-1]
        end = splitN[i]
        splitOutputs.append(outputs[start:end])
 
    # Output strings to files
    for fileIndex, f in enumerate(files):
        for index, output in enumerate(splitOutputs[fileIndex]):
            text = (''.join(output)).encode('utf-8')
            f.write(str(index+1) + "\n" + text + "\n\n")
