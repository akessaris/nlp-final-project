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

def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='4',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

if __name__ == '__main__':
    #Initialize file to write article number, mesh codes, and abstract to test/training/dev files
    test = open("data.test", "w")    
    train = open("data.train", "w")
    dev = open("data.dev", "w")
    files = [test, train, dev]

    # Store outputs
    outputs = []

    results = search('skin cancer')
    id_list = results['IdList']
    papers = fetch_details(id_list)

    # Store dictionary of abstracts
    abstract_dict = {}

    # Keep array of all mesh codes from corresponding articles
    meshList = []

    for i, paper in enumerate(papers['PubmedArticle']):
      # Get PubMed ID of paper and article object
      pmid = int(str(paper['MedlineCitation']['PMID']))
      article = paper['MedlineCitation']['Article']

      # If article doesn't have an abstract, don't process it
      if 'Abstract' not in article:
        continue

      # Get the abstract
      abstract = article['Abstract']['AbstractText'][0]
      abstract_dict[pmid] = abstract

      # Store dictionary of mesh codes for paper
      meshDict = {}

      # Used for output
      meshString = ""
      
      # Get mesh codes for paper (if paper contains them)
      if (u'MeshHeadingList' in paper['MedlineCitation'].keys()):
        meshHeadingList = paper['MedlineCitation'][u'MeshHeadingList'] 
      # Otherwise, move on to next paper
      else:
        continue 
      for  head in enumerate(meshHeadingList):
        qualifierNames = head[1]
        for j in enumerate(qualifierNames[u'QualifierName']):
          for k in enumerate(j):
            if (k[0] == 1):
              # Add mesh term and id to meshDict
              meshTerm = k[1]
              meshUniqueID =  meshTerm.__getattribute__("attributes")[u'UI']
              meshDict[meshTerm] = meshUniqueID

              # Add mesh term to meshString as well
              meshString += meshTerm + "\t" 

      #Append  mesh dict to meshList
      meshList.append({pmid: meshDict})

      # Add to output string
      outputs.append(meshString + "\n" + abstract)
      
    # Shuffle the list of outputs before writing to files
    shuffle(outputs)

    # Split list of outputs evenly among files
    n = int(math.ceil(len(outputs)/3.0))
    splitOutputs = [outputs[i * n:(i + 1) * n] for i in range((len(outputs) + n - 1) // n )]  
 
    # Output strings to files
    for fileIndex, f in enumerate(files):
      for index, output in enumerate(splitOutputs[fileIndex]):
         text = (''.join(output)).encode('utf-8')
         f.write(str(index+1) + "\n" + text)
