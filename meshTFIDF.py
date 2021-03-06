#!/cygdrive/c/python27/python.exe
# TODO: CHANGE THIS

from stop_list import closed_class_stop_words
from operator import itemgetter
import math
import sys
# from nltk.stem import *
import argparse
import gensim
from gensim.utils import simple_preprocess
from gensim.parsing.preprocessing import STOPWORDS
from nltk.stem import WordNetLemmatizer, SnowballStemmer
from nltk.stem.porter import *
import numpy as np
from sets import Set
import Queue

# import nltk
# nltk.download('wordnet')

def main():

    print "starting script..."

    # Handle arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('trainingFileName', help='Name/Path of the training file')
    parser.add_argument('testFileName', help = 'Name/Path of the testing file')
    parser.add_argument('--o', help = 'Optional name of the output file')
    args = parser.parse_args()
    outputFileName = args.o if args.o else args.testFileName + '.output'

    trainingFile = open(args.trainingFileName,'r')
    testFile = open(args.testFileName, 'r')
    outputFile = open(outputFileName,'w')

    stopWords = closed_class_stop_words

    print "Processing MeSH training data..."
    meshInfo = processTraining(trainingFile,stopWords)
    meshWords = meshInfo["meshWords"]
    meshTFC = meshInfo["termFrequencyCount"]
    meshTD = meshInfo["totalDocuments"]


    print "Processing test data..."
    abstractInfo = processDocuments(testFile,stopWords)
    abstractTerms = abstractInfo["termDocumentCount"]
    abstractTermFrequencyCount = abstractInfo["termFrequencyCount"]
    totalAbstractCount = abstractInfo["totalDocuments"]


    # print totalAbstractCount
    # return
    print "TFIDFing MeSH data..."
    meshVectors = getTFIDFVector(meshWords,meshTFC, meshTD)

    print "TFIDFing test data..."
    abstractVectors = getTFIDFVector(abstractTerms,abstractTermFrequencyCount,totalAbstractCount)
    

    numClosest = 1

    print "Finding best Mesh codes..."
    meshMatches = {}

    for abstractNum in abstractVectors:
        # if abstractNum > 10:
        #     break
        # print abstractNum

        meshMatches[abstractNum] = []
        abstractVector = abstractVectors[abstractNum]
        for meshCode in meshVectors:
            meshVector = meshVectors[meshCode]
            meshMatches[abstractNum].append((meshCode, inverseCosine(abstractVector,meshVector)))

        meshMatches[abstractNum] = sorted(meshMatches[abstractNum], key = itemgetter(1),reverse=True)


        # Two options here:
        # 1)  Take the numClosest best mesh codes
        # meshMatches[abstractNum] = meshMatches[abstractNum][0:numClosest]
        best = meshMatches[abstractNum][0][1]
        cutoff = best - best*0.16
        # 2)  Get all with a similarity greater than cutoff
        meshMatches[abstractNum] = getMeshWithCutoff(meshMatches[abstractNum],cutoff)


        meshList = [tup[0] for tup in meshMatches[abstractNum]]

        # This line was used when we wanted to find the best cutoff
        # meshList = [tup[0] + "\t" + str(tup[1]) for tup in meshMatches[abstractNum]]

        outputFile.write("\t".join(meshList))
        outputFile.write("\n")






#########################################
#           MATH FUNCTIONS              #
#########################################

def getTFIDFVector(termDocumentCount, termFrequencyCount, totalDocuments):
    vector = {}
    for documentName in termDocumentCount:
        vector[documentName] = {}
        document = termDocumentCount[documentName]
        for word in document:
            vector[documentName][word] = TFIDF(document[word], termFrequencyCount[word],totalDocuments)
    return vector

# TODO: Play with normalizing termCount
def TFIDF(termCount, documentFrequency, totalDocuments):
    return (termCount * math.log(totalDocuments/documentFrequency))

def inverseCosine(abstractVector, queryVector):
    numerator = 0
    abstractSum = 0
    querySum = 0
    denominator = 0
    for word in queryVector:
        numerator+= queryVector[word] * abstractVector.get(word,0)
        querySum += queryVector[word] * queryVector[word]
        abstractSum += abstractVector.get(word,0) * abstractVector.get(word,0)
    
    if numerator == 0:
        return 0

    denominator = math.sqrt(querySum*abstractSum)
    return numerator / denominator

def getMeshWithCutoff(meshTuples, cutoff):
    for i in range(len(meshTuples)):
        meshTup = meshTuples[i]
        # print meshTup
        if meshTup[1] < cutoff:
            if not i == 0:
                return meshTuples[0:i]
            else:
                return meshTuples[0:1]



#########################################
#             PREPROCESSING             #
#########################################

# Trial preprocessing
def lemmatize_stemming(text):
    stemmer = SnowballStemmer('english')
    return stemmer.stem(WordNetLemmatizer().lemmatize(text, pos='v'))

def preprocess(text):
    result = []
    for token in gensim.utils.simple_preprocess(text):
        if token not in gensim.parsing.preprocessing.STOPWORDS and len(token) > 3:
            result.append(lemmatize_stemming(token))
    return result

# disabled for now
def clean(word):
    # return word
    # word = word.lower()
    # word = word.replace(r"([\W\d])+","")
    return word





#########################################
#               PROCESSING              #
#########################################

# Returns a dictionary of the following form:
"""
{
    "meshWords" : meshWords,                      # meshWords[meshCode] gets you the vectorized word count for meshCode
    "termFrequencyCount" : termFrequencyCount,    # termFrequencyCount[word] gets you the number of different mesh codes word appears in
    "totalDocuments" : totalDocuments,            # number of mesh codes that showed up
}
"""
def processTraining(file,stopWords):
    # stemmer = PorterStemmer()

    # Worth playing with this value to be either every time the word comes up or how many different mesh codes its in.
    termDocumentCount = {}
    termFrequencyCount = {}
    totalDocuments = 0
    currentDocument = 0
    lastDocumentAddedForWord = {}

    meshWords = {}

    currentCodes = []

    currentMode = 'num'

    

    for line in file:


        # Starting a new abstract / Empty line case
        if len(line.split()) == 0:
            currentMode = 'num'
            continue

        # Handling docNum, currently just ignore
        if currentMode == 'num':
            currentMode = 'codes'
            continue

        # Processing MeSH codes
        if currentMode == 'codes':
            currentCodes = line.split()
            # usedCodes = {}
            for i in range(len(currentCodes)):
                code = currentCodes[i]
                currentCodes[i] =  code[0:11]
            # remove duplicates
            currentCodes = set(currentCodes)
            currentMode = 'abstract'
            continue

        # No preprocessing currently, might be worth looking into!
        # words = preprocess(line)

        # Handle abstracts
        words = line.split()

        bigram = ["",""]

        replaceWordInBigram = 0
        for word in words:
            word = clean(word)
            if word in stopWords:
                continue

            # Keep track of the previous words (done in this fashion so that we can easily implement trigrams if needed)
            useWordInBigram = (replaceWordInBigram + 1) % 2
            bigramWord = bigram[useWordInBigram] + " " + word
            bigram[replaceWordInBigram] = word
            replaceWordInBigram = useWordInBigram




            if word not in termFrequencyCount:
                termFrequencyCount[word] = 0
            if bigramWord not in termFrequencyCount:
                termFrequencyCount[bigramWord] = 0


            for mesh in currentCodes:

                # Handle appearance of new mesh code
                if mesh not in meshWords:
                    meshWords[mesh] = {}
                    totalDocuments+=1

                # A document is the vector associated with a Mesh code, not an abstract
                document = meshWords[mesh]

                # This is the first time the word appears associated with the mesh code
                if word not in document:
                    document[word] = 0
                    # We keep track of how many different mesh codes the word shows up in
                    termFrequencyCount[word]+=1

                if bigramWord not in document:
                    document[bigramWord] = 0
                    termFrequencyCount[bigramWord]+=1

                document[word] += 1
                document[bigramWord] += 1


    return {
        "meshWords" : meshWords,
        "termFrequencyCount" : termFrequencyCount,
        "totalDocuments" : totalDocuments, 
    }

def processDocuments(file,stopWords):
    # stemmer = PorterStemmer()
    termDocumentCount = {}
    termFrequencyCount = {}
    totalDocuments = 0
    currentDocument = 0

    for line in file:
        # words = preprocess(line)
        words = line.split()

        termDocumentCount[currentDocument] = {}
        document = termDocumentCount[currentDocument]
        replaceWordInBigram = 0
        bigram = ["",""]
        for word in words:
            word = clean(word)
            if word in stopWords:
                continue
            
            # Keep track of the previous words (done in this fashion so that we can easily implement trigrams if needed)
            useWordInBigram = (replaceWordInBigram + 1) % 2
            bigramWord = bigram[useWordInBigram] + " " + word
            bigram[replaceWordInBigram] = word
            replaceWordInBigram = useWordInBigram

            # init termfrequency
            if word not in termFrequencyCount:
                termFrequencyCount[word] = 0
            if bigramWord not in termFrequencyCount:
                termFrequencyCount[bigramWord] = 0

            # if a word has not shown up in this document before, increment termfrequencycount and init document
            if word not in document:
                document[word] = 0
                termFrequencyCount[word]+=1
            if bigramWord not in document:
                document[bigramWord] = 0
                termFrequencyCount[bigramWord]+=1

            # finally, increment document count
            document[word] += 1
            document[bigramWord]+=1

        totalDocuments += 1
        currentDocument += 1


    return {
        "termDocumentCount" : termDocumentCount,
        "termFrequencyCount" : termFrequencyCount,
        "totalDocuments" : totalDocuments, 
    }




if __name__ == '__main__':
    main()