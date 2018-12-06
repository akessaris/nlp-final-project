# Clinical Article Classification and MeSH Tagging

## Group
* Alexander Kessaris
* Josh Mendelsohn
* Sabrina de Silva

## Problem Statement
Medical literature contains an abundance of useful information, but when it comes to classifying the exact data it holds, it can be overwhelming for any content management system. Thankfully, there are ways of grouping content under the standardized classification system known as the Medical Subject Headings (MeSH). Using vector representations of clinical documents based on MeSH terms would allow for rapid, automatic, and standardized classification and tagging of other clinical documents.

## Data Collection
We are using Biopython to extract the abstracts and MeSH codes of PubMed articles.

* collectData.py
  * Specify how many articles we want to be returned
  * Pass query to db
  * From each PubMed article:
    * Parse in all available article text
    * Parse in the MeSH unique ids and convert to MeSH tree numbers
  * Shuffle articles and their corresponding metadata
    * Since articles returned are sorted by relevance, this is to prevent one data set from containing radically more relevant data than another set
  * Output articles in a 80-10-10 split between training, tset and development data sets respectively  

* getMeshTree.py
  * Creates JSON file with all MeSH terms and their corresponding MeSH tree numbers
  * Works with 2017 MeSH data set
  * Accessed by collectData to convert MeSH unique ids into tree numbers

* buildTestSet.py
  * Creates data.xxx.unlabeled and data.xxx.answers files to train and score with.
  * Run by python buildTestSet.py [dataFileName]

* meshTFIDF.py
  * Uses data found from collectData to train, and data.xxx.unlabeled from buildTestSet to test.
  * Builds a vector for each mesh code in training using TFIDF, words, and bigrams
  * Compares each document to each mesh code and sorts them based off likeliness
  * Outputs a data.xxx.unlabeled.output file which has the most relevant mesh codes
  * Current best "precision" (quotes because we are only taking the most relevant mesh) is 77%
  * Run by python meshTFIDF.py [trainingDataFileName] [unlabeledDataFileName]

* score.py
  * Currently only calculates average precision
  * Will be used to calculate f-measure soon as well
  * Run by python score.py [generatedMeshFileName] [answerMeshFileName]
