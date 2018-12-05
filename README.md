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
