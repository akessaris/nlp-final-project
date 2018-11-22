# Clinical Article Classification and MeSH Tagging

## Group

* Alexander Kessaris
* Josh Mendelsohn
* Sabrina de Silva

## Problem Statement
Medical literature contains an abundance of useful information, but when it comes to classifying the exact data it holds, it can be overwhelming for any content management system. Thankfully, there are ways of grouping content under the standardized classification system known as the Medical Subject Headings (MeSH). Using vector representations of clinical documents based on MeSH terms would allow for rapid, automatic, and standardized classification and tagging of other clinical documents.


## Data Collection

We are using Biopython to extract the abstracts and MeSH codes of PubMed articles.

* collectData.py works as follows:
  * Specifiy how many articles we want to be returned (line 18)
  * Pass query to db (line 43)
  * From each PubMed article:
    * Parse in the abstract (line 60)
    * Parse in the mesh codes (line 68)
  * Shuffle articles and their corresponding metadata (line 95)
    * Since articles returned are sorted by relevance, this is to prevent one data set from containing radically more relevant data than another set
  * Output articles evenly to test, training, and development data sets (line 102)  
