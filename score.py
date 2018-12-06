from __future__ import division
import math
import sys
import argparse
from itertools import izip



def main():

    # Handle arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('generatedFileName', help='Name/Path of generated MeSH code file')
    parser.add_argument('answerFileName', help = 'Name/Path of the correct MeSH code file')
    args = parser.parse_args()

    genFile = open(args.generatedFileName,'r')
    answerFile = open(args.answerFileName,'r')

    totalPrecision = 0
    totalRecall = 0
    totalLines = 0
    for genLine, answerLine in izip(genFile,answerFile):
        if genLine == None or answerLine == None:
            break
        genCodes = genLine.split()
        answerCodes = answerLine.split()
        codeDict = {}
        for code in answerCodes:
            codeDict[code] = True

        genCodeDict = {}
        correctCodes = 0
        for code in genCodes:
            if code in codeDict:
                correctCodes += 1
            genCodeDict[code] = True

        recalledCodes = 0
        for code in answerCodes:
            if code in genCodeDict:
                recalledCodes += 1


        precision = correctCodes / len(genCodes)
        recall = recalledCodes / len(answerCodes)


        # print precision
        totalPrecision += precision
        totalRecall += recall
        totalLines += 1

    averagePrecision = totalPrecision/totalLines
    averageRecall = totalRecall/totalLines

    fmeasure = 2 / ((1 / averagePrecision) + (1 / averageRecall))

    print "Average precision: " + str(averagePrecision)
    print "Average recall: " + str(averageRecall)
    print "F-measure of the averages: " + str(fmeasure)


if __name__ == '__main__':
    main()