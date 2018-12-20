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
    outputFile = open('graph.data','w')
    allCodes = []
    allCodes.append([])
    docNum = 0
    for genLine in genFile:
        line = genLine.split()
        if len(line) == 0:
            docNum+=1
            allCodes.append([])
        else:
            allCodes[docNum].append(line)

    answerLines = []
    for answerLine in answerFile:
        if answerLine == None:
            break
        answerLines.append(answerLine)



    bestScalar = 0.5
    bestFMeasure = 0
    finalPrecision = 0
    finalRecall = 0
    scalar = 0.5

    while scalar > 0:
        totalPrecision = 0
        totalRecall = 0
        totalLines = 0
        totalFMeasure= 0
        totalBestOne = 0
        curDoc = 0
        for answerLine in answerLines:
            if answerLine == None:
                break
            preCutoffGenTuples = allCodes[curDoc]
            best = float(preCutoffGenTuples[0][1])
            cutoff = best - best * scalar
            genTuples = getMeshWithCutoff(preCutoffGenTuples,cutoff)
            genCodes = [tup[0] for tup in genTuples]
            # print len(preCutoffGenTuples)- len(genCodes) 
            answerCodes = answerLine.split()



            codeDict = {}
            for code in answerCodes:
                codeDict[code] = True

            genCodeDict = {}
            correctCodes = 0
            firstCodeDone = False
            for code in genCodes:
                if not firstCodeDone:
                    if code in codeDict:
                        totalBestOne += 1
                if code in codeDict:
                    correctCodes += 1
                genCodeDict[code] = True

            recalledCodes = 0
            for code in answerCodes:
                if code in genCodeDict:
                    recalledCodes += 1


            precision = correctCodes / len(genCodes)
            recall = recalledCodes / len(answerCodes)
            fmeasure = 0
            if precision != 0 and recall != 0:
                fmeasure = 2 / ((1 / precision) + (1 / recall))

            # print precision
            totalPrecision += precision
            totalRecall += recall
            totalFMeasure += fmeasure
            totalLines += 1
            curDoc+=1

        averagePrecision = totalPrecision/totalLines
        averageRecall = totalRecall/totalLines
        averageFmeasure = totalFMeasure / totalLines
        avgFmeasure = 2 / ((1 / averagePrecision) + (1 / averageRecall))
        avgBestOne = totalBestOne/totalLines

        if averageFmeasure > bestFMeasure:
            bestScalar = scalar
            bestFMeasure = averageFmeasure
            finalRecall = averageRecall
            finalPrecision = averagePrecision
            # print scalar
        # print str(scalar) + " " + str(avgFmeasure)
        outputFile.write(str(scalar) + "\t" + str(averagePrecision) + "\t" + str(averageRecall) + "\t" + str(avgFmeasure)+"\n")
        scalar -= 0.001


    print "bestFMeasure = " + str(bestFMeasure)
    print "bestScalar = " + str(bestScalar)
    print "final recall = " + str(finalRecall)
    print "final precision = " + str(finalPrecision)
    # print "Average precision: " + str(averagePrecision)
    # print "Average recall: " + str(averageRecall)
    # print "Average f-measure: " + str(averageFmeasure)
    # print "F-measure of the averages: " + str(avgFmeasure)
    # print "Average precision of best-fit MeSH  " + str(avgBestOne)


def getMeshWithCutoff(meshTuples, cutoff):
    for i in range(len(meshTuples)):
        meshTup = meshTuples[i]
        if float(meshTup[1]) < cutoff:
            if not i == 0:
                return meshTuples[0:i]
            else:
                return meshTuples[0:1]

    return meshTuples


if __name__ == '__main__':
    main()