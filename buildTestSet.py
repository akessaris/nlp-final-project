import argparse
def main():

    # Handle arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dataFileName', help='Name/Path of the data file')
    args = parser.parse_args()

    unlabeledFileName = args.dataFileName + '.unlabeled'
    answersFileName = args.dataFileName + '.answers'

    inputFile = open(args.dataFileName,'r')

    unlabeledFile = open(unlabeledFileName,'w')
    answersFile = open(answersFileName,'w')

    currentMode = 'num'
    lineNumber = 1
    for line in inputFile:

        # Starting a new abstract / Empty line case
        if len(line.split()) == 0:
            if currentMode == 'abstract':
                unlabeledFile.write("\n")

            currentMode = 'num'
            continue

        # Handling docNum, currently just ignore
        if currentMode == 'num':
            currentMode = 'codes'
            continue

        # Processing MeSH codes
        if currentMode == 'codes':
            codes = line.split()
            codestr = ""
            for i in range(len(codes)):
                code = codes[i]
                codestr += code[0:11] + "\t"
            answersFile.write(codestr + "\n")
            currentMode = 'abstract'
            continue

        # Continuing adding to the line until the next abstract starts
        unlabeledFile.write(line.rstrip() + " ")
        lineNumber+=1


if __name__ == '__main__':
    main()