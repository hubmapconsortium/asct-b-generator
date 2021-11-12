#!/usr/bin/python3

import sys
import re

DEBUG = False

ASLevels = 10
CTLevels = 3
BGLevels = 10

asDict = {}

def outputStruct(parents):
    # we want to go from the top level structure to the bottom level
    parents.reverse()

    level = 0
    for key in parents:
        if DEBUG: print(key, end="\t")
        
        sType = asDict[key]['sType']
        if sType == "AS":
            # process anatomical structures
            if level == 0:
                outString = key + "\t" + asDict[key]['label'] + "\t" + asDict[key]['id']
            else:
                outString = "\t" + key + "\t" + asDict[key]['label'] + "\t" + asDict[key]['id']
            level += 1

        elif sType == "CT":
            # process cell types
            outString = ""
            while level < ASLevels:
                outString += "\t\t\t"
                level += 1
            outString += "\t" + key + "\t" + asDict[key]['label'] + "\t" + asDict[key]['id']
            level += 1
                
        elif sType == "BG":
            # process genes
            outString = ""
            while level < (ASLevels + CTLevels):
                outString += "\t\t\t"
                level += 1
            outString += "\t" + key + "\t" + asDict[key]['label'] + "\t" + asDict[key]['id']
            level += 1

        # outFile.write(key + "\t" + asDict[key]['label'] + "\t" + asDict[key]['id'])
        outFile.write(outString)
                
    if DEBUG: print()
    outFile.write("\n")

    # restore list order
    parents.reverse()

# this is very inefficient but optimized for easy inputs and assumed processing speed is irrelevant.
def getParents(leaf, parents):
    # loop through all stuctures and print any that contain this leaf as a child
    for key in asDict.keys():
        children = asDict[key]["children"]
        if not children:
            continue
        if leaf in children:
            parents.append(key)
            getParents(key, parents)

    # if at the top level, then output the list
    if leaf == "ovary": outputStruct(parents)

    # remove latest structure as we back up the recursion.
    parents.pop()
            
# execute script
if __name__ == "__main__":
    # open input file. argv[0] is the program name
    inFilename = sys.argv[1]
    inFile = open(inFilename, "r")
    
    # this will overwrite any existing file
    outFilename = sys.argv[2]
    outFile = open(outFilename, "w")
    
    contents = inFile.readlines()
    for line in contents:
        # the first line might be a header line. That's intrisically resolved as the header "shouldn't" match to any children
        name, label, reference, sType, children = re.split(r'\t', line.rstrip('\n'))

        # remove the quotes from the TSV file and convert from a string into a list
        if children:
            # doing the replace as a separate step to allow for ', ' and ',' in the source file
            children = children.replace('"','').replace(', ',',')
            kids = children.split(',')
        else:
            kids = []
            
        # add anatomical structure to our dictionary
        asDict.update({name:{"label":label,"id":reference,"sType":sType,"children":kids}})

    # all structures loaded so now we step though each structure and print out all leaves
    for key in asDict.keys():
        children = asDict[key]["children"]
        if not children:
            # no children so end point and need to print upstream structures
            parents = [key]
            getParents(key, parents)
            
    inFile.close()
    outFile.close()
