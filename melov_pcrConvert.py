#!/usr/bin/env python

from optparse import OptionParser
import sys
import csv

# ------------------------------------------------------------------
def readCSV(inputFile = []):
    
    if inputFile:
        inputStream = open(inputFile)
    else:
        inputStream = sys.stdin
    
    # read the standard input stream, using csv reader
    line = 0
    reader = csv.reader(inputStream)
    for row in reader:
        line += 1
        
        # The PCR output header is 11 lines; do we need any information
        # from the header?
        if 'FAM-MGB Ct' in row:

            # The next line has the 'well' numbers, so fill the 'Well'
            # dictionary with a list of these numbers.
            row = reader.next()
            wellN = int(row[-1]) + 1

            # The next line is the gene labels
            row = reader.next()
            geneLabels = row[2:]        # skip the first two blank columns
            # Create a gene dictionary to hold the Ct values
            geneDict = {}
            for gene in geneLabels:
                geneDict[gene] = []

            # Read the sample and Ct data from each well for all genes
            geneDict["samples"] = []
            for well in range(wellN - 1):
                row = reader.next()
                if len(row) > 1:
                    geneDict["samples"].append(row[1])
                    Ct = [float(i) for i in row[2:]]
                    # Allocate these Ct values into geneDict
                    for i in range(len(Ct)):
                        gene = geneLabels[i]
                        geneDict[gene].append(Ct[i])
            
            # This has now read all the data required, unless we also need
            # the 'Quality Results' data
            break
    
    inputStream.close()
    return geneDict


# ------------------------------------------------------------------
def writeCSV(geneDict):
    """Output to STDOUT a .csv format of the processed data."""
    
    print "%s," % "Well",
    print "%s," % "Type",
    print "%s," % "Sample",
    print "%s," % "Gene",
    print "%s," % "Ct",
    print "%s," % "Quantity",
    print "%s\n" % "Exclusion",

    # The well numbers in the output stream should be continuous,
    # rather than repeat the values in the well array.  The items
    # variable is just an iteration counter.
    items = 0
    
    samples = geneDict["samples"]
    #samples = [s.replace("hr ", "hrs") for s in samples]
    #samples = [s.replace("hrs", " hrs") for s in samples]
    #samples = [s.replace("Days", "days") for s in samples]
    
    geneLabels = geneDict.keys()
    geneLabels.sort()
    geneLabels.remove("samples")
    for gene in geneLabels:
        Ct = geneDict[gene]
        gene = gene.replace(" (house keeping)", "_hk")
        
        for well in range(len(Ct)):
            items += 1
            #print "%04d, " % (well + 1),
            print "%04d," % items,
            print "%s," % "UNKN",
            print "%s," % samples[well],
            print "%s," % gene,
            if Ct[well] == 999:
                print ", , %d\n" % (1),
            else:
                print "%5.2f, , \n" % Ct[well],
    
    return None



# ------------------------------------------------------------------
def findAllIndices(thisList, value):
 	return [i for i in xrange(len(thisList)) if thisList[i] == value]

# ------------------------------------------------------------------
def checkSampleRange(geneDict, rangeThreshold = 0.5):
    """Check the sample replicates for excess variations.
    
    The default range criterion is 0.5.  If the range of values in a
    sample set is greater than 0.5, the function will exclude the min
    or max value from the set (whichever reduces the range more).
    This is done by assignment of 999 to exclude the value.  This
    process continues until the range < criterion or there is only 1
    value left.
    """
    
    samples = geneDict["samples"]
    uniqueSamples = dict.fromkeys(samples).keys()
    
    geneLabels = geneDict.keys()
    geneLabels.sort()
    geneLabels.remove("samples")
    
    for gene in geneLabels:
        Ct = geneDict[gene]
        # Find all the matching sample values
        for uniqueSample in uniqueSamples:
            CtSamplesIndex = findAllIndices(samples, uniqueSample)
            # Extract the sample Ct values
            CtSamples = [Ct[i] for i in CtSamplesIndex]
            # Copy the data and sort it for analysis, keeping the
            # original to facilitate substitutions for outliers.
            CtSorted = CtSamples[:]
            CtSorted.sort()
            # Exclude any 999 values
            while CtSorted.count(999.0):
                 CtSorted.remove(999.0)
            # Test the range against rangeThreshold
            if len(CtSorted) > 1:
                CtRange = max(CtSorted) - min(CtSorted)
            while len(CtSorted) > 2 and CtRange > rangeThreshold:
                
                # Calculate successive differences.  All differences
                # must be >=0 as the array is sorted.
                CtPaired = zip(CtSorted[:-1], CtSorted[1:])
                CtDif = [i[1] - i[0] for i in CtPaired]
                
                # Remove the min or max value from CtSorted.
                if CtDif[0] > CtDif[-1]:
                    CtResetValue = CtSorted[0]
                    del CtSorted[0]
                else:
                    CtResetValue = CtSorted[-1]
                    del CtSorted[-1]
                
                # Reset the extreme Ct value to 999.
                CtResetIndices = findAllIndices(CtSamples, CtResetValue)
                if len(CtResetIndices) == 0:
                    raise ValueError, "No CtReset values"
                for i in CtResetIndices:
                    # CtSamples and CtSamplesIndex are the same size,
                    # so we can use the index values from
                    # CtResetIndices to find the original value
                    # indices in geneDict.
                    CtIndex = CtSamplesIndex[i]
                    geneDict[gene][CtIndex] = 999.0
                
                if len(CtSorted) > 1:
                    CtRange = max(CtSorted) - min(CtSorted)
            
    return geneDict

# ------------------------------------------------------------------
def excludeSamples(geneDict, samplesExclude = ["Water"]):
    """Set all the Ct values to 999 to exclude a sample."""
    
    samples = geneDict["samples"]
    samplesUnique = dict.fromkeys(samples).keys()
    
    for sample in samplesExclude:
        if not sample in samplesUnique:
            msg = "The '%s' sample is not in this dataset" % (sample,)
            raise ValueError, msg
    
    geneLabels = geneDict.keys()
    geneLabels.sort()
    geneLabels.remove("samples")
    
    for sampleExclude in samplesExclude:
        # Find all the matching sample values
        CtSamplesIndex = findAllIndices(samples, sampleExclude)
        for gene in geneLabels:
            # Reset the Ct values to 999
            for i in CtSamplesIndex:
                geneDict[gene][i] = 999.0
    
    return geneDict



# ------------------------------------------------------------------
# Main

def main():
    # Parse command line options
    p = OptionParser()
    p.add_option("-i", "--inputfile", dest="inputfile",
                 action="store", type="string",
                 help="read input \"FILE.csv\"", metavar="FILE.csv")
    p.add_option("-o", "--outputfile", dest="outputfile",
                 action="store", type="string",
                 help="write output \"FILE.csv\"", metavar="FILE.csv")
    p.add_option("-r", "--sampleRange", dest="sampleRange", default=0.5,
                 action="store", type="float",
                 help="screen and exclude samples with range > SAMPLE_RANGE", metavar="SAMPLE_RANGE")
    p.add_option("-s", "--sampleExclude", dest="sampleExclude", default="",
                 action="store", type="string",
                 help="exclude sample names: \"SampleNameA SampleNameB\"", metavar="SAMPLE_NAME(S)")
    (options, args) = p.parse_args()
    #print options
    #sys.exit()
    
    # read input file
    if options.inputfile:
        geneDict = readCSV(options.inputfile)
    else:
        geneDict = readCSV()
    
    # Exclude some samples
    samples = ["Water"]
    for s in list(options.sampleExclude):
        samples.append(s)
    geneDict = excludeSamples(geneDict, samples)
    
    # Check the sample data range
    geneDict = checkSampleRange(geneDict, options.sampleRange)
    
    # write output file
    if options.outputfile:
        writeCSV(geneDict, options.outputfile)
    else:
        writeCSV(geneDict)
    # Use pyExcelerate to write out an .xls file?

if __name__ == "__main__":
    main()
    sys.exit()
