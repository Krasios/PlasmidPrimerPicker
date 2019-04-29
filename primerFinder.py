from model import Sample,OrientedSample

MIN_PRIMER_LENGTH =18
segmentTable = {}

def findGCClamps(seq):
    GCClampStarts = []
    for i in range(len(seq)-5):
        gc=sum([1 for x in seq[i:i+5] if x == 'G' or x =='C'])
        if (gc < 4 and gc > 0):
            GCClampStarts.append(i)
    return GCClampStarts

def findPossibleInsertPrimers(sample):
    mainSample = sample.gDNA1
    commonInsertRegions = [mainSample.sequence[mainSample.insertPos[0]:mainSample.insertPos[1]]]
    if (len(mainSample.introns) > 0):
        commonInsertRegions = [mainSample.sequence[mainSample.insertPos[0]:mainSample.introns[0][0]]]
        commonInsertRegions.append(mainSample.sequence[mainSample.introns[-1][1]:mainSample.insertPos[1]])
        commonInsertRegions.extend([mainSample.sequence[mainSample.introns[x][1]:mainSample.introns[x+1][0]] for x in range(len(mainSample.introns)-1)])

    ##filter out fragments too short to match a primer
    commonInsertRegions = [x for x in commonInsertRegions if len(x) > MIN_PRIMER_LENGTH]
    ##print(commonInsertRegions)
    possibleStarts = [(x,findGCClamps(x)) for x in commonInsertRegions]
    populateSegTable(possibleStarts)
    print(segmentTable)
    return []

##make segments
def generateComplement(segment):
    comp = {'A':'T','C':'G','T':'A','G':'C'}
    return ''.join([comp[i] for i in segment[::-1]])

def populateSegTable(segStarts):
    for starts in segStarts:
        segLen = len(starts[0])
        for j in starts[1]:
            if j >= 13:
                segmentTable[generateComplement(starts[0][j-13:j+5])] = calculateTM(starts[0][j-13:j+5])
                if j >= 14:
                    segmentTable[generateComplement(starts[0][j-14:j+5])] = calculateTM(starts[0][j-14:j+5])
                    if j >= 15:
                        segmentTable[generateComplement(starts[0][j-15:j+5])] = calculateTM(starts[0][j-15:j+5])
                        if j >= 16:
                            segmentTable[generateComplement(starts[0][j-16:j+5])] = calculateTM(starts[0][j-16:j+5])
                            if j >= 17:
                                segmentTable[generateComplement(starts[0][j-17:j+5])] = calculateTM(starts[0][j-17:j+5])

            if j < segLen - 13:
                segmentTable[starts[0][j:j+18]] = calculateTM(starts[0][j:j+18])
                if j < segLen - 14:
                    segmentTable[starts[0][j:j+19]] = calculateTM(starts[0][j:j+19])
                    if j < segLen - 15:
                        segmentTable[starts[0][j:j+20]] = calculateTM(starts[0][j:j+20])
                        if j < segLen - 16:
                            segmentTable[starts[0][j:j+21]] = calculateTM(starts[0][j:j+21])
                            if j < segLen - 17:
                                segmentTable[starts[0][j:j+22]] = calculateTM(starts[0][j:j+22])



def calculateTM(primer):
    gcCount = sum([1 for x in primer if x == 'G' or x =='C'])
    atCount = sum([1 for x in primer if x == 'A' or x =='T'])
    return 64.9+(41*(gcCount-16.4)/(gcCount+atCount))