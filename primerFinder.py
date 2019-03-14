from model import Sample,OrientedSample

MIN_PRIMER_LENGTH =18
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
    print(commonInsertRegions)
    possibleStarts = [(x,findGCClamps(x)) for x in commonInsertRegions]
    print(possibleStarts)
    return []