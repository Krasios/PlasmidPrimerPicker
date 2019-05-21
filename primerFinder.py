from model import Sample,OrientedSample

MIN_PRIMER_LENGTH =18


def findGCClamps(seq):
    GCClampStarts = []
    for i in range(len(seq)-5):
        gc=sum([1 for x in seq[i:i+5] if x == 'G' or x =='C'])
        if (gc < 4 and gc > 0):
            GCClampStarts.append(i)
    return GCClampStarts

def findPossiblePrimers(sample):
    mainSample = sample.gDNA1
    commonInsertRegions = [mainSample.sequence[mainSample.insertPos[0]:mainSample.insertPos[1]]]
    if (len(mainSample.introns) > 0):
        commonInsertRegions = [mainSample.sequence[mainSample.insertPos[0]:mainSample.introns[0][0]]]
        commonInsertRegions.append(mainSample.sequence[mainSample.introns[-1][1]:mainSample.insertPos[1]])
        commonInsertRegions.extend([mainSample.sequence[mainSample.introns[x][1]:mainSample.introns[x+1][0]] for x in range(len(mainSample.introns)-1)])

    ##filter out fragments too short to match a primer
    commonInsertRegions = [x for x in commonInsertRegions if len(x) > MIN_PRIMER_LENGTH]
    possibleStarts = [(x,findGCClamps(x)) for x in commonInsertRegions]
    segmentList = []
    
    populateTMList(possibleStarts,segmentList,segmentList)   
    segmentList = [x for x in segmentList if GCInBounds(x[0])]
    plasmidPairs = findPlasmidPrimerPairs(sample)
    plasmidPairs = [x for x in plasmidPairs if GCInBounds(x[0][0]) and GCInBounds(x[1][0])] 
    return [[x[0],y,x[1]] for x in plasmidPairs for y in segmentList if abs(x[0][1]-y[1])<5 and abs(x[1][1]-y[1])<5]

def findPlasmidPrimerPairs(sample):
    possiblePlasmidStarts = [(x,findGCClamps(x)) for x in sample.plasmid]
    if (sample.circular):
        plasmidForwardList = []
        plasmidReverseList = []
        populateTMList(possiblePlasmidStarts,plasmidForwardList,plasmidReverseList)
        return [[x,y] for x in plasmidForwardList for y in plasmidReverseList if abs(x[1]-y[1])<10]
    plasmidForwardList = [[],[]]
    plasmidReverseList = [[],[]]
    populateTMList([possiblePlasmidStarts[0]],plasmidForwardList[0],plasmidReverseList[0])
    populateTMList([possiblePlasmidStarts[1]],plasmidForwardList[1],plasmidReverseList[1])
    retList = [[x,y] for x in plasmidForwardList[0] for y in plasmidReverseList[1] if abs(x[1]-y[1])<10]
    retList.extend([[x,y] for x in plasmidForwardList[1] for y in plasmidReverseList[0] if abs(x[1]-y[1])<10])
    return retList

##make segments
def generateComplement(segment):
    comp = {'A':'T','C':'G','T':'A','G':'C'}
    return ''.join([comp[i] for i in segment[::-1]])

def populateTMList(segStarts,forwardList,reverseList):
    for starts in segStarts:
        segLen = len(starts[0])
        for j in starts[1]:
            if j >= 13:
                forwardList.append((generateComplement(starts[0][j-13:j+5]),calculateTM(starts[0][j-13:j+5])))
                if j >= 14:
                    forwardList.append((generateComplement(starts[0][j-14:j+5]),calculateTM(starts[0][j-14:j+5])))
                    if j >= 15:
                        forwardList.append((generateComplement(starts[0][j-15:j+5]),calculateTM(starts[0][j-15:j+5])))
                        if j >= 16:
                            forwardList.append((generateComplement(starts[0][j-16:j+5]),calculateTM(starts[0][j-16:j+5])))
                            if j >= 17:
                                forwardList.append((generateComplement(starts[0][j-17:j+5]),calculateTM(starts[0][j-17:j+5])))

            if j < segLen - 13:
                reverseList.append((starts[0][j:j+18],calculateTM(starts[0][j:j+18])))
                if j < segLen - 14:
                    reverseList.append((starts[0][j:j+19],calculateTM(starts[0][j:j+19])))
                    if j < segLen - 15:
                        reverseList.append((starts[0][j:j+20],calculateTM(starts[0][j:j+20])))
                        if j < segLen - 16:
                            reverseList.append((starts[0][j:j+21],calculateTM(starts[0][j:j+21])))
                            if j < segLen - 17:
                                reverseList.append((starts[0][j:j+22],calculateTM(starts[0][j:j+22])))
##todo:check for selfdimers
##may switch out for other formulas to better match neb
def calculateTM(primer):
    gcCount = sum([1 for x in primer if x == 'G' or x =='C'])
    atCount = sum([1 for x in primer if x == 'A' or x =='T'])
    return 64.9+(41*(gcCount-16.4)/(gcCount+atCount))
def GCInBounds(primer):
    gcCount = sum([1 for x in primer if x == 'G' or x =='C'])
    gcPercent = gcCount/len(primer)*100
    return gcPercent>40 and gcPercent<60