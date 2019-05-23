from model import Sample,OrientedSample
import re
MIN_PRIMER_LENGTH =18
#approximated using ladders
MIN_DIFF=[(100,50),(1000,100),(10000,1000)]


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
    segmentList = [set(),set()]
    
    populateTMList(possibleStarts,segmentList[0],segmentList[1])   
    segmentList = [[x for x in y if GCInBounds(x[0])] for y in segmentList]
    return segmentList
    
def findCompatiblePlasmidPrimerPair(sample,primers,reverse):
    possiblePlasmidStarts = [(x,findGCClamps(x)) for x in sample.plasmid]
    if (sample.circular):
        plasmidLen = len(sample.plasmid)
        plasmidForwardSet = set()
        plasmidReverseSet = set()
        populateTMList(possiblePlasmidStarts,plasmidForwardSet,plasmidReverseSet)
        for primerF in plasmidForwardSet:
            for primerR in plasmidReverseSet:
                for primerI in primers:
                    if abs(primerI[1]-primerF[1])<5 and abs(primerI[1]-primerR[1])<5:
                        amplicons = getAmpliconSizes(sample,primerI[0],primerF[0],primerR[0],reverse)
                        if (distinguishable(amplicons)):
                            return (primerF,(generateComplement(primerR[0]),primerR[1]),primerI if not reverse else (generateComplement(primerI[0]),primerI[1]),amplicons)
    return None
    #plasmidForwardList = [[],[]]
    #plasmidReverseList = [[],[]]
    #populateTMList([possiblePlasmidStarts[0]],plasmidForwardList[0],plasmidReverseList[0])
    #populateTMList([possiblePlasmidStarts[1]],plasmidForwardList[1],plasmidReverseList[1])
    #retList = [[x,y] for x in plasmidForwardList[0] for y in plasmidReverseList[1] if abs(x[1]-y[1])<10]
    #retList.extend([[x,y] for x in plasmidForwardList[1] for y in plasmidReverseList[0] if abs(x[1]-y[1])<10])
    #return retList

def getAmpliconSizes(sample,insertPrimer,plasmidPrimer,plasmidPrimer2,reverse):
    insertPrimerComp = generateComplement(insertPrimer) if reverse else insertPrimer
    insertPrimerComp = '(?=('+''.join(insertPrimerComp)+'))'
    insertPrimerAligned = insertPrimer if reverse else generateComplement(insertPrimer)
    insertPrimerAligned = '(?=('+''.join(insertPrimerAligned)+'))'
    plasmidPrimer = '(?=('+''.join(plasmidPrimer)+'))'
    plasmidPrimer2 = '(?=('+''.join(plasmidPrimer2)+'))'
    insertMatches0 = [i.start(1) for i in re.finditer(insertPrimerComp,str(sample.gDNA1.sequence))]
    insertMatches1 = [i.start(1) for i in re.finditer(insertPrimerAligned,str(sample.gDNA2.sequence))]
    insertMatches2 = [i.start(1) for i in re.finditer(insertPrimerComp,str(sample.cDNA1.sequence))]
    insertMatches3 = [i.start(1) for i in re.finditer(insertPrimerAligned,str(sample.cDNA2.sequence))]
    plasmidMatches0 = [i.start(1) for i in re.finditer(plasmidPrimer,str(sample.gDNA1.sequence))]
    plasmidMatches1 = [i.start(1) for i in re.finditer(plasmidPrimer,str(sample.gDNA2.sequence))]
    plasmidMatches2 = [i.start(1) for i in re.finditer(plasmidPrimer,str(sample.cDNA1.sequence))]
    plasmidMatches3 = [i.start(1) for i in re.finditer(plasmidPrimer,str(sample.cDNA2.sequence))]
    plasmidMatches20 = [i.start(1) for i in re.finditer(plasmidPrimer2,str(sample.gDNA1.sequence))]
    plasmidMatches21 = [i.start(1) for i in re.finditer(plasmidPrimer2,str(sample.gDNA2.sequence))]
    plasmidMatches22 = [i.start(1) for i in re.finditer(plasmidPrimer2,str(sample.cDNA1.sequence))]
    plasmidMatches23 = [i.start(1) for i in re.finditer(plasmidPrimer2,str(sample.cDNA2.sequence))]
    amplicons = []
    amplicons.append(plasmidMatches0[-1]-insertMatches0[0] if reverse else insertMatches0[-1]-plasmidMatches20[0])
    amplicons.append(insertMatches1[-1]-plasmidMatches21[0] if reverse else plasmidMatches1[-1]-insertMatches1[0])
    amplicons.append(plasmidMatches2[-1]-insertMatches2[0] if reverse else insertMatches2[-1]-plasmidMatches22[0])
    amplicons.append(insertMatches3[-1]-plasmidMatches23[0] if reverse else plasmidMatches3[-1]-insertMatches3[0])
    return amplicons

def distinguishable(amplicons):
    minFrag = min(amplicons)
    for x in MIN_DIFF:
        if minFrag <= x[0] and abs(amplicons[0]-amplicons[2])>x[1] and abs(amplicons[1]-amplicons[3])>=x[1]:
            return True
    return False

##make segments
def generateComplement(segment):
    comp = {'A':'T','C':'G','T':'A','G':'C'}
    return ''.join([comp[i] for i in segment[::-1]])

def populateTMList(segStarts,forwardList,reverseList):
    for starts in segStarts:
        segLen = len(starts[0])
        for j in starts[1]:
            if j >= 13:
                forwardList.add((tuple(starts[0][j-13:j+5]),calculateTM(starts[0][j-13:j+5])))
                if j >= 14:
                    forwardList.add((tuple(starts[0][j-14:j+5]),calculateTM(starts[0][j-14:j+5])))
                    if j >= 15:
                        forwardList.add((tuple(starts[0][j-15:j+5]),calculateTM(starts[0][j-15:j+5])))
                        if j >= 16:
                            forwardList.add((tuple(starts[0][j-16:j+5]),calculateTM(starts[0][j-16:j+5])))
                            if j >= 17:
                                forwardList.add((tuple(starts[0][j-17:j+5]),calculateTM(starts[0][j-17:j+5])))

            if j < segLen - 13:
                reverseList.add((tuple(starts[0][j:j+18]),calculateTM(starts[0][j:j+18])))
                if j < segLen - 14:
                    reverseList.add((tuple(starts[0][j:j+19]),calculateTM(starts[0][j:j+19])))
                    if j < segLen - 15:
                        reverseList.add((tuple(starts[0][j:j+20]),calculateTM(starts[0][j:j+20])))
                        if j < segLen - 16:
                            reverseList.add((tuple(starts[0][j:j+21]),calculateTM(starts[0][j:j+21])))
                            if j < segLen - 17:
                                reverseList.add((tuple(starts[0][j:j+22]),calculateTM(starts[0][j:j+22])))

##may switch out for other formulas to better match neb
def calculateTM(primer):
    gcCount = sum([1 for x in primer if x == 'G' or x =='C'])
    atCount = sum([1 for x in primer if x == 'A' or x =='T'])
    return 64.9+(41*(gcCount-16.4)/(gcCount+atCount))

def GCInBounds(primer):
    gcCount = sum([1 for x in primer if x == 'G' or x =='C'])
    gcPercent = gcCount/len(primer)*100
    return gcPercent>40 and gcPercent<60