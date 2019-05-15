class Sample:

    def __init__(self,sequence,insertPos,intronPos,circular):
        self.circular = circular
        self.plasmid = self.circularPlasmid(sequence,insertPos) if circular else [sequence[:insertPos[0]],sequence[insertPos[1]:]]
        self.gDNA1 = OrientedSample(sequence,insertPos,intronPos,True,True)
        self.gDNA2 = OrientedSample(sequence,insertPos,intronPos,True,False)
        self.cDNA1 = OrientedSample(sequence,insertPos,intronPos,False,True)
        self.cDNA2 = OrientedSample(sequence,insertPos,intronPos,False,False)

    def circularPlasmid(self,sequence,insertPos):
        plasmid = list(sequence[insertPos[1]:])
        plasmid.extend([sequence[:insertPos[0]]])
        return plasmid

class OrientedSample:

    def __init__(self,sequence,insertPos,intronPos,hasIntron,isForward):
        insertArea = sequence[insertPos[0]-1:insertPos[1]]
        newInsertPos = insertPos
        comp = {'A':'T','C':'G','T':'A','G':'C'}
        if not hasIntron:
            beginEx = insertPos[0]-1
            insertArea = ""
            for intron in intronPos:
                insertArea += sequence[beginEx:intron[0]-1]
                beginEx = intron[1]
            insertArea += sequence[beginEx:insertPos[1]]
            newInsertPos = (insertPos[0],insertPos[0]+len(insertArea)-1)
        if not isForward:
            insertArea = ''.join([comp[x] for x in insertArea[::-1]])
        introns = intronPos
        if hasIntron and not isForward:
            introns = [(insertPos[1]+insertPos[0]-x[1],insertPos[1]+insertPos[0]-x[0]) for x in intronPos]
        self.introns = introns
        self.insertPos = newInsertPos
        self.sequence = sequence[:insertPos[0]-1] + insertArea + sequence[insertPos[1]:]
        self.hasIntron = hasIntron
        self.isForward = isForward