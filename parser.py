from Bio import SeqIO
from model import Sample
from primerFinder import *

arch = "pcr4-gDNA1.gb"
record = SeqIO.parse(arch, "genbank")
rec = next(record)
##print(rec.annotations['topology'])   
intronPositions = []                    # there is only one record
for f in rec.features:
    if f.type.lower() == 'insert' or 'insert' in [l.lower() for l in f.qualifiers['label']]:
        ##print("Insert sequence begins at {} and ends at {}".format(f.location.start,f.location.end))
        insertPosition = (f.location.start.position,f.location.end.position)
    if f.type.lower() == 'intron' or 'intron' in [l.lower() for l in f.qualifiers['label']]:
        ##print("Intron sequence begins at {} and ends at {}".format(f.location.start,f.location.end))
        intronPositions.append((f.location.start.position,f.location.end.position))
sample = Sample(rec.seq,insertPosition,intronPositions,rec.annotations['topology']=="circular")
insertPrimers = findPossiblePrimers(sample)
#print(insertPrimers)
primers = findCompatiblePlasmidPrimerPair(sample,insertPrimers[0],False)
if primers == None:
    primers = findCompatiblePlasmidPrimerPair(sample,insertPrimers[1],True)
#website uses primers object
print(primers)

