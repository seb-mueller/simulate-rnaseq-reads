# Sebastian Y. Mueller 2017, sm934@cam.ac.uk
# /applications/anaconda/anaconda3/bin/python3 simDropSeqSeqs.py --T Drosophila_melanogaster.BDGP6.cdna.all.fa Mus_musculus.GRCm38.cdna.all.fa Homo_sapiens.GRCh38.cdna.all.fa --seqlength 50 --nrcells 40 --nrseqs 10000 --Outputname melcells > melcells.cells

# O(1) instead O(N) complexity to draw random entry from a dict
from randomdict import RandomDict
from Bio import SeqIO
# for generating random sequences
from random import randint
import random
from Bio.SeqRecord import SeqRecord
import argparse
parser = argparse.ArgumentParser()
# args = parser.parse_args(['--T', 'test', 'sim.cells'])
parser.add_argument("--T", type=argparse.FileType('r'), nargs='+',
                    help="Transcriptome(s) [fasta].\
                    If multiple files, seperated by space.")
parser.add_argument("--Outputname", type=str,
                    help="Base name of output files [fastq].\
                    Two files *_R1.fastq and *_R2.fastq will be generated.")
parser.add_argument("--nrcells", type=int, nargs='?',
                    help="Number of simulated cells", default=50)
parser.add_argument("--nrseqs", type=int, nargs='?',
                    help="Number of simulated sequences per cell",
                    default=5000)
parser.add_argument("--seqlength",
                    nargs='?',
                    help="Length of generated sequences (in base pairs)",
                    default=100,
                    type=int)
args = parser.parse_args()

myLength = args.seqlength
nrCells = args.nrcells
nrSeqs = args.nrseqs
nrSpecies = len(args.T)

print(str(myLength) + "|" + str(nrCells) + "|" + str(nrSeqs))


def randomseq(length):
    return ''.join(random.choice("ATCG") for i in range(length))


# converting to fastq: we need to add made up qualities to the SeqRecord
# object in using the .letter_annotations attribute (see chap 4 biopy tutorial)

# record_dict = SeqIO.index("mmsubset.fa", "fasta")
# record_dict = SeqIO.index("Mus_musculus.GRCm38.68.cdna.all.fa", "fasta")
# record_dict = SeqIO.index(args.T[0].name, "fasta")
# coecing into a RandomDict for O(1) random access
# record_rand = RandomDict(record_dict)

count = 1
with open(args.Outputname + "_R1.fastq", "w") as handleR1, open(args.Outputname + "_R2.fastq", "w") as handleR2:
    for myfile in args.T:
        print("created dict for " + myfile.name)
        record_dict = SeqIO.index(myfile.name, "fasta")
        record_rand = RandomDict(record_dict)
        for cell in range(0, nrCells):
            cellSeq = randomseq(12)
            print(cellSeq)
            for gene in range(0, nrSeqs):
                UMIseq = randomseq(8)
                myItem = record_rand.random_value()
                limit = len(myItem)
                myItem += "A" * 102
                start = randint(0, limit)
                # print(myItem.id + " length: " + str(limit))
                end = start + myLength
                mySub = myItem[start:end]
                seqid = myItem.id + ":" + cellSeq + ":" + UMIseq + ":" + str(count)
                mySub.id = seqid + "/2"
                mySub.letter_annotations["phred_quality"] = [24] * len(mySub)
                tmp = SeqIO.write(mySub, handleR2, "fastq")
                mySub2 = SeqRecord(cellSeq + UMIseq,
                                   id=seqid + "/1",
                                   description=mySub.description)
                mySub2.letter_annotations["phred_quality"] = [24] * len(mySub2)
                tmp = SeqIO.write(mySub2, handleR1, "fastq")
                count += 1
