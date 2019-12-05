

# Sebastian Y. Mueller 2017, sm934@cam.ac.uk
# /applications/anaconda/anaconda3/bin/python3 simDropSeqSeqs.py --T Drosophila_melanogaster.BDGP6.cdna.all.fa Mus_musculus.GRCm38.cdna.all.fa Homo_sapiens.GRCh38.cdna.all.fa --seqlength 50 --nrcells 40 --nrseqs 10000 --Outputname melcells > melcells.cells

# O(1) instead O(N) complexity to draw random entry from a dict
from randomdict import RandomDict
from Bio import SeqIO
from Bio.Seq import Seq
# for generating random sequences
from random import randint
import random
from Bio.SeqRecord import SeqRecord
import argparse
import Levenshtein
import pandas as pd
from copy import deepcopy
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
parser.add_argument("--cell_bc_length",
                    help="Length of cell barcode",
                    default=8,
                    type=int)
parser.add_argument("--umi_bc_length",
                    help="Length of umi barcode",
                    default=12,
                    type=int)
parser.add_argument("--whitelist_out",
                    help="Cell barcode list output file",
                    default='whitelist_barcodes.csv',
                    type=str,
                    dest='whitelist')
parser.add_argument("--input-matrix",
                    help="Input matrix of data to produce",
                    default=None,
                    type=str,
                    dest='input_matrix')
parser.add_argument("--amplification-biais",
                    help="Duplication of reads biais",
                    default=0,
                    type=float,
                    dest='amplification')

# parser.add_argument("--reads",
#                     required=True,
#                     help="Number of reads to produce",
#                     default=None,
#                     type=int,
#                     dest='reads')


args = parser.parse_args()


myLength = args.seqlength
nrCells = args.nrcells
nrSeqs = args.nrseqs
nrSpecies = len(args.T)

#print(str(myLength) + "|" + str(nrCells) + "|" + str(nrSeqs))


def randomseq(length, letters):
    """
    Returns a sequence of a certain length with a composition of letters given as input.

    Args:
        length (int): Length of the sequence
        letters (string): A string containing all the possible letters

    Returns:
        Sequence (string): A sequence in uppercase

    """
    sequence = ''.join(random.choice(letters) for i in range(length))
    return(sequence.upper())


def mutate(sequence, mutation_rate=0.0066, letters=set('NATCG')):
    """
    Mutate a sequence with a given mutation rate and letters to pick from.

    Args:
        sequence (string): Sequence to mutate
        mutation_rate (float): Mutation rate
        letters (string): Letter to pick from to mutate

    Returns:
        mutated_sequence (string): Mutated sequence
    """
    result = []
    mutations = []
    for base in sequence:
        if random.random() <= mutation_rate:
            new_letters = deepcopy(letters)
            new_letters.remove(base)
            new_base = random.sample(new_letters, 1)[0] # negatives are OK
            result.append(new_base)
        else:
            result.append(base)
    mutated_sequence = "".join(result)
    return(mutated_sequence)

def generate_umi_set(length, n, mutation_rate):
    umis = set()
    for i in range(n):
        seq = mutate(randomseq(length, 'ATGC'), mutation_rate=mutation_rate)
        umis.add(seq)
    return(umis)


# converting to fastq: we need to add made up qualities to the SeqRecord
# object in using the .letter_annotations attribute (see chap 4 biopy tutorial)

# record_dict = SeqIO.index("mmsubset.fa", "fasta")
# record_dict = SeqIO.index("Mus_musculus.GRCm38.68.cdna.all.fa", "fasta")
# record_dict = SeqIO.index(args.T[0].name, "fasta")
# coecing into a RandomDict for O(1) random access
# record_rand = RandomDict(record_dict)
# if args.input_matrix:
#     input_data = pd.read_csv(args.input_matrix, index_col=0)
# print(input_data)
# exit()


# Create a cell pool to sample from.

def generate_cell_pool(length, letters, min_lev_dist):
    cell_pool = set([randomseq(args.cell_bc_length, letters)])

    while len(cell_pool) != nrCells:
        new_cell = randomseq(args.cell_bc_length, letters)
        good = True
        for cell in cell_pool:
            if Levenshtein.hamming(cell, new_cell) < min_lev_dist:
                good = False
        if good:
            cell_pool.add(new_cell)

    return(cell_pool)

def generate_mutants(cell_pool, base_mutation_rate, number_of_mutants):
    mutants = set()
    print('Make mutants for cell pool')
    for cell in cell_pool:
        for i in range(number_of_mutants):
            mutants.add(mutate(cell, base_mutation_rate))
    cell_pool.update(mutants)

    return(cell_pool)

cell_pool = generate_cell_pool(length=args.cell_bc_length, letters='ATGC', min_lev_dist=5)

#Write the whitelist to disk
with open(args.whitelist, 'w') as whitelist_file:
    for cell_barcode in cell_pool:
        whitelist_file.write('{}\n'.format(cell_barcode))

#Generate some mutants
mutants = generate_mutants(cell_pool=cell_pool, base_mutation_rate=0.001, number_of_mutants=10)



print('Writing sim data to file.')
count = 1
with open(args.Outputname + "_R1.fastq", "w") as handleR1, open(args.Outputname + "_R2.fastq", "w") as handleR2:
    for myfile in args.T:
        record_dict = SeqIO.index(myfile.name, "fasta")
        record_rand = RandomDict(record_dict)
        for cell in range(0, nrCells):
            cellSeq = random.sample(cell_pool, 1)[0]
            cell_pool.discard(cellSeq)
            umi_set = generate_umi_set(args.umi_bc_length, 100, 0.2)
            for record in range(0, nrSeqs):
                UMIseq = random.sample(umi_set, 1)[0]
                myItem = record_rand.random_value()
                if random.random() < 0.01:
                    myItem.seq = Seq(mutate(sequence=myItem.seq, mutation_rate=0.001))
                myItem += "A" * 150
                start = 0
                #print(myItem.id + " length: " + str(limit))
                end = start + myLength
                mySub = myItem[start:end]
                seqid = myItem.id + ":" + cellSeq + ":" + UMIseq + ":" + str(count)
                mySub.id = seqid + "/2"
                mySub.letter_annotations["phred_quality"] = [24] * len(mySub)
                tmp = SeqIO.write(mySub, handleR2, "fastq")
                    
                mySub2 = SeqRecord(Seq(cellSeq + UMIseq),
                                   id=seqid + "/1",
                                   description=mySub.description)
                mySub2.letter_annotations["phred_quality"] = [24] * len(mySub2)

                tmp = SeqIO.write(mySub2, handleR1, "fastq")
                count += 1
