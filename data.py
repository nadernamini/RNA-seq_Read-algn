from Bio import SeqIO
import numpy as np
import csv
from shared import *
from project import Aligner
from time import time


def read_reads():
    """
    Reads all the RNA reads from reads.fa and returns a list of sequences as strings.

    :return: a list of RNA sequence reads as strings
    """
    max_id_len = 60
    reads = np.zeros((1575, 2), dtype='S{:d}'.format(max_id_len))
    i = 0
    for seq_record in SeqIO.parse("reads.fa", "fasta"):
        reads[i, 0] = seq_record.id
        reads[i, 1] = str(seq_record.seq)
        i += 1
    return reads


def read_genome():
    """
    Reads the genome sequence from genome.fa and return the genome sequence as a string.

    :return: a string of the genome sequence
    """
    genome = None
    for seq_record in SeqIO.parse("genome.fa", "fasta"):
        genome = str(seq_record.seq)
    return genome


def read_known_genes():
    """
    Reads the un/known genes, isoforms, and exons from genes.tab and constructs objects for each
    and return the list of constructed genes.

    :return: a list of known Gene objects
    """
    genes, isoforms, exons = {}, {}, {}
    with open("genes.tab") as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            name = line[0].split('_')[-1]
            if name == 'gene':
                genes[line[1]] = line[2].split(';')
            elif name == 'isoform':
                isoforms[line[1]] = line[2].split(';')
            elif name == 'exon':
                exons[line[1]] = Exon(line[1], int(line[2]), int(line[3]))
    # Create the Isoform objects
    for k in isoforms:
        isoforms[k] = Isoform(k, [exons[key] for key in isoforms[k]])

    # Create the Genes objects
    for k in genes:
        genes[k] = Gene(k, [isoforms[key] for key in genes[k]])
    return genes


if __name__ == "__main__":
    reads = read_reads()
    genome_sequence = read_genome()
    known_genes = read_known_genes()

    t = -time()
    aligner = Aligner(genome_sequence, known_genes)
    t += time()

    print("time to run Aligner.__init__: " + str(t))
