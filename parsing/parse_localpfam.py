"""
Takes sequences from pfam_localpfam_nomax50.pair and takes their fasta sequences from
pfam_nomax50.fasta and puts them into a new file.

__author__ = "Ben Iovino"
__date__ = "8/15/2023"
"""

from Bio import SeqIO
from Bio.SeqIO import FastaIO


def main():


    # Get all sequences from pfam_localpfam_nomax50.info
    seqs, pairs = [], []
    with open('pfam_data/pfam_localpfam_nomax50.info', 'r', encoding='utf8') as pfile:
        for line in pfile:
            line = line.split()
            if line[0] not in seqs:
                seqs.append(line[0])
            if line[1] not in seqs:
                seqs.append(line[1])
            pairs.append(f'{line[0]} {line[1]} hom')

    # Get pairs from pfam_localpfam_nomax50-dctsim.txt
    sim_pairs = []
    with open('pfam_data/pfam_localpfam_nomax50-dctsim.txt', 'r', encoding='utf8') as pfile:
        for line in pfile:
            sim_pairs.append((line.split()[0], line.split()[1]))

    # Write pairs to file
    with open('pfam_data/pfam_localpfam_nomax50.pair', 'w', encoding='utf8') as pfile:
        pfile.write('#prot1 prot2 label\n')
        for pair in pairs:
            if (pair.split()[0], pair.split()[1]) in sim_pairs:
                pfile.write(f'{pair}\n')

    # Get their fasta seqs from pfam_nomax50.fasta
    with open('pfam_data/pfam_nomax50.fasta', 'r', encoding='utf8') as ffile:
        with open('pfam_data/pfam_localpfam_nomax50.fasta', 'w', encoding='utf8') as outfile:
            for record in SeqIO.parse(ffile, 'fasta'):
                if record.id in seqs:
                    FastaIO.FastaWriter(outfile, wrap=None).write_record(record)


if __name__ == '__main__':
    main()
