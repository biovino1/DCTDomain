"""
Takes sequences from pfam_localpfam_nomax50.pair and takes their fasta sequences from
pfam_nomax50.fasta and puts them into a new file.

__author__ = "Ben Iovino"
__date__ = "8/15/2023"
"""

from Bio import SeqIO
from Bio.SeqIO import FastaIO

def main():


    # Get all sequences from pfam_localpfam_nomax50.pair
    seqs = []
    with open('pfam_data/pfam_localpfam_nomax50.pair', 'r', encoding='utf8') as pfile:
        for line in pfile:
            line = line.split()
            if line[0] not in seqs:
                seqs.append(line[0])
            if line[1] not in seqs:
                seqs.append(line[1])

    # Get their fasta seqs from pfam_nomax50.fasta
    with open('pfam_data/pfam_nomax50.fasta', 'r', encoding='utf8') as ffile:
        with open('pfam_data/pfam_localpfam_nomax50.fasta', 'w', encoding='utf8') as outfile:
            for record in SeqIO.parse(ffile, 'fasta'):
                if record.id in seqs:
                    FastaIO.FastaWriter(outfile, wrap=None).write_record(record)


if __name__ == '__main__':
    main()
