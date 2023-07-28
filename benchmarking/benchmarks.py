"""================================================================================================
This script takes all pairs from pfam_max50.pair and performs homology searchs on each one, with
the first sequence acting as the 'sequence database' and the second sequence acting as the query
sequence. E-value and bit scores are saved to a file.

Ben Iovino  07/27/23   DCTDomain
================================================================================================"""

import os
import pickle
import subprocess
from Bio import SeqIO
from sklearn import metrics


def get_pairs(file: str) -> list:
    """=============================================================================================
    This function takes a file of pairs and returns a list of tuples.

    :param file: file of pairs
    :return list: list of tuples
    ============================================================================================="""

    pairs = []
    with open(file, 'r', encoding='utf8') as f:
        for line in f:
            pairs.append(tuple(line.strip().split()))

    return pairs[1:]  # Skip header


def get_seqs(file: str) -> dict:
    """=============================================================================================
    This function takes a fasta file and returns a dictionary of sequences and their IDs.

    :param file: fasta file
    :return dict: dictionary of sequences
    ============================================================================================="""

    # Read each line and add to dictionary
    seqs = {}
    with open(file, 'r', encoding='utf8') as f:
        for seq in SeqIO.parse(f, 'fasta'):
            seqs[seq.id] = str(seq.seq)

    return seqs


def phmmer_search(pairs: list, seqs: dict):
    """=============================================================================================
    This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'sequence database' and 
    the second sequence is used as the query sequence in a phmmer search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    ============================================================================================="""

    results = {}
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open('benchmarking/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open('benchmarking/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from phmmer search
        # 14th line of stdout, will be blank if no hits detected
        result = subprocess.getoutput('phmmer benchmarking/db_seq.fa benchmarking/query_seq.fa')
        result_line = result.split('\n')[14].split()

        # If there is a hit, check if the hit is the same as the query sequence
        if result_line == [] or result_line[-1] != pair[1]:
            result_line = 0
        else:
            result_line = [result_line[0], result_line[1], result_line[3], result_line[4]]
        results[(pair[0], pair[1])] = result_line

    # Save so we can load for later parsing, much faster than running search each time
    with open('benchmarking/phmmer_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def phmmer_results(pairs: list):
    """=============================================================================================
    This function loads the results from phmmer_search() and compares the results to the ground
    truth to determine the accuracy of the search.

    :param pairs: protein pairs and their label (hom/nonhom)
    ============================================================================================="""

    with open('benchmarking/phmmer_results.pkl', 'rb') as f:
        results = pickle.load(f)

    print(results)


def blast_search(pairs: list, seqs: dict):
    """=============================================================================================
    This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'blast database' and 
    the second sequence is used as the query sequence in a blastp search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    ============================================================================================="""

    results, db_seq = {}, ''
    for pair in pairs:

        # If first seq is different, make new blast database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            os.system('makeblastdb -in benchmarking/db_seq.fa '
                '-dbtype prot -parse_seqids -out benchmarking/blastdb/db_seq')
            with open('benchmarking/db_seq.fa', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open('benchmarking/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from blastp search
        result = subprocess.getoutput('blastp -query benchmarking/query_seq.fa '
                                      '-db benchmarking/blastdb/db_seq')

        # 29th line of stdout, will be blank if no hits detected
        result_line = result.split('\n')[29].split()[-2:]
        if result_line == []:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    with open('benchmarking/blast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def main():

    # Get all pairs and seqs from max50 files
    max50_pairs = 'embedding/pfam_max50.pair'
    max50_seqs = 'embedding/pfam_max50.fasta'
    pairs = get_pairs(max50_pairs)
    seqs = get_seqs(max50_seqs)

    # Perform phmmer search on each pair
    phmmer_search(pairs, seqs)
    phmmer_results(pairs)
    #blast_search(pairs, seqs)


if __name__ == '__main__':
    main()
