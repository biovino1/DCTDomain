"""This script takes all pairs from pfam_max50.pair and performs homology searchs on each one, with
the first sequence acting as the 'sequence database' and the second sequence acting as the query
sequence. E-value and bit scores are saved to a file.

__author__ = 'Benjamin Iovino'
__date__ = '07/27/23'
"""

import argparse
import os
import pickle
import subprocess
from Bio import SeqIO


def get_pairs(file: str) -> list:
    """This function takes a file of pairs and returns a list of tuples.

    :param file: file of pairs
    :return list: list of tuples
    """

    pairs = []
    with open(file, 'r', encoding='utf8') as f:
        for line in f:
            pairs.append(tuple(line.strip().split()))

    return pairs[1:]  # Skip header


def get_seqs(file: str) -> dict:
    """This function takes a fasta file and returns a dictionary of sequences and their IDs.

    :param file: fasta file
    :return dict: dictionary of sequences
    """

    # Read each line and add to dictionary
    seqs = {}
    with open(file, 'r', encoding='utf8') as f:
        for seq in SeqIO.parse(f, 'fasta'):
            seqs[seq.id] = str(seq.seq)

    return seqs


def phmmer_search(pairs: list, seqs: dict):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'sequence database' and 
    the second sequence is used as the query sequence in a phmmer search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results = {}
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open('bm_data/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open('bm_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from phmmer search
        # 14th line of stdout, will be blank if no hits detected
        result = subprocess.getoutput('phmmer bm_data/db_seq.fa bm_data/query_seq.fa')
        try:
            result_line = result.split('\n')[14].split()
        except IndexError:
            print(result)
            print(result_line)
            break

        # If there is a hit, check if the hit is the same as the query sequence
        if result_line == [] or result_line[-1] != pair[1]:
            result_line = 0
        else:
            result_line = [result_line[0], result_line[1]]
        results[(pair[0], pair[1])] = result_line

    # Save so we can load for later parsing, much faster than running search each time
    with open('benchmarking/results/phmmer_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def blast_search(pairs: list, seqs: dict):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'blast database' and 
    the second sequence is used as the query sequence in a blastp search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results, db_seq = {}, ''
    for i, pair in enumerate(pairs):
        print(i)
        # If first seq is different, make new blast database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            os.system('rm -rf bm_data/')
            os.system('mkdir bm_data')
            with open('bm_data/db_seq.fa', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')
            os.system('makeblastdb -in bm_data/db_seq.fa '
                '-dbtype prot -parse_seqids -out bm_data/blastdb/db_seq')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open('bm_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from blastp search
        result = subprocess.getoutput('blastp -query bm_data/query_seq.fa '
                                      '-db bm_data/blastdb/db_seq')

        # 29th line of stdout, will be blank if no hits detected
        result_line = result.split('\n')[29].split()[-2:]
        if result_line == []:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    os.system('rm -rf bm_data/')
    with open('benchmarking/results/blast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def fasta_search(pairs: list, seqs: dict):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'fasta database' and 
    the second sequence is used as the query sequence in a fasta search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results = {}
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open('bm_data/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open('bm_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get bit score, SW score, and E-value from fasta search
        result = subprocess.getoutput('fasta36 bm_data/query_seq.fa bm_data/db_seq.fa')
        result_line = result.split('\n')[22].split()
        result_line = [result_line[9], result_line[11]]
        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    with open('benchmarking/results/fasta_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def ublast_search(pairs: list, seqs: dict):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a ublast search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results = {}
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open('bm_data/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open('bm_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from ublast, unfortunately have to save results to file
        os.system('usearch -ublast bm_data/query_seq.fa '
                  '-db bm_data/db_seq.fa -evalue 1e-9 '
                  '-userout bm_data/hits.txt -userfields bits+evalue')

        # Get results from file
        with open('bm_data/hits.txt', 'r', encoding='utf8') as f:
            result_line = f.read().split()
        if result_line == []:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    with open('benchmarking/results/ublast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def usearch_search(pairs: list, seqs: dict):
    """
    This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a usearch global search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results = {}
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open('bm_data/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open('bm_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from usearch, unfortunately have to save results to file
        os.system('usearch -usearch_global bm_data/query_seq.fa '
                  '-db bm_data/db_seq.fa -id 0.8 '
                  '-userout bm_data/hits.txt -userfields bits+evalue')

        # Get results from file
        with open('bm_data/hits.txt', 'r', encoding='utf8') as f:
            result_line = f.read().split()
        if result_line == []:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    with open('benchmarking/results/usearch_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def csblast_search(pairs: list, seqs: dict):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a csblast search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results, db_seq = {}, ''
    for i, pair in enumerate(pairs):
        print(i)
        # If first seq is different, make new blast database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            with open('bm_data/db_seq.fa', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')
            os.system('formatdb -t bm_data/db -i bm_data/db_seq.fa -p T -l bm_data/formatdb.log')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open('bm_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from blastp search
        result = subprocess.getoutput('csblast -i bm_data/query_seq.fa '
                                       '-d bm_data/db_seq.fa '
                                       '-D /home/ben/anaconda3/data/K4000.lib '
                                       '--blast-path /home/ben/anaconda3/envs/benchmarking/bin')

        result_line = result.split('\n')
        score_line = [s.find('Score') for s in result_line]
        # Look for non -1 value in score_line
        for j, score in enumerate(score_line):
            if score not in (1, -1):
                result_line = result_line[j+3].split()
                result_line = [result_line[1], result_line[2]]
        if len(result_line) > 3:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    os.system('rm -rf bm_data/')
    with open('benchmarking/results/csblast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def hhsearch_search(pairs: list, seqs: dict):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a hhsearch search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    results, db_seq = {}, ''
    for i, pair in enumerate(pairs):
        if i != 1846:
            continue
        # If first seq is different, make new hhsearch database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            os.system('rm -rf hhs_data/')
            os.system('mkdir hhs_data/')
            with open('hhs_data/db.fas', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')

            # Lots of commands here to make a database
            os.system('ffindex_from_fasta -s hhs_data/db_fas.ffdata '
                      'hhs_data/db_fas.ffindex hhs_data/db.fas')

            # DATABASE GOES IN BENCHMARKING FOLDER
            os.system('hhblits_omp -i hhs_data/db_fas -d benchmarking/scop40_01Mar17/scop40 '
                      '-oa3m hhs_data/db_a3m_wo_ss -n 2 -cpu 1 -v 0')
            os.system('mv hhs_data/db_a3m_wo_ss.ffindex hhs_data/db_a3m.ffindex')
            os.system('mv hhs_data/db_a3m_wo_ss.ffdata hhs_data/db_a3m.ffdata')
            os.system('ffindex_apply hhs_data/db_a3m.ffdata hhs_data/db_a3m.ffindex '
                      '-i hhs_data/db_hmm.ffindex -d hhs_data/db_hmm.ffdata '
                      '-- hhmake -i stdin -o stdout -v 0')
            os.system('cstranslate -f -x 0.3 -c 4 -I a3m '
                      '-i hhs_data/db_a3m -o hhs_data/db_cs219')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open('hhs_data/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')
        result = subprocess.getoutput('hhsearch -i hhs_data/query_seq.fa -d hhs_data/db')
        print(result)
        result_line = result.split('\n')
        score_line = [s.find('No Hit') for s in result_line]
        for j, score in enumerate(score_line):
            if score != -1:
                result_line = result_line[j+1].split()
                if result_line != []:
                    result_line = [result_line[5], result_line[3]]

        results[(pair[0], pair[1])] = result_line

    # Save for later parsing
    os.system('rm -rf hhs_data/')
    with open('benchmarking/results/hhsearch_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str,  default='hhsearch')
    args = parser.parse_args()

    # Get all pairs and seqs from max50 files
    max50_pairs = 'embedding/pfam_max50.pair'
    max50_seqs = 'embedding/pfam_max50.fasta'
    pairs = get_pairs(max50_pairs)
    seqs = get_seqs(max50_seqs)

    # Run search
    fxn = f'{args.s}_search'
    globals()[fxn](pairs, seqs)


if __name__ == '__main__':
    main()
