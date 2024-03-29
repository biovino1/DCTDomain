"""This script takes all pairs from pfam pair datasets and performs homology searchs on each one,
with the first sequence acting as the 'sequence database' and the second sequence acting as the
query sequence. E-value and bit scores are saved to a file.

__author__ = 'Ben Iovino'
__date__ = '07/27/23'
"""

import argparse
import datetime
import os
import logging
import pickle
import pyprost
import subprocess
from Bio import SeqIO

log_filename = 'logs/benchmarks.log'  #pylint: disable=C0103
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, filemode='w',
                     level=logging.INFO, format='%(message)s')


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


def prost_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'sequence database' and 
    the second sequence is used as the query sequence in a prost search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    results, db_seq = {}, ''
    total_time = 0
    for pair in pairs:

        # Get sequences and write each to file
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            db_emb = pyprost.quantSeq(db_seq)
        query_seq = seqs[pair[1]]
        query_emb = pyprost.quantSeq(query_seq)

        # Get E-value and bit score from prost search
        start = datetime.datetime.now()
        dist = pyprost.prostDistance(db_emb, query_emb)
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)
        results[(pair[0], pair[1])] = dist

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    with open(f'benchmarking/results/{dataset}/prost_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def phmmer_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'sequence database' and 
    the second sequence is used as the query sequence in a phmmer search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    direc = f'phm_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results = {}
    total_time = 0
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open(f'{direc}/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from phmmer search
        # 14th line of stdout, will be blank if no hits detected
        start = datetime.datetime.now()
        result = subprocess.getoutput(f'phmmer --max -E 1000000000 '
                                      f'{direc}/db_seq.fa {direc}/query_seq.fa')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)

        result_line = result.split('\n')
        score_line = [s.find(pair[1]) for s in result_line]
        for j, score in enumerate(score_line):
            if score > 10:  # db seq is last in line reporting bitscore/eval
                result_line = result_line[j].split()
                result_line = [result_line[1], result_line[0]]
                break

        if len(result_line) > 3:  # No hits detected
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)  #NOSONAR
    # Save so we can load for later parsing, much faster than running search each time
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/phmmer_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def blast_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'blast database' and 
    the second sequence is used as the query sequence in a blastp search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    direc = f'bls_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results, db_seq = {}, ''
    total_time = 0
    for pair in pairs:
        # If first seq is different, make new blast database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            os.system(f'rm -rf {direc}')
            os.system(f'mkdir -p {direc}')
            with open(f'{direc}/db_seq.fa', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')
            os.system(f'makeblastdb -in {direc}/db_seq.fa '
                f'-dbtype prot -parse_seqids -out {direc}/blastdb/db_seq')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from blastp search
        start = datetime.datetime.now()
        result = subprocess.getoutput(f'blastp -query {direc}/query_seq.fa '
                                      f'-db {direc}/blastdb/db_seq -evalue 1000000000')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)

        # 29th line of stdout, will be blank if no hits detected
        result_line = result.split('\n')[29].split()[-2:]
        if result_line == []:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/blast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def fasta_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'fasta database' and 
    the second sequence is used as the query sequence in a fasta search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    direc = f'fas_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results = {}
    total_time = 0
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open(f'{direc}/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get bit score, SW score, and E-value from fasta search
        start = datetime.datetime.now()
        result = subprocess.getoutput(f'fasta36 {direc}/query_seq.fa '
                                      f'{direc}/db_seq.fa -b 1000000000')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)
        result_line = result.split('\n')[22].split()
        result_line = [result_line[9], result_line[11]]
        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/fasta_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def ublast_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a ublast search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    direc = f'ubl_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results, db_seq = {}, ''
    total_time = 0
    for pair in pairs:

        # If first seq is different, make new ublast database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            os.system(f'rm -rf {direc}')
            os.system(f'mkdir -p {direc}')
            with open(f'{direc}/db_seq.fa', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')
            os.system(f'usearch -makeudb_ublast {direc}/db_seq.fa -output {direc}/db.udb')

        # Query sequence is always different
        query_seq = seqs[pair[1]]
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from ublast, unfortunately have to save results to file
        start = datetime.datetime.now()
        os.system(f'usearch -ublast {direc}/query_seq.fa '
                  f'-db {direc}/db.udb -evalue 1000000000 '
                  f'-userout {direc}/hits.txt -userfields bits+evalue')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)

        # Get results from file
        with open(f'{direc}/hits.txt', 'r', encoding='utf8') as f:
            result_line = f.read().split()
        if result_line == []:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/ublast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def usearch_search(pairs: list, seqs: dict, dataset: str):
    """
    This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a usearch global search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    direc = f'use_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results = {}
    total_time = 0
    for pair in pairs:

        # Get sequences and write each to file
        db_seq = seqs[pair[0]]
        query_seq = seqs[pair[1]]
        with open(f'{direc}/db_seq.fa', 'w', encoding='utf8') as db:
            db.write(f'>{pair[0]}\n{db_seq}')
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from usearch, unfortunately have to save results to file
        start = datetime.datetime.now()
        os.system(f'usearch -search_local {direc}/query_seq.fa '
                  f'-db {direc}/db_seq.fa -evalue 1000000000 '
                  f'-userout {direc}/hits.txt -userfields bits+evalue')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)

        # Get results from file
        with open(f'{direc}/hits.txt', 'r', encoding='utf8') as f:
            result_line = f.read().split()
        if result_line == []:
            result_line = [0]
        else:
            result_line = [result_line[0], result_line[1]]
        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/usearch_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def csblast_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a csblast search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    """

    direc = f'csb_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results, db_seq = {}, ''
    total_time = 0
    for pair in pairs:

        # If first seq is different, make new blast database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            with open(f'{direc}/db_seq.fa', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')
            os.system(f'formatdb -t {direc}/db -i {direc}/db_seq.fa -p T -l {direc}/formatdb.log')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Get E-value and bit score from blastp search
        start = datetime.datetime.now()
        result = subprocess.getoutput(f'csblast -i {direc}/query_seq.fa '
                                       f'-d {direc}/db_seq.fa '
                                       '-D /home/ben/anaconda3/data/K4000.lib '
                                       '--blast-path $CONDA_PREFIX/bin '
                                       '-e 1000000000')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)

        # Find line with bit score and E-value
        result_line = result.split('\n')
        try:  # No hits occured and result_line[j] is out of range
            score_line = [s.find(pair[0]) for s in result_line]
            for j, score in enumerate(score_line):
                if score == 0:  # db seq is first index in line reporting bitscore/eval
                    result_line = result_line[j].split()
                    result_line = [result_line[1], result_line[2]]
        except IndexError:
            result_line = [0]
        try:  # No hits occured and wrong line was selected
            if len(result_line) > 3 or 'e' in result_line[0] or result_line[0] == '[blastpgp]':
                result_line = [0]
        except TypeError:
            result_line = [0]
        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/csblast_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def hhsearch_search(pairs: list, seqs: dict, dataset: str):
    """This function takes a list of protein pairs, each pair being used to get their respective
    sequences from a dictionary of seqs. The first sequence is used as the 'database' and the second
    sequence is used as the query sequence in a hhsearch search.

    :param pairs: list of pairs
    :param seqs: dictionary where seq is ID and value is sequence
    :param dataset: dataset name
    """

    direc = f'hh_data/{dataset}'
    os.system(f'rm -rf {direc}')
    os.system(f'mkdir -p {direc}')
    results, db_seq = {}, ''
    total_time = 0
    for pair in pairs:
        # If first seq is different, make new hhsearch database
        if db_seq != seqs[pair[0]]:
            db_seq = seqs[pair[0]]
            os.system(f'rm -rf {direc}')
            os.system(f'mkdir -p {direc}')
            with open(f'{direc}/db.fas', 'w', encoding='utf8') as db:
                db.write(f'>{pair[0]}\n{db_seq}')

            # Lots of commands here to make a database
            os.system(f'ffindex_from_fasta -s {direc}/db_fas.ffdata '
                      f'{direc}/db_fas.ffindex {direc}/db.fas')

            # DATABASE GOES IN BENCHMARKING FOLDER
            os.system(f'hhblits_omp -i {direc}/db_fas -d benchmarking/scop40_01Mar17/scop40 '
                      f'-oa3m {direc}/db_a3m -n 2 -cpu 1 -v 0')
            os.system(f'ffindex_apply {direc}/db_a3m.ffdata {direc}/db_a3m.ffindex '
                      f'-i {direc}/db_hmm.ffindex -d {direc}/db_hmm.ffdata '
                      '-- hhmake -i stdin -o stdout -v 0')
            os.system('cstranslate -f -x 0.3 -c 4 -I a3m '
                      f'-i {direc}/db_a3m -o {direc}/db_cs219')

        # Query sequence is always different so write to file
        query_seq = seqs[pair[1]]
        with open(f'{direc}/query_seq.fa', 'w', encoding='utf8') as query:
            query.write(f'>{pair[1]}\n{query_seq}')

        # Record total time to run search
        start = datetime.datetime.now()
        result = subprocess.getoutput(f'hhsearch -i {direc}/query_seq.fa '
                                       f'-d {direc}/db -E 1000000000')
        end = datetime.datetime.now()
        total_time += (end-start).total_seconds()
        logging.info('%s %s %s', pair[0], pair[1], end-start)

        result_line = result.split('\n')
        score_line = [s.find('No Hit') for s in result_line]
        for j, score in enumerate(score_line):
            if score != -1:
                result_line = result_line[j+1].split()
                if result_line != []:
                    result_line = [result_line[5], result_line[3]]

        results[(pair[0], pair[1])] = result_line

    logging.info('Total time: %s', total_time)
    # Save for later parsing
    os.system(f'rm -rf {direc}')
    with open(f'benchmarking/results/{dataset}/hhsearch_results.pkl', 'wb') as f:
        pickle.dump(results, f)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, default='pfam_max50')
    parser.add_argument('-s', type=str,  default='prost')
    args = parser.parse_args()

    # Get all pairs and seqs from max50 files
    pairs = f'pfam_data/{args.d}.pair'
    seqs = f'pfam_data/{args.d}.fasta'
    pairs = get_pairs(pairs)
    seqs = get_seqs(seqs)

    # Run search
    fxn = f'{args.s}_search'
    globals()[fxn](pairs, seqs, args.d)


if __name__ == '__main__':
    main()
