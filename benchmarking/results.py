"""This script takes the results from benchmarks.py and creates csv files with 2 fields:
    1. Homolog/NonHomologous label
    2. Bit Score

__author__ = 'Ben Iovino'
__date__ = '08/04/23'
"""

import pickle
from benchmarks import get_pairs


def main():
    """Main
    """

    # Read each dict from results file and create csv file
    dataset = 'pfam_localpfam_nomax50'
    max50_pairs = f'pfam_data/{dataset}.pair'
    pairs = get_pairs(max50_pairs)
    results = ['blast', 'csblast', 'fasta', 'phmmer', 'ublast', 'usearch']
    for result in results:
        with open(f'benchmarking/results/{dataset}/{result}_results.pkl', 'rb') as f:
            result_dict = pickle.load(f)
        with open(f'benchmarking/results/{dataset}/{result}_results.csv',
                   'w', encoding='utf8') as f:
            f.write('Homolog, Bit Score\n')
            for pair in pairs:
                if pair[2].startswith('hom'):
                    label = 1
                else:
                    label = 0
                try:
                    f.write(f'{label}, {result_dict[(pair[0], pair[1])][0]}\n')
                except IndexError:
                    f.write(f'{label}, 0\n')


if __name__ == '__main__':
    main()
