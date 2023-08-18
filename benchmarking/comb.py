"""This script combines localpfam homologous pairs with non-hom pairs from both
the nomax50 dataset.

__author__ = 'Ben Iovino'
__date__ = '08/17/23'
"""


def clean_nomax_results(nomax50_files: list):
    """I messed up and ran benchmarks.py on pfam_nomax50.pair instead of pfam_nomax50.pair.found.
    This bad boy will fix it right up.

    :param nomax50_files: list of files to clean
    """

    # Find differences between .pair and .pair.found
    with open(nomax50_files[0], 'r', encoding='utf8') as f:
        pair = f.readlines()
    with open(nomax50_files[0].replace('.found', ''), 'r', encoding='utf8') as f:
        pair_found = f.readlines()
    diff = set(pair_found) - set(pair)

    # Find positions of differences
    diff_pos = []
    for line in diff:
        diff_pos.append(pair_found.index(line))

    diff_pos = sorted(diff_pos)

    # Remove differences from pair
    #for i, pos in enumerate(diff_pos):
        #pair_found.pop(pos-i)

    # Open results files and remove positions
    for i in range(2, len(nomax50_files)):
        with open(nomax50_files[i], 'r', encoding='utf8') as f:
            results = f.readlines()
        for j, pos in enumerate(diff_pos):
            results.pop(pos-j)
        with open(nomax50_files[i], 'w', encoding='utf8') as f:
            f.writelines(results)


def main():

    localpfam_files = ['pfam_data/pfam_localpfam_nomax50.pair',
                        'pfam_data/pfam_localpfam_nomax50-dctsim.txt',
                        'benchmarking/results/pfam_localpfam_nomax50/blast_results.csv',
                        'benchmarking/results/pfam_localpfam_nomax50/fasta_results.csv',
                        'benchmarking/results/pfam_localpfam_nomax50/hhsearch_results.csv',
                        'benchmarking/results/pfam_localpfam_nomax50/phmmer_results.csv',
                        'benchmarking/results/pfam_localpfam_nomax50/ublast_results.csv',
                        'benchmarking/results/pfam_localpfam_nomax50/usearch_results.csv',
                        'benchmarking/results/pfam_localpfam_nomax50/csblast_results.csv']

    nomax50_files = ['pfam_data/pfam_nomax50.pair.found',
                      'pfam_data/pfam_nomax50-dctsim.txt',
                      'benchmarking/results/pfam_nomax50/blast_results.csv',
                      'benchmarking/results/pfam_nomax50/fasta_results.csv',
                      'benchmarking/results/pfam_nomax50/hhsearch_results.csv',
                      'benchmarking/results/pfam_nomax50/phmmer_results.csv',
                      'benchmarking/results/pfam_nomax50/ublast_results.csv',
                      'benchmarking/results/pfam_nomax50/usearch_results.csv',
                      'benchmarking/results/pfam_nomax50/csblast_results.csv']

    #clean_nomax_results(nomax50_files)

    # Combine localpfam and nomax50 pairs
    for i in range(len(localpfam_files)):  #pylint: disable=C0200
        with open(nomax50_files[i], 'r', encoding='utf8') as f:
            nomax50 = f.readlines()[36279:]
        with open(localpfam_files[i], 'r', encoding='utf8') as f:
            localpfam = f.readlines()

            # Combine lists, starting with localpfam
            combined = localpfam + nomax50

        # Write combined list to file
        with open(localpfam_files[i], 'w', encoding='utf8') as f:
            f.writelines(combined)


if __name__ == '__main__':
    main()
