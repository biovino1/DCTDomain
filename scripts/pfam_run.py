"""================================================================================================
This script runs cluster_dct.py with set parameters and then runs compare_dfs.py to compare the
clustering results to the ground truth pfam clans with compare_dfs.py.

Ben Iovino  06/20/23   DCTDomain
================================================================================================"""

import os


def main():
    """=============================================================================================
    Main checks for linkage data and calls cluster_dct.py and compare_dfs.py for each threshold.
    ============================================================================================="""

    # Check if there is linkage data
    if not os.path.exists('pfam_data/linkage_data.pkl'):
        linkage_data = 'None'
    else:
        linkage_data = 'pfam_data/linkage_data.pkl'

    linkage = 'ward'
    metric = 'euclidean'
    thresholds = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    criterion = 'distance'

    for t in thresholds:
        os.system('python scripts/pfam_cluster.py ' \
                  f'-m {linkage} -p {metric} -l {linkage_data} -t {t} -c {criterion}')
        os.system('python scripts/pfam_compare.py')


if __name__ == '__main__':
    main()
