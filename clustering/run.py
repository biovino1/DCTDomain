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

    data = 'scop_data'

    # Check if there is linkage data
    if not os.path.exists(f'{data}/linkage_data.pkl'):
        linkage_data = 'None'
    else:
        linkage_data = f'{data}/linkage_data.pkl'

    linkage = 'ward'
    metric = 'euclidean'
    thresholds = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
    criterion = 'distance'
    comparison = 'fams_df.pkl'

    for t in thresholds:
        os.system('python scripts/cluster.py ' \
                  f'-d {data} -m {linkage} -p {metric} -l {linkage_data} -t {t} -c {criterion}')
        os.system('python scripts/compare.py ' \
                  f'-d {data} -c {comparison}')


if __name__ == '__main__':
    main()
