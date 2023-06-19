"""================================================================================================
This script clusters the vectors in dct_df.pkl and returns a dataframe with the clusters
assignments.

Ben Iovino  06/15/23   DCTDomain
================================================================================================"""

import argparse
import logging
import pickle
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

logging.basicConfig(filename='data/clustering.log',
                     level=logging.INFO, format='%(message)s')


def get_linkage(dct_df: pd.DataFrame, args: argparse.Namespace):
    """=============================================================================================
    This function accepts a dataframe of DCT vectors and returns the linkage data.

    :param dct_df: dataframe of DCT vectors
    :param args: argparse namespace
    ============================================================================================="""

    # Load linkage file if given as argument
    if args.l != 'None':
        with open(args.l, 'rb') as file:
            linkage_data = pickle.load(file)
    else:
        linkage_data = linkage(dct_df, method=args.m, metric=args.p)
        with open('data/linkage_data.pkl', 'wb') as file:
            pickle.dump(linkage_data, file)

    return linkage_data


def get_clusters(linkage_data: np.ndarray, args: argparse.Namespace, dct_df: pd.DataFrame):
    """=============================================================================================
    This function accepts linkage data from scipy's linkage function and returns a dataframe with
    the cluster assignments for each family in dct_df.

    :param linkage_data: linkage data
    :param dct_df: dataframe of DCT vectors
    :param args: argparse namespace
    ============================================================================================="""

    clusters = fcluster(linkage_data, args.t, criterion=args.c)

    # Rename dct_families to be in increasing order so we can assign values to an array
    dct_families = dct_df.index
    dct_fam_ints = {}
    for i, fam in enumerate(dct_families):
        dct_fam_ints[fam] = i

    # Indices of clusters should match the indices of the families
    # Make a dictionary of clusters with cluster as key and list of families as value
    cluster_dict = {}
    for i in range(len(clusters)):  #pylint: disable=C0200
        if clusters[i] not in cluster_dict:
            cluster_dict[clusters[i]] = [list(dct_fam_ints.values())[i]]
        else:
            cluster_dict[clusters[i]].append(list(dct_fam_ints.values())[i])
    logging.info('Number of clusters: %s', len(cluster_dict))

    # For each cluster, we want to assign a 1 to the cell between each family in the cluster
    cluster_df = pd.DataFrame(0, index=dct_families, columns=dct_families)
    for cluster in cluster_dict:  #pylint: disable=C0206
        for family1 in cluster_dict[cluster]:
            row = [0 for i in range(len(dct_df))]  # Family1 relations to others
            for family2 in cluster_dict[cluster]:
                row[family2] = 1  # If families in same cluster, set to 1

            # Convert family integer to family name after we indexed the row
            # Place whole row in df
            family1_name = list(dct_fam_ints.keys())[list(dct_fam_ints.values()).index(family1)]
            cluster_df.loc[family1_name] = row

    # Save cluster_df as pickle file
    cluster_df.to_pickle('data/cluster_df.pkl')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, help='linkage method', default='ward')
    parser.add_argument('-p', type=str, help='linkage metric', default='euclidean')
    parser.add_argument('-l', type=str, help='linkage data file', default='data/linkage_data.pkl')
    parser.add_argument('-t', type=float, help='fcluster threshold', default=500)
    parser.add_argument('-c', type=str, help='fcluster criterion', default='distance')
    args = parser.parse_args()
    logging.info(args)

    # Open dct_df.pkl and get linkage_data
    with open('data/dct_df.pkl', 'rb') as file:
        dct_df = pickle.load(file)
    linkage_data = get_linkage(dct_df, args)

    # Get clusters from linkage data and save as pickle file
    get_clusters(linkage_data, args, dct_df)


if __name__ == '__main__':
    main()
