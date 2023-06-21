"""================================================================================================
This script takes two dataframes and compares their values.

Ben Iovino  06/14/23   DCTDomain
================================================================================================"""

import argparse
import logging
import pickle
import pandas as pd

logging.basicConfig(filename='clustering.log',
                     level=logging.INFO, format='%(message)s')


def load_dfs(df1: str, df2: str) -> pd.DataFrame:
    """=============================================================================================
    This function accepts two filenames and loads them as dataframes. It also makes sure they
    have matching rows and columns.

    :param df1: filename of first dataframe
    :param df2: filename of second dataframe
    :return pd.DataFrame: clans dataframe, cluster dataframe
    ============================================================================================="""

    # Load dataframes
    with open(df1, 'rb') as file:
        clans_df = pickle.load(file)
    with open(df2, 'rb') as file:
        cluster_df = pickle.load(file)

    # Make sure they have matching rows and columns
    # Compare column names from df and clans dataframes to see which family is missing in df
    cluster_families = cluster_df.index
    clan_families = clans_df.index
    clan_missing = [i for i in clan_families if i not in cluster_families]

    # Remove these families from clans dataframe, both indices and columns
    clans_df = clans_df.drop(clan_missing, axis=1)
    clans_df = clans_df.drop(clan_missing, axis=0)

    return clans_df, cluster_df


def compare_dfs(df1: pd.DataFrame, df2: pd.DataFrame):
    """=============================================================================================
    This function accepts two dataframes of equal shape and compares their values. The first df
    is used as the ground truth, and the second df is used as the test df.

    :param df1: filename of ground truth dataframe
    :param df2: filename of test dataframe
    ============================================================================================="""

    # Counts for ground truth df and counts for comparison
    # 19621 fams compared to themselves so subtract that from counts that include 1
    all_zeros, all_ones = 0, 0
    true_zeros, true_ones, false_zeros, false_ones = 0, 0, 0, 0

    # Iterate through each family
    for i, family in enumerate(df1):
        clans_arr = df1[family].values  # Ground truth array
        cluster_arr = df2[family].values  # Test array
        for j in range(len(clans_arr)):  #pylint: disable=C0200
            if i == j: # Ignore self-comparison
                continue

            # By checking clans_df values to be either 0 or 1, we are ignoring NaN values
            if clans_arr[j] == 1:
                all_ones += 1
                if cluster_arr[j] == 1:  # True positive
                    true_ones += 1
                if cluster_arr[j] == 0:  # False negative
                    false_zeros += 1
            if clans_arr[j] == 0:
                all_zeros += 1
                if cluster_arr[j] == 0:  # True negative
                    true_zeros += 1
                if cluster_arr[j] == 1:  # False positive
                    false_ones += 1

    # Report confusion matrix
    logging.info('True Positives: %s', true_ones)
    logging.info('True Negatives: %s', true_zeros)
    logging.info('False Positives: %s', false_ones)
    logging.info('False Negatives: %s', false_zeros)

    # Calculate accuracy, precision, recall, and f1 score
    accuracy = (true_ones + true_zeros) / (all_ones + all_zeros)
    precision = true_ones / (true_ones + false_ones)
    recall = true_ones / (true_ones + false_zeros)
    f1 = 2 * (precision * recall) / (precision + recall)

    # Report accuracy, precision, recall, and f1 score
    logging.info('Accuracy: %s', round(accuracy, 3))
    logging.info('Precision: %s', round(precision, 3))
    logging.info('Recall: %s', round(recall, 3))
    logging.info('F1 Score: %s\n', round(f1, 3))


def main():
    """=============================================================================================
    Main initializes path to two dataframes and loads them with load_dfs(). It then compares their
    values with compare_dfs().
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='data used', default='scop_data')
    parser.add_argument('-c', type=str, help='comparison data', default='fams_df.pkl')
    args = parser.parse_args()

    # Comparison df is used for comparing to cluster assignments in cluster_df
    df1, df2 = f'{args.d}/{args.c}', f'{args.d}/cluster_df.pkl'
    clans_df, cluster_df = load_dfs(df1, df2)

    # Compare values
    compare_dfs(clans_df, cluster_df)


if __name__ == '__main__':
    main()
