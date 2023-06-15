"""================================================================================================
This script takes two dataframes and compares their values.

Ben Iovino  06/14/23   DCTDomain
================================================================================================"""

import pickle
import pandas as pd
import numpy as np


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

    zero_count, one_count, total_count = 0, -19621, -19621
    for family in df1:
        clans_arr = df1[family].values  # Ground truth array
        cluster_arr = df2[family].values  # Test array
        for i in range(len(clans_arr)):  #pylint: disable=C0200
            if clans_arr[i] == 1 and cluster_arr[i] == 1:
                one_count += 1
            if clans_arr[i] == 0 and cluster_arr[i] == 0:
                zero_count += 1

        # For accuracy calculation, ignore NaN values because they don't belong to clans
        total_count += (len(clans_arr) - np.count_nonzero(pd.isna(clans_arr)))

    # Remove the number of families compared to themselves
    similarity = (one_count) / (total_count)
    print(zero_count, one_count, total_count)
    print(similarity)


def main():
    """=============================================================================================
    Main initializes path to two dataframes and loads them with load_dfs(). It then compares their
    values with compare_dfs().
    ============================================================================================="""

    # Clans df is used for comparing to cluster df
    df1, df2 = 'data/clans_df.pkl', 'data/cluster_df.pkl'
    clans_df, cluster_df = load_dfs(df1, df2)

    # Compare values
    compare_dfs(clans_df, cluster_df)


if __name__ == '__main__':
    main()
