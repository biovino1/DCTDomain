"""================================================================================================
This script takes two dataframes and compares their values.

Ben Iovino  06/14/23   DCTDomain
================================================================================================"""

import pickle
import pandas as pd


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


def main():

    # Clans df is used for comparing to cluster df
    df1, df2 = 'data/clans_df.pkl', 'data/cluster_df.pkl'
    clans_df, cluster_df = load_dfs(df1, df2)
    print(clans_df.shape, cluster_df.shape)



if __name__ == '__main__':
    main()
