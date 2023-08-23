"""================================================================================================
This script takes a .npz file containing DCT vectors and returns a dataframe of the domain ids and
their DCT vectors.

Ben Iovino  06/21/23   DCTDomain
================================================================================================"""

import pickle
import pandas as pd
import numpy as np


def make_matrix(filename: str):
    """=============================================================================================
    This function accepts a .npz file, parses for SCOP domain ids and DCT vectors, and then saves
    it as a dataframe.

    :param filename: .npz file containing DCT vectors
    ============================================================================================="""

    # Load data from .npz file
    data_all0 = np.load(filename)
    data_all = data_all0['arr_0']  # np.array of np.arrays with 3 elements: domain, range, and DCT vector
    data_all = data_all[np.argsort(data_all[:, 0])]  # Sort to match with family dfs
    tot = len(data_all)  # 36544 domains

    # Get the DCT vectors and domain IDs
    domains = []
    features = []
    for i in range(tot):
        data = data_all[i]
        domains.append(data[0])
        dct = np.fromstring(data[2], sep=",")  # pull out DCT vector
        features.append(dct)

    # Make pandas dataframe
    dct_df = pd.DataFrame(features,
                            index=[str(i.split('.')[0]) for i in domains],
                            columns=[str(i) for i in range(475)])

    # Save dataframe
    with open('scop_data/dct_df.pkl', 'wb') as file:
        pickle.dump(dct_df, file)


def main():
    """=============================================================================================
    Main initializes a filename for .npz file and calls make_matrix to create a dataframe.
    ============================================================================================="""

    filename = 'scop_data/scop_fa_represeq_lib_latest.npz'
    make_matrix(filename)


if __name__ == '__main__':
    main()
