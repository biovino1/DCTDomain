"""This script takes a .npy file containing DCT vectors and returns a dataframe of the Pfam families
and their DCT vectors

__author__ = "Ben Iovino"
__date__ = "06/14/23"
"""

import pickle
import pandas as pd
import numpy as np


def make_matrix(filename: str):
    """=============================================================================================
    This function accepts a .npz file, parses for Pfam family names and DCT vectors, and then saves
    it as a dataframe.

    :param filename: .npz file containing DCT vectors
    ============================================================================================="""

    # Load data from .npy file
    dcts = np.load(filename, allow_pickle=True)

    # Get the DCT vectors and domain IDs
    families = []
    features = []
    for i in range(len(dcts)):  #pylint: disable=C0200
        dct = dcts[i]
        families.append(dct[0])
        features.append(dct[1])

    # Make pandas dataframe
    dct_df = pd.DataFrame(features,
                            index=[str(i.split('.')[0]) for i in families],
                            columns=[str(i) for i in range(600)])
    dct_df = dct_df.sort_index()  # Reorder by family ID

    # Save dataframe
    with open('pfam_data/searchemb_df.pkl', 'wb') as file:
        pickle.dump(dct_df, file)


def main():
    """=============================================================================================
    Main initializes a filename for .npz file and calls make_matrix to create a dataframe.
    ============================================================================================="""

    filename = "pfam_data/esm2_17_875_avg.npy"
    make_matrix(filename)


if __name__ == '__main__':
    main()
