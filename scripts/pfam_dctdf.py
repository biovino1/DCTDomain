"""================================================================================================
This script takes a .npz file containing DCT vectors and returns a dataframe of the Pfam families
and their DCT vectors

Ben Iovino  06/14/23   DCTDomain
================================================================================================"""

import pickle
import pandas as pd
import numpy as np


def make_matrix(filename: str):
    """=============================================================================================
    This function accepts a .npz file, parses for Pfam family names and DCT vectors, and then saves
    it as a dataframe.

    :param filename: .npz file containing DCT vectors
    ============================================================================================="""

    # Load data from .npz file
    data_all0 = np.load(filename)
    data_all = data_all0['arr_0']  # np.array of np.arrays with 3 elements: domain, range, and DCT vector
    tot = len(data_all)  # 19621 domains

    # Get the DCT vectors and domain IDs
    families = []
    features = []
    for i in range(tot):
        data = data_all[i]

        #si: pfam domain ID, domstr: range, qstr: DCT vector
        (si, domstr, qstr) = data[:] #pylint: disable=W0612
        ai = np.fromstring(qstr, sep=",")[-475:] # vector of 475 numbers
        families.append(si)
        features.append(ai)

    # Make pandas dataframe
    dct_df = pd.DataFrame(features, index=[str(i.split('.')[0]) for i in families], columns=[str(i) for i in range(475)])
    dct_df = dct_df.sort_index()  # Reorder by family ID

    # Save dataframe
    with open('pfam_data/dct_df.pkl', 'wb') as file:
        pickle.dump(dct_df, file)


def main():
    """=============================================================================================
    Main initializes a filename for .npz file and calls make_matrix to create a dataframe.
    ============================================================================================="""

    filename = "pfam_data/Pfam-A-cons-domain.npz"
    make_matrix(filename)


if __name__ == '__main__':
    main()
