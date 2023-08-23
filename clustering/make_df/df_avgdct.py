"""This script takes a .npy file containing DCT vectors and returns a dataframe of the Pfam families
and their DCT vectors

__author__ = "Ben Iovino"
__date__ = "06/14/23"
"""

import pickle
import pandas as pd
import numpy as np


def read_pfam(file: str) -> dict:
    """Returns a dict of pfam families as key and AC as value.

    :param file: filename of pfam seed database
    :return dict: pfam families as key and AC as value
    """

    fams = {}
    with open(file, 'r', encoding='utf8', errors='replace') as file:
        for line in file:
            if line.startswith('#=GF ID'):
                fams[line.split()[2]] = ''
            elif line.startswith('#=GF AC'):
                fams[list(fams.keys())[-1]] = line.split()[2]

    return fams


def change_names(dct_file: str, fams: dict) -> np.ndarray:
    """Changes the names of the families in the dataframe to match the AC names.

    :param df: filename of dataframe
    :param fams: dict of pfam families as key and AC as value
    return: array of DCT vectors with renamed families
    """

    # Load array and rename families
    dcts = np.load(dct_file, allow_pickle=True)
    for dct in dcts:
        dct[0] = fams[dct[0]]

    return dcts


def make_matrix(dcts: np.ndarray):
    """This function accepts an array, parses for Pfam family names and DCT vectors, and then
    saves it as a dataframe.

    :param filename: .npz file containing DCT vectors
    """

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
    with open('pfam_data/searchemb_df.pkl', 'wb') as file:
        pickle.dump(dct_df, file)


def main():
    """Changes names in .npy file and saves it as a pandas df.
    """

    seed = 'pfam_data/Pfam-A.seed'
    dcts = "pfam_data/esm2_17_875_avg.npy"
    fams = read_pfam(seed)
    dcts = change_names(dcts, fams)
    make_matrix(dcts)


if __name__ == '__main__':
    main()
