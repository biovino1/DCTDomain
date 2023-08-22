"""Parses pfam seed database and gets the ID and AC for each family.

__author__ = "Ben Iovino"
__date__ = "8/21/2023"
"""

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


def change_names(dct_file: str, fams: dict):
    """Changes the names of the families in the dataframe to match the AC names.

    :param df: filename of dataframe
    :param fams: dict of pfam families as key and AC as value
    """

    # Load array and rename families
    dcts = np.load(dct_file, allow_pickle=True)
    for dct in dcts:
        dct[0] = fams[dct[0]]

    # Save array
    np.save(dct_file, dcts)


def main():
    """Saves renamed families to .npy file
    """

    seed = 'pfam_data/Pfam-A.seed'
    dcts = 'pfam_data/esm2_17_875_avg.npy'
    fams = read_pfam(seed)
    change_names(dcts, fams)


if __name__ == '__main__':
    main()
