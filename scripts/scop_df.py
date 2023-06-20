"""================================================================================================
This script takes the scop-cla-latest.txt file and returns two matrices indicating which domains
are in the same family and superfamily.

Ben Iovino  06/20/23   DCTDomain
================================================================================================"""

import csv
import os
import pandas as pd


def get_families(classifications: str) -> tuple:
    """=============================================================================================
    This function accepts a scop classification file and returns a list of domain names, a dict of
    families as keys with a list of domains in that family as the value, and a dict of superfamilies
    as keys with a list of domains in that superfamily as the value.

    :param superfams: SCOP superfamilies file
    :return tuple: list of domain ids, family dict, and superfamily dict
    ============================================================================================="""

    domids, families, superfamilies = [], {}, {}
    with open(classifications, 'r', encoding='utf8') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue

            # Get representative domain ID and add to respective superfam ID and fam ID
            row = ''.join(row).split()
            domid, superfam, fam = row[0], row[10].split(',')[3], row[10].split(',')[4]
            domids.append(domid)
            if fam not in families:  # Update fam dict or append to list
                families[fam] = [domid]
            else:
                families[fam].append(domid)
            if superfam not in superfamilies:  # Update superfam dict or append to list
                superfamilies[superfam] = [domid]
            else:
                superfamilies[superfam].append(domid)

    return domids, families, superfamilies


def make_matrix(domids: list, classifier: dict, label: str):
    """=============================================================================================
    This function accepts a list of scop domain representatives and a dictionary that classifies
    them into families or superfamilies. It returns a matrix that indicates which domains are in
    the same class

    :param domids: list of domain names
    :param classifier: dict of families or superfamilies as keys with list of domains in that class
    :param label: string indicating whether classifier is family or superfamily
    ============================================================================================="""

    # Create len(domids) x len(domids) matrix
    matrix = pd.DataFrame(0, index=domids, columns=domids)

    # Fill matrix by going through each member in dict
    for key, value in classifier.items():  #pylint: disable=W0612
        for family1 in value:
            for family2 in value:
                matrix.loc[family1, family2] = 1

    # Save matrix as pickle file
    matrix.to_pickle(f'scop_data/{label}_df.pkl')


def main():
    """=============================================================================================
    Main downloads scop classification file if it does not exist, then calls get_families() to
    get pdbids, families, and superfamilies. Then, it calls make_matrix() to create matrices of
    pdbid relations.
    ============================================================================================="""

    # Initialize scop classification file if it exists
    if os.path.exists('scop_data/scop-cla-latest.txt'):
        cla = 'scop_data/scop-cla-latest.txt'
    else:
        print('SCOP classification file not found. Downloading from SCOP...')
        if not os.path.exists('scop_data'):
            os.mkdir('scop_data')
        os.system('wget -P scop_data ' \
            'https://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt')
        cla = 'scop-cla-latest.txt'

    # Read cla file and create matrices
    domids, families, superfamilies = get_families(cla)
    make_matrix(domids, families, 'fams')
    make_matrix(domids, superfamilies, 'superfams')


if __name__ == '__main__':
    main()
