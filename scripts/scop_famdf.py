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

    domids, folds, superfams, fams = [], {}, {}, {}
    with open(classifications, 'r', encoding='utf8') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue

            # Get representative domain ID and add to respective superfam ID and fam ID
            row = ''.join(row).split()
            row10 = row[10].split(',')
            domid, fold, superfam, fam = row[0], row10[2], row10[3], row10[4]
            if domid in domids:  # Skip if domid already in list, its a diff region of same domain
                continue
            domids.append(domid)
            if fold not in folds:  # Update fold dict or append to list
                folds[fold] = [domid]
            else:
                folds[fold].append(domid)
            if superfam not in superfams:  # Update superfam dict or append to list
                superfams[superfam] = [domid]
            else:
                superfams[superfam].append(domid)
            if fam not in fams:  # Update fam dict or append to list
                fams[fam] = [domid]
            else:
                fams[fam].append(domid)

    # Sort domids and create new sorted dicts for fams and superfams
    domids.sort()

    print(folds)
    return domids, folds, superfams, fams


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

    # Rename domids to be in increasing order so we can assign values to an array
    domid_int = {}
    for i, ind in enumerate(domids):
        domid_int[ind] = i

    # Fill matrix by going through each member in dict
    for fam, members in classifier.items():  #pylint: disable=W0612
        for domid1 in members:
            row = [0 for i in range(len(matrix))]  #domid1 relation to all other domids
            for domid2 in members:
                row[domid_int[domid2]] = 1  # Assign 1 to corresponding domid2 index in row
            matrix.loc[domid1] = row  # Update row in matrix

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
    domids, folds, superfams, fams = get_families(cla)
    make_matrix(domids, folds, 'folds')
    make_matrix(domids, superfams, 'superfams')
    make_matrix(domids, fams, 'fams')
    


if __name__ == '__main__':
    main()
