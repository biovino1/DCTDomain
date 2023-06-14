"""================================================================================================
This script takes the Pfam-A.clans file and returns a matrix that indicates which sequences
are in the same clan. 

Ben Iovino  06/12/23   DCTDomain
================================================================================================"""

import csv
import pandas as pd
import os


def get_families(pfam: str) -> list:
    """=============================================================================================
    This function accepts a pfam database file and parses it into a list of all families in the
    database. It also returns a list of families that belong to a clan.

    :param pfam: Pfam database file
    :return list: all families in the database
    ============================================================================================="""

    # Read clans db with csv reader
    families, clan_fams = [], []
    with open(pfam, 'r', encoding='utf8', errors='replace') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            families.append(row[0])
            if row[1]:
                clan_fams.append(row[0])

    return families, clan_fams


def read_clans(pfam: str) -> dict:
    """=============================================================================================
    This function accepts a pfam database file and parses it into a dictionary. Each key is a
    clan and its value is a list of families in that clan. The dict is saved as a pickle file.

    :param pfam: Pfam database file
    :return dict: clans as keys with list of families in that clan as the value
    ============================================================================================="""

    # Read clans db with csv reader
    clans = {}
    with open(pfam, 'r', encoding='utf8', errors='replace') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            clan = row[1]
            if clan:  # If this family is a part of a clan, add family to clan dict
                if clan in clans:
                    clans[clan].append(row[0])
                else:  # Initialize new key
                    clans[clan] = [row[0]]

    return clans


def make_matrix(families, clan_fams, clans):
    """=============================================================================================
    This function accepts a list of pfam families and a dictionary of clans and their families. It
    creates a len(families) x len(families) matrix where each value indicates if the two families
    belong to the same clan (1) or not (0). If one of the families does not belong to a clan, it
    will be indicated by a NaN.

    :param families: list of families in Pfam clans db
    :param clan_fams: list of families that belong to a clan
    :param clans: dict of clans as keys with list of families in that clan as the value
    ============================================================================================="""

    # Create len(fam) x len(fam) matrix
    matrix = pd.DataFrame(index=families, columns=families)

    # Fill matrix with 0s where both families are in a clan
    for family1 in clan_fams:
        for family2 in clan_fams:
            matrix.loc[family1, family2] = 0

    # Fill matrix by going through each clan
    for key, value in clans.items():  #pylint: disable=W0612
        for family1 in value:
            for family2 in value:
                matrix.loc[family1, family2] = 1

    # Save matrix as pickle file
    matrix.to_pickle('data/clans_df.pkl')


def main():
    """=============================================================================================
    Main initialies pfam db file, gets all families and clans with get_families() and read_clans(),
    then creates a matrix of families with values indicating if they are in the same clan with
    make_matrix().
    ============================================================================================="""

    # Read Pfam-A.seed if it exists
    if os.path.exists('Data/Pfam-A.clans.tsv'):
        pfam = 'Data/Pfam-A.clans.tsv'
    else:
        print('Pfam-A.seed not found. Downloading from Pfam...')
        if not os.path.exists('Data'):
            os.mkdir('Data')
        os.system('wget -P Data ' \
            'https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.clans.tsv.gz')
        os.system('gunzip Data/Pfam-A.clans.tsv.gz')
        pfam = 'Data/Pfam-A.clans.tsv'

    # Get families and clans
    families, clan_fams = get_families(pfam)
    clans = read_clans(pfam)

    # Create matrix of families with values indicating if they are in the same clan
    make_matrix(families, clan_fams, clans)


if __name__ == '__main__':
    main()
