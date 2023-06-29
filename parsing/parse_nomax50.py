"""================================================================================================
This script takes a file with protein sequence ID's and their pfam domains. It returns the ID's
which share at least one domain, but not all domains.

Ben Iovino  06/28/23   DCTDomain
================================================================================================"""

import pickle


def read_file(file: str) -> dict:
    """=============================================================================================
    This function accepts a file path and parses each line. The first element is the sequence ID and
    the second is the domains associated with that sequence.

    :param file: file path
    :return dict: sequence ID is key, value is a list containing the domains.
    ============================================================================================="""

    # Read each line and add to dict
    prots = {}
    with open(file, 'r', encoding='utf8') as f:
        for line in f:
            line = line.split()
            prots[line[0]] = line[1].split(';')

    return prots


def compare_dicts(count1: dict, count2: dict, prot1: str, prot2: str):
    """=============================================================================================
    This function accepts two dictionaries of domains and their counts and compares them.

    :param count1: dict of domains and counts for protein 1
    :param count2: dict of domains and counts for protein 2
    :param prot1: protein 1 ID
    :param prot2: protein 2 ID
    ============================================================================================="""

    # Count number of matching domains
    pairs = {}
    matches = 0
    for dom in count1:
        #if dom.startswith('CL'):
            #continue
        if dom in count2:

            # Add number of matches as key to dict, value is list of matching prots and domains
            matches = min(count1[dom], count2[dom])

            # Sum number of total domains in each count dict
            total1, total2 = sum(count1.values()), sum(count2.values())
            if matches < max(total1, total2):  # Don't include pair if all domains match
                pairs[matches] = pairs.get(matches, []) + [prot1, prot2, dom]

    return pairs


def compare_domains(prots: dict) -> dict:
    """=============================================================================================
    This function accepts a dict of prot ID's and domains and returns all pairs that share at
    least one domain, but not all domains.

    :param prots: dict of prot ID's and domains
    :return dict: key is two prot ID's, value is a list of matching domains and number of matches
    ============================================================================================="""

    pairs = {}
    for prot1, doms1 in prots.items():

        # Keep count of each domain in this protein
        count1 = {}
        for dom in doms1:
            count1[dom] = count1.get(dom, 0) + 1

        # Get domain counts for all other proteins
        for prot2, doms2 in prots.items():
            if f'{prot2};{prot1}' in pairs:  # Ignore repeat pairs
                continue
            count2 = {}
            for dom in doms2:
                count2[dom] = count2.get(dom, 0) + 1

            # Compare the two dictionaries
            if prot1 != prot2 and count1 != count2:
                compare = compare_dicts(count1, count2, prot1, prot2)

                # Add to pairs dict
                for match, pair in compare.items():
                    key = f'{pair[0]};{pair[1]}'
                    value = pair[2]
                    pairs[key] = pairs.get(key, []) + [match, value]

    return pairs


def main():
    """=============================================================================================
    Main initializes file path, reads the nomax50 file for proteins and their domains with
    read_file(), and then returns the proteins which share at least one domain, but not all domains
    with compare_domains(). The resulting dict of pairs is saved as a pickle file.
    ============================================================================================="""

    file = 'parsing/pfam_nomax50.info'
    prots = read_file(file)
    pairs = compare_domains(prots)

    # Save to pickle
    with open('parsing/pfam_nomax50_pairs.pkl', 'wb') as file:
        pickle.dump(pairs, file)


if __name__ == '__main__':
    main()
