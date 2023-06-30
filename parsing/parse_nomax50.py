"""================================================================================================
This script takes a file with protein sequence ID's and their pfam domains. It returns the ID's
which share at least one domain, but not all domains.

Ben Iovino  06/28/23   DCTDomain
================================================================================================"""

import argparse
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


def compare_dicts(count1: dict, count2: dict, args: argparse.Namespace) -> list:
    """=============================================================================================
    This function accepts two dictionaries of domains and their counts and compares them.

    :param count1: dict of domains and counts for protein 1
    :param count2: dict of domains and counts for protein 2
    :param args: argparse namespace
    :return dict: list of all matching domains (repeats are included)
    ============================================================================================="""

    # Ignore if all domains match
    if count1 == count2:
        return

    # Get all matching domains between the two proteins
    matches = []
    for dom1, co1 in count1.items():
        if args.c == 'no' and dom1.startswith('CL'):
            continue
        for dom2, co2 in count2.items():
            if dom1 == dom2:
                matches += [dom1] * min(co1, co2)

    return matches


def number_domains(prots: dict, args: argparse.Namespace) -> dict:
    """=============================================================================================
    This function accepts a dict of prot ID's and domains and returns all pairs that share at
    least one domain, but not all domains, in any order

    :param prots: dict of prot ID's and domains
    :param args: argparse namespace
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

            # Ignore repeat pairs
            if f'{prot2};{prot1}' in pairs:
                continue
            count2 = {}
            for dom in doms2:
                count2[dom] = count2.get(dom, 0) + 1

            # Compare the two dictionaries
            matches = compare_dicts(count1, count2, args)

            if matches:
                ids = f'{prot1};{prot2}'
                pairs[ids] = pairs.get(ids, [0]) + matches
                pairs[ids][0] += len(matches)

    return pairs


def match_domains(prots: dict, args = argparse.Namespace) -> dict:
    """=============================================================================================
    This function accepts a dict of prot ID's and domains and returns all pairs that share at
    least one domain, but not all domains, in the same order.

    :param prots: dict of prot ID's and domains
    :param args: argparse namespace
    :return dict: key is two prot ID's, value is a list of matching domains and number of matches
    ============================================================================================="""

    pairs = {}
    for prot1, doms1 in prots.items():
        for prot2, doms2, in prots.items():
            if prot1 == prot2:
                continue

            if f'{prot2};{prot1}' in pairs:  # Ignore repeat pairs
                continue

            # Find matching domains in order
            matches = []
            if doms1 == doms2:
                continue
            for i in range(min(len(doms1), len(doms2))):
                if doms1[i] == doms2[i]:
                    if args.c == 'no' and doms1[i].startswith('CL'):
                        continue
                    matches.append(doms1[i])

            # Add matches to dict if not all domains match
            if matches and len(matches) < max(len(doms1), len(doms2)):
                key = f'{prot1};{prot2}'
                value = matches
                pairs[key] = pairs.get(key, []) + value

    return pairs


def main():
    """=============================================================================================
    Main initializes file path, reads the nomax50 file for proteins and their domains with
    read_file(), and then returns the proteins which share at least one domain, but not all domains
    with compare_domains(). The resulting dict of pairs is saved as a pickle file.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, help='match/nomatch locations', default='nomatch')
    parser.add_argument('-c', type=str, help='include clans or not', default='no')
    args = parser.parse_args()

    # Get proteins and domains from file
    file = 'parsing/pfam_nomax50.info'
    prots = read_file(file)

    # Get pairs of proteins
    if args.m == 'match':
        pairs = match_domains(prots, args)
    if args.m == 'nomatch':
        pairs = number_domains(prots, args)

    # Save to pickle
    with open(f'parsing/pfam_nomax50_{args.m}_pairs.pkl', 'wb') as file:
        pickle.dump(pairs, file)


if __name__ == '__main__':
    main()
