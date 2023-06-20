"""================================================================================================
This script takes files with domains and returns how many exact and close matches there are between
the two files.

Ben Iovino  05/24/23   DCTDomain
================================================================================================"""


def read_file(file: str) -> dict:
    """=============================================================================================
    This function accepts a file path and parses each line. The first element is the sequence, the
    second is the number of domains, and the third is the domain positions.

    :param file: file path
    :return dict: sequence is key, value is a list containing number of domains and their positions.
    ============================================================================================="""

    # Read each line and add to dict
    domains = {}
    with open(file, 'r', encoding='utf8') as file:
        for line in file:
            line = line.split()
            if line[1] != '0':  # Check if there are domains
                dom_list = line[2].split(';')[:-1]  # Split domains, last element is empty
                domains[line[0]] = [line[1], dom_list]
            else:
                domains[line[0]] = [line[1]]

    return domains


def mult_domain(dom1: str, dom2: str) -> int:
    """=============================================================================================
    This function accepts two strings, one or both of which contain multiple domains. It returns
    how many exact matches between domains there are, as well as how many are close, meaning that
    the beginning and end of every domain is within 5 positions of each other.

    :param dom1: domain string from first file
    :param dom2: domain string from second file
    return int: number of exact matches, number of close matches
    ============================================================================================="""

    # If only one domain is multiple, skip
    match_dom, close_dom = 0, 0
    split_d1, split_d2 = dom1.split(','), dom2.split(',')
    if len(split_d1) != len(split_d2):
        return match_dom, close_dom

    # Split domains and check if they are the same
    if split_d1 == split_d2:
        match_dom += 1

    # Check if BOTH beginning and end positions are both within 5 positions of each other
    else:
        match = len(split_d1)  # Number of domains that have to match
        for i, dom in enumerate(split_d1):

            # Check if corresponding domain in second string is within 5 positions
            d1_start, d1_end = int(dom.split('-')[0]), int(dom.split('-')[1])
            d2_start, d2_end = int(split_d2[i].split('-')[0]), int(split_d2[i].split('-')[1])
            if abs(d1_start - d2_start) <= 5 and abs(d1_end - d2_end) <= 5:
                match -= 1

        # If all domains match, add to close_dom
        if match == 0:
            close_dom += 1

    return match_dom, close_dom


def check_range(dom1: list, dom2: list) -> int:
    """=============================================================================================
    This function accepts two lists of strings and returns how many exact matches between strings
    there are between the lists, as well as how many are close.

    :param dom1: list of domains from first file
    :param dom2: list of domains from second file
    return int: number of exact matches, number of close matches
    ============================================================================================="""

    # Move through first list, find matching domain in second list if it exists
    match_dom, close_dom = 0, 0
    for d1 in dom1:
        for d2 in dom2:

            # Check for domains with multiple regions
            if ',' in d1 or ',' in d2:
                match, close = mult_domain(d1, d2)
                match_dom += match
                close_dom += close
                continue

            # Check for exact match
            if d1 == d2:
                match_dom += 1
            else:

                # Get beginning and end of domains
                split_d1, split_d2 = d1.split('-'), d2.split('-')
                d1_start, d1_end = int(split_d1[0]), int(split_d1[1])
                d2_start, d2_end = int(split_d2[0]), int(split_d2[1])

                # Check if beginning and end positions are both within 5 positions of each other
                if abs(d1_start - d2_start) <= 5 and abs(d1_end - d2_end) <= 5:
                    close_dom += 1

    return match_dom, close_dom


def compare_domains(domains1: dict, domains2: dict):
    """=============================================================================================
    This function accepts two dictionaries of domains and compares them.

    :param domains1: dictionary of domains from file 1
    :param domains2: dictionary of domains from file 2
    ============================================================================================="""

    num_dom1, num_dom2 = 0, 0 # Total number of domains
    more_dom1, more_dom2 = 0, 0  # Number of sequences that have more domains in file 1 and 2
    match_dom, close_dom = 0, 0  # Number of exact and close matches

    # Move through first dict, find matching sequence in second dict if it exists
    for seq, domains in domains1.items():
        num_dom1 += int(domains[0])
        if seq in domains2:
            num_dom2 += int(domains2[seq][0])

            # Check if there are more domains in file 1 or 2
            if int(domains[0]) > int(domains2[seq][0]):
                more_dom1 += 1
            elif int(domains[0]) < int(domains2[seq][0]):
                more_dom2 += 1

            # Check how many domains in first dict match second dict
            if int(domains2[seq][0]) != 0:
                match, close = check_range(domains[1], domains2[seq][1])
                match_dom += match
                close_dom += close

    print(f'Number of domains in test.info = {num_dom1}\n')
    print(f'Number of domains in test.FUpred = {num_dom2}\n')
    print(f'Number of matching domains = {match_dom}\n')
    print(f'Number of close (within 5 pos) domains = {close_dom}\n')
    print(f'Number of sequences with more domains in test.info = {more_dom1}\n')
    print(f'Number of sequences with more domains in test.FUpred = {more_dom2}\n')


def main():
    """=============================================================================================
    Main initializes file paths, reads them both for their domains with read_file(), and then
    compares their domains with compare_domains().
    ============================================================================================="""

    # File paths
    file1 = 'pfam_data/test.info'
    file2 = 'pfam_data/test.FUpred'

    # Parse files and compare domains
    f1_domains = read_file(file1)
    f2_domains = read_file(file2)
    compare_domains(f1_domains, f2_domains)


if __name__ == '__main__':
    main()
