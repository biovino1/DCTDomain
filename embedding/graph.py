"""================================================================================================
This script graphs the AUROC for each layer of ESM-2, found in auroc.txt which was extracted
from the output of distance.py.

Ben Iovino  07/18/23   DCTDomain
================================================================================================"""

import argparse
import matplotlib.pyplot as plt


def read_auroc(file: str) -> list:
    """=============================================================================================
    This function reads the AUROC values from a file and returns them as a list.

    :param file: file containing AUROC values
    :return auroc: list of AUROC values and the layers they correspond to
    ============================================================================================="""

    auroc = []
    with open(file, 'r', encoding='utf8') as file:
        for line in file:
            line = line.split()
            auroc.append((line[3].strip(':'), float(line[4])))
    return auroc


def graph_auroc(auroc: list):
    """=============================================================================================
    This function graphs the AUROC values for each layer of ESM-2.

    :param auroc: list of AUROC values
    ============================================================================================="""

    aurange = range(len(auroc))
    plt.figure(figsize=(12, 5))
    plt.plot(list(aurange), [auroc[i][1] for i in aurange], 'o-')
    plt.xlabel('Layer')
    plt.ylabel('AUROC')
    plt.title('AUROC vs. Layer')
    plt.xticks(list(aurange), [auroc[i][0] for i in aurange])
    plt.grid(axis='y', linestyle='--')
    plt.grid(axis='x', linestyle='--')
    plt.savefig('embedding/auroc.png')


def main():
    """=============================================================================================
    Main takes a file containing AUROC values and graphs them.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='embedding/auroc.txt')
    args = parser.parse_args()

    # Get auroc values and graph
    auroc = read_auroc(args.f)
    graph_auroc(auroc)


if __name__ == '__main__':
    main()
