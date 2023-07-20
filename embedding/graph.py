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
    plt.savefig('embedding/graphs/auroc.png')


def heat_map(auroc: list):
    """=============================================================================================
    This function graphs the AUROC values from one layer of ESM-2 as a heat map.

    :param auroc: list of AUROC values
    ============================================================================================="""

    # dct dimensions
    layer = auroc[0][0]
    i = [3, 4, 5, 6, 7, 8]
    j = [20, 30, 40, 50, 60, 70, 80]

    # Create a dictionary of AUROC values
    auroc_dict = {}
    for s1 in i:
        for s2 in j:
            auroc_dict[(s1, s2)] = auroc.pop(0)[1]

    # Create a heat map
    plt.figure(figsize=(12, 5))
    plt.imshow([[auroc_dict[(s1, s2)] for s2 in j] for s1 in i], cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.xlabel('col')
    plt.ylabel('rows')
    plt.title('AUROC Heat Map')
    plt.xticks(range(len(j)), j)
    plt.yticks(range(len(i)), i)
    plt.savefig(f'embedding/graphs/auroc{layer}_heat.png')


def main():
    """=============================================================================================
    Main takes a file containing AUROC values and graphs them.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', type=str, default='embedding/graphs/auroc.txt')
    parser.add_argument('-f2', type=str, default='embedding/graphs/auroc17_dim.txt')
    parser.add_argument('-f3', type=str, default='embedding/graphs/auroc23_dim.txt')
    args = parser.parse_args()

    # Get auroc values and graph
    auroc = read_auroc(args.f1)
    auroc17 = read_auroc(args.f2)
    auroc23 = read_auroc(args.f3)
    graph_auroc(auroc)
    heat_map(auroc17)
    heat_map(auroc23)


if __name__ == '__main__':
    main()
