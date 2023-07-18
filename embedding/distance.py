"""================================================================================================
This script computes the distance between all DCT vectors in a file and saves the results.

Ben Iovino  07/17/23   DCTDomain
================================================================================================"""

import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cityblock
from sklearn import metrics

logging.basicConfig(filename='max50_data/distance.log',
                     level=logging.INFO, format='%(message)s')


def comp_dist(transforms: np.ndarray) -> tuple:
    """=============================================================================================
    This function takes an array DCT vectors and computes the distance between each one.

    :param transforms: array of DCT vectors, first index is label and second is DCT vector
    :return:
        pairs: list of tuples containing the pair of labels and their relationship (hom, nonhom)
        distances: list of tuples containing the pair of labels and the distance between them
    ============================================================================================="""

    # Load pfam_max50 pairs
    pairs = []
    with open('embedding/pfam_max50.pair', 'r', encoding='utf8') as file:
        for line in file:
            if line.startswith('#'):  # Skip first line
                continue
            pairs.append(line.split())
    prot_pairs = [(pairs[i][0], pairs[i][1]) for i in range(len(pairs))]

    # Compute distances
    distances = []
    for i in range(len(transforms)):  #pylint: disable=C0200
        for j in range(len(transforms)):  #pylint: disable=C0200
            if i == j:
                continue

            # Only if pair is in pfam_max50
            if (transforms[i][0], transforms[j][0]) in prot_pairs:
                dist = 1/cityblock(transforms[i][1], transforms[j][1])
                logging.info('%s %s %s', transforms[i][0], transforms[j][0], dist)
                distances.append((transforms[i][0], transforms[j][0], dist))

    return pairs, distances


def comp_auroc(pairs: list, distances: list):
    """=============================================================================================
    This function takes two lists of protein pairs, one with a label and the other with the distance
    between their DCT vectors. It then computes AUROC.

    :param pairs: list of tuples containing the pair of labels and their relationship (hom, nonhom)
    :param distances: list of tuples containing the pair of labels and the distance between them
    :return: list of lists containing the fpr, tpr, and auc
    ============================================================================================="""

    # Sort pairs and distances by first label in pair
    pairs.sort(key=lambda x: (x[0], x[1]))
    distances.sort(key=lambda x: (x[0], x[1]))

    # Convert to dataframes
    label = pd.DataFrame(pairs, columns=['prot1', 'prot2', 'label'])
    dist = pd.DataFrame(distances, columns=['prot1', 'prot2', 'dist'])

    # Add labels to distances and change to be 0 or 1 based on homology
    dist['label'] = label['label']
    dist['label'] = dist['label'].apply(lambda x: 0 if x == 'nonhom' else 1)

    actual, preds, rocs = dist['label'], dist['dist'], []
    fpr, tpr, thresholds = metrics.roc_curve(actual, preds, pos_label=1)  #pylint: disable=W0612
    auc = metrics.auc(fpr, tpr)
    rocs.append([fpr, tpr, auc])

    return rocs


def graph_auroc(rocs: list, layer: int):
    """=============================================================================================
    This function takes a list of lists containing the fpr and tpr for each threshold and graphs
    them.

    :param rocs: list of lists containing the fpr, tpr, and auc
    :param l: layer of ESM-2 used
    ============================================================================================="""

    fig = plt.figure()
    ax = fig.add_subplot()
    auc = rocs[0][2]
    ax.plot([0,1], [0,1], "k--", label="Random (AUC = 0.500)")
    ax.plot(rocs[0][0], rocs[0][1], label=f'AUROC: {auc:.3f}')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'AUROC (ESM2 Layer {layer}')
    ax.legend(loc='lower right')
    plt.savefig(f'max50_data/auroc{layer}.png')


def main():
    """=============================================================================================
    Main takes a file of DCT vectors and computes the distance between each one. It then computes
    AUROC and graphs the results.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='max50_data/transforms.npy')
    parser.add_argument('-l', type=int, default=35)
    args = parser.parse_args()

    # Load file and compute distances
    transforms = np.load(args.f, allow_pickle=True)
    pairs, distances = comp_dist(transforms)

    # Compute auroc and save results
    rocs = comp_auroc(pairs, distances)
    logging.info('AUROC For Layer %s: %s\n', args.l, rocs[0][2])
    graph_auroc(rocs, args.l)


if __name__ == '__main__':
    main()
