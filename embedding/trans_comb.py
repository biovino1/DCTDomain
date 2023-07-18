"""================================================================================================
This script takes DCT transforms from npy files and combines them into a single numpy array.

Ben Iovino  07/14/23   DCTDomain
================================================================================================"""

import argparse
import numpy as np


def load_transforms(layers: list) -> dict:
    """=============================================================================================
    This function takes a list of layers used to embed sequences and loads the corresponding
    transforms.

    :param layers: list of layers to load
    :return: dictionary of transforms with sequence as key and list of transforms as value
    ============================================================================================="""

    trdict = {}
    for layer in layers:
        with open(f'max50_data/embeddings_{layer}/transforms.npy', 'rb') as emb:
            transforms = np.load(emb, allow_pickle=True)

            # Add sequence as dict key with list as value if not exist
            for seq, transform in transforms:
                trdict[seq] = trdict.get(seq, []) + [transform]

    return trdict


def comb_transforms(trdict: dict):
    """=============================================================================================
    This function takes a dictionary of transforms and combines them into a single transform and
    then saves them as their original data structure.
    ============================================================================================="""

    # Concatenate all transforms for each sequence
    for key, val in trdict.items():
        trdict[key] = np.concatenate(val)

    # Save as original data structure
    transforms = []
    for seq, transform in trdict.items():
        transforms.append(np.array([seq, transform], dtype=object))
    with open('max50_data/transforms.npy', 'wb') as emb:
        np.save(emb, transforms)


def main():
    """=============================================================================================
    Main
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', type=list, default=[17, 23], help='list of layers to combine')
    args = parser.parse_args()

    # Load each transform, concatenate, and save
    trdict = load_transforms(args.l)
    comb_transforms(trdict)



if __name__ == '__main__':
    main()
