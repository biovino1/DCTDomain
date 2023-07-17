"""================================================================================================
This script takes embeddings from numpy arrays and transforms them with DCT.

Ben Iovino  07/14/23   DCTDomain
================================================================================================"""

import argparse
import os
import numpy as np
from scipy.fft import dct, idct


def scale(vec: np.ndarray) -> np.ndarray:
    """=============================================================================================
    Scale from protsttools. Takes a vector and returns it scaled between 0 and 1.

    :param v: vector to be scaled
    :return: scaled vector
    ============================================================================================="""

    maxi = np.max(vec)
    mini = np.min(vec)
    return (vec - mini) / float(maxi - mini)


def iDCTquant(vec: np.ndarray, num: int) -> np.ndarray:
    """=============================================================================================
    iDCTquant from protsttools. Takes a vector and returns its inverse discrete cosine transform.

    :param v: vector to be transformed
    :param n: number of coefficients to keep
    :return: transformed vector
    ============================================================================================="""

    f = dct(vec.T, type=2, norm='ortho')
    trans = idct(f[:,:num], type=2, norm='ortho')  #pylint: disable=E1126
    for i in range(len(trans)):  #pylint: disable=C0200
        trans[i] = scale(trans[i])  #pylint: disable=E1137
    return trans.T  #pylint: disable=E1101


def quant2D(emb: np.ndarray, n: int, m: int) -> np.ndarray:
    """=============================================================================================
    quant2D from protsttools. Takes an embedding and returns its inverse discrete cosine transform
    on both axes.

    :param emb: embedding to be transformed (n x m array)
    :param n: number of coefficients to keep on first axis
    :param m: number of coefficients to keep on second axis
    :return: transformed embedding (n*m 1D array)
    ============================================================================================="""

    dct = iDCTquant(emb[1:len(emb)-1],n)  #pylint: disable=W0621
    ddct = iDCTquant(dct.T,m).T
    ddct = ddct.reshape(n*m)
    return (ddct*127).astype('int8')


def transform_embeds(embeds: list, args: argparse.Namespace):
    """=============================================================================================
    This function takes a list of embedding files and transforms them with DCT, saving the entire
    collection to a single file.

    :param embeds: list of embedding files
    :param args: DCT dimensions
    ============================================================================================="""

    dcts = []
    for embed in embeds:
        embed = np.load(embed, allow_pickle=True)
        embed[1] = quant2D(embed[1], args.s1, args.s2)
        dct_vec = np.array([embed[0], embed[1]], dtype=object)
        dcts.append(dct_vec)

    with open('nomax_data/transforms.npy', 'wb') as emb:
        np.save(emb, dcts)


def main():
    """=============================================================================================
    Main takes a directory of embeddings and transforms them with DCT. The resulting transforms are
    saved to a single file.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, default='nomax_data/embeddings')
    parser.add_argument('-s1', type=int, default=5)
    parser.add_argument('-s2', type=int, default=44)
    args = parser.parse_args()

    # Make directory for transforms
    if not os.path.exists('nomax_data/transforms'):
        os.mkdir('nomax_data/transforms')

    # Get list of all embedding files and transform all of them
    embeds = [f'{args.d}/{file}' for file in os.listdir(args.d)]
    transform_embeds(embeds, args)


if __name__ == '__main__':
    main()
