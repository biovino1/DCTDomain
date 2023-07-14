"""================================================================================================
This script takes embeddings from a numpy array and transforms them with DCT.

Ben Iovino  07/14/23   DCTDomain
================================================================================================"""

import argparse
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


def transform_embeds(args: argparse.Namespace):
    """=============================================================================================
    This function takes an array of embeddings and transforms them using DCT.

    :param args: file with embeddings and DCT dimensions
    ============================================================================================="""

    embeds = np.load(args.f, allow_pickle=True)
    for i, embed in enumerate(embeds):  # transform each embedding in array
        embeds[i][1] = quant2D(embed[1], args.s1, args.s2)
    with open(f'{args.f.split(".")[0]}.dct', 'wb') as emb:
        np.save(emb, embeds)


def main():
    """=============================================================================================
    Main
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='embedding/pfam_max50.emb')
    parser.add_argument('-s1', type=int, default=5)
    parser.add_argument('-s2', type=int, default=44)
    args = parser.parse_args()

    transform_embeds(args)


if __name__ == '__main__':
    main()
