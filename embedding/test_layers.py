"""================================================================================================
This script embeds all protein sequences in pfam_max50.fasta using a particular ESM-2 layer,
transforms them using DCT, and then finds the AUROC of the distance between all pairs of proteins
in pfam_max50.pair.

Ben Iovino  07/17/23   DCTDomain
================================================================================================"""

import logging
import os

logging.basicConfig(filename='embedding/test_layers.log',
                     level=logging.INFO, format='%(message)s')


def main():
    """=============================================================================================
    Main
    ============================================================================================="""

    # Embed, transform, and compute distances
    for i in range(1, 37):
        logging.info('Embedding sequences (layer %s)...', i)
        os.system(f'python embedding/embed.py -f embedding/pfam_max50.fasta -l {i}')
        logging.info('Transforming embeddings (layer %s)...', i)
        os.system('python embedding/transform.py')
        logging.info('Computing distances (layer %s)...\n', i)
        os.system(f'python embedding/distance.py -l {i}')


if __name__ == '__main__':
    main()