"""================================================================================================
This script embeds sequences from pfam_max50.fasta and transforms them with DCT dimensions
of 5x44. The distance between the vectors are compared to each other and the AUROC is computed
using pfam_max50.pair as the ground truth. AUROC is written in distance.log and used to determine
the best layer of ESM-2 to use for embedding the protein sequences.

Ben Iovino  07/18/23   DCTDomain
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
    for i in range(1, 36):
        logging.info('Embedding sequences (layer %s)...', i)
        os.system(f'python embedding/embed.py -f embedding/pfam_max50.fasta -l {i}')
        logging.info('Transforming embeddings (layer %s)...', i)
        os.system(f'python embedding/transform.py -l {i}')
        logging.info('Computing distances (layer %s)...\n', i)
        os.system(f'python embedding/distance.py -l {i}')
        os.system(f'rm -rf max50_data/embeddings_{i}')


if __name__ == '__main__':
    main()
