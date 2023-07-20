"""================================================================================================
Test script embeds sequences from pfam_max50.fasta using layers 17 and 23 from ESM-2 and transforms
them with DCT using various dimensions. The distance between the vectors are compared to each other
and the AUROC is computed using pfam_max50.pair as the ground truth. AUROC is written in
distance.log and used to determine the best DCT dimensions to use for reducing the representation
size of the protein embeddings.

Ben Iovino  07/19/23   DCTDomain
================================================================================================"""

import logging
import os

logging.basicConfig(filename='embedding/test_transforms.log',
                     level=logging.INFO, format='%(message)s')


def main():
    """=============================================================================================
    Main
    ============================================================================================="""

    # DCT transforms to test
    i = [3, 4, 5, 6, 7, 8]
    j = [20, 30, 40, 50, 60, 70, 80]

    # Embed using layers 17 and 23
    for lay in [17, 23]:
        logging.info('Embedding sequences (layer %s)...', i)
        os.system(f'python embedding/embed.py -f embedding/pfam_max50.fasta -l {i}')

        # For every combination of i and j, transforms embeddings from layers 17 and 23 and
        # compute distances to determine best dimensions
        for s1 in i:
            for s2 in j:

                # If there exist transforms in the directory, delete them
                if os.path.exists(f'max50_data/embeddings_{lay}/transforms.npy'):
                    os.remove(f'max50_data/embeddings_{lay}/transforms.npy')

                logging.info('Transforming embeddings (layer %s) with dims %s x %s...', lay, s1, s2)
                os.system(f'python embedding/transform.py -l {lay} -s1 {s1} -s2 {s2}')

                logging.info('Computing distances (layer %s)...\n', lay)
                os.system(f'python embedding/distance.py -l {lay}')


if __name__ == '__main__':
    main()
