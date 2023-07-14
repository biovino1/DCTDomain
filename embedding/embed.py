"""================================================================================================
This script takes sequences from a fasta file and returns their embeddings.

Ben Iovino  07/14/23   DCTDomain
================================================================================================"""

import argparse
import esm
import numpy as np
import torch
from Bio import SeqIO


def load_seqs(file: str) -> list:
    """=============================================================================================
    This function takes a fasta file and returns a list of sequences and their IDs.

    :param file: fasta file
    :return list: list of sequences
    ============================================================================================="""

    # Read each line and add to list
    seqs = []
    with open(file, 'r', encoding='utf8') as f:
        for seq in SeqIO.parse(f, 'fasta'):
            seqs.append((seq.id, str(seq.seq)))

    return seqs


def embed_seqs(seqs: list, args: argparse.Namespace):
    """=============================================================================================
    This function takes a filename and a list of sequences. It saves a numpy array of the
    embeddings.

    :param seqs: list of sequences
    :param args: args for filename and encoder layer
    ============================================================================================="""

    # Load ESM-2 model
    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    # Load to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  # pylint: disable=E1101
    model = model.to(device)

    # Embed sequences
    batch_labels, batch_strs, batch_tokens = batch_converter(seqs)  #pylint: disable=W0612
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[args.l])
    token_representations = results["representations"][args.l]

    # Make an array of each label and its embedding, add to list
    embeds = []
    for i in range(len(seqs)):  #pylint: disable=C0200
        embeds.append(np.array([seqs[i][0], token_representations[i].numpy()], dtype=object))

    # Save embeddings
    with open(f'{args.f}.npy', 'wb') as emb:
        np.save(emb, embeds)


def main():
    """=============================================================================================
    Main
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='embedding/pfam_max50.fasta')
    parser.add_argument('l', type=int, default=36)
    args = parser.parse_args()

    # Load sequences and embed
    seqs = load_seqs(args.f)
    seqs = seqs[:19]
    embed_seqs(seqs, args)


if __name__ == '__main__':
    main()
