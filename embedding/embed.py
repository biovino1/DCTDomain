"""================================================================================================
This script takes sequences from a fasta file and embeds each one individiually.

Ben Iovino  07/14/23   DCTDomain
================================================================================================"""

import argparse
import os
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


def embed_seq(seq: tuple, args: argparse.Namespace, model: torch.nn.Module,
               batch_converter: esm.data.BatchConverter, device: str):
    """=============================================================================================
    This function takes a filename and a protein seq and its id. It saves a numpy array of the
    embedding.

    :param seqs: list of sequences
    :param args: args for filename and encoder layer
    ============================================================================================="""

    # Embed sequences
    batch_labels, batch_strs, batch_tokens = batch_converter([seq])  #pylint: disable=W0612
    batch_tokens = batch_tokens.to(device)  # send tokens to gpu

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[args.l])
    token_representations = results["representations"][args.l]

    # Make an array of label and its embedding, save to file
    embed = np.array([seq[0], token_representations[0].cpu().numpy()], dtype=object)
    with open(f'max50_data/embeddings/{batch_labels[0]}.npy', 'wb') as emb:
        np.save(emb, embed)


def main():
    """=============================================================================================
    Main loads sequences from input file and embeds them with ESM-2. Embeddings are saved as
    individual numpy arrays because they are too large to save as a single file with the model
    loaded on the GPU (lame).
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='embedding/pfam_max50.fasta')
    parser.add_argument('-l', type=int, default=35)
    args = parser.parse_args()

    # Load sequences and model
    seqs = load_seqs(args.f)
    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    # Load to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  # pylint: disable=E1101
    model = model.to(device)

    # Make directory for embeddings if it doesn't exist
    if not os.path.exists('max50_data'):
        os.mkdir('max50_data')
        os.mkdir('max50_data/embeddings')

    # Embed each sequence
    for seq in seqs:
        embed_seq(seq, args, model, batch_converter, device)


if __name__ == '__main__':
    main()
