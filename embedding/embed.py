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
from transformers import T5EncoderModel, T5Tokenizer


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


def prot_t5xl_embed(seq: tuple, args: argparse.Namespace,
                    tokenizer: T5Tokenizer, encoder: T5EncoderModel, device:str) -> list:
    """=============================================================================================
    This function accepts a protein sequence and returns a list of vectors, each vector representing
    a single amino acid using RostLab's ProtT5_XL_UniRef50 model.

    :param seq: protein sequence
    :param tokenizer: tokenizer model
    :param encoder: encoder model
    :param device: gpu/cpu
    :param layer: layer to extract features from
    return: list of vectors
    ============================================================================================="""

    # Tokenize, encode, and load sequence
    ids = tokenizer.batch_encode_plus(seq[1], add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids']).to(device)  # pylint: disable=E1101
    attention_mask = torch.tensor(ids['attention_mask']).to(device)  # pylint: disable=E1101

   # Foward hook for extracting features from layer n
    activation = {}
    def get_activation(name):
        def hook(model, input, output):  # pylint: disable=W0613, W0622
            output = output.numpy()
            activation[name] = output
        return hook
    encoder.encoder.block[args.l].layer[1].layer_norm.register_forward_hook(
        get_activation(f'layer_{args.l}'))

    # Extract sequence features
    with torch.no_grad():
        embedding = encoder(input_ids=input_ids,attention_mask=attention_mask)
    embedding = activation[f'layer_{args.l}']

    # Remove padding and special tokens
    features = []
    for seq_num in range(len(embedding)):  # pylint: disable=C0200
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len-1]
        features.append(seq_emd)

    # Make an array of label and its embedding, save to file
    embed = np.array([seq[0], features[0]], dtype=object)
    with open(f'max50_data/embeddings_{args.l}/{seq[0]}.npy', 'wb') as emb:
        np.save(emb, embed)


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
    with open(f'max50_data/embeddings_{args.l}/{batch_labels[0]}.npy', 'wb') as emb:
        np.save(emb, embed)


def main():
    """=============================================================================================
    Main loads sequences from input file and embeds them with ESM-2. Embeddings are saved as
    individual numpy arrays because they are too large to save as a single file with the model
    loaded on the GPU (lame).
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, default='embedding/pfam_max50.fasta')
    parser.add_argument('-l', type=int, default=23)
    args = parser.parse_args()

    # Load sequences and model
    seqs = load_seqs(args.f)
    #model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    #batch_converter = alphabet.get_batch_converter()
    #model.eval()

    # Load the tokenizer
    # Load to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  # pylint: disable=E1101
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").to(device)

    # Make directory for embeddings if it doesn't exist
    if not os.path.exists('max50_data'):
        os.mkdir('max50_data')
    if not os.path.exists('max50_data/embeddings'):
        os.mkdir(f'max50_data/embeddings_{args.l}')

    # Embed each sequence
    for seq in seqs:
        prot_t5xl_embed(seq, args, tokenizer, model, device)


if __name__ == '__main__':
    main()
