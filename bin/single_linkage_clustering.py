#!/usr/bin/env python

## Agglomerative clustering (equvalent of `usearch -cluster_aggd -linkage min`)

## Example usage:
# single_linkage_clustering.py \
#  --input_file  Cluster0_dist.txt.gz \
#  --output_file ./calc_distm_out/Cluster0_out_03 \
#  --cutoff 0.97

## Input:
# Tab-delimited (square) distance matrix (optionally, gzip-compressed),
# with header containg sequence labels

## Output:
# Tab-delimted text file (without header) with one line per sequence,
# and two columns: cluster number (zero-based) and sequence label

## TO DO:
# - use condenced form of dissim matrix? (but how to keep track the lables?)

import fastcluster
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster
import pandas as pd
import numpy as np
import gzip
import argparse

def load_dissimilarity_matrix(filename):
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            df = pd.read_csv(f, sep='\t', header=0, index_col=None, dtype=float)
    else:
        df = pd.read_csv(filename, sep='\t', header=0, index_col=None, dtype=float)

    print("DataFrame shape:", df.shape)
    seq_names = df.columns.tolist()
    matrix = df.values
    return matrix, seq_names

def single_linkage_clustering(matrix):
    assert matrix.shape[0] == matrix.shape[1], "Matrix must be square."
    condensed_matrix = squareform(matrix, checks=False)
    linkage_matrix = fastcluster.linkage(condensed_matrix, method="single")
    return linkage_matrix

def assign_cluster_labels(linkage_matrix, cutoff):
    labels = fcluster(linkage_matrix, cutoff, criterion='distance')
    labels_zero_based = labels - 1
    return labels_zero_based

def save_clusters(filename, labels, seq_names):
    df = pd.DataFrame({'Cluster': labels, 'SeqID': seq_names})
    df_sorted = df.sort_values(by=['Cluster', 'SeqID'])
    df_sorted.to_csv(filename, sep='\t', index=False, header=False)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform single linkage clustering on a tab-separated distance matrix with header.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file containing the dissimilarity matrix.')
    parser.add_argument('--output_file', type=str, required=True, help='Output file to save the clustering results.')
    parser.add_argument('--cutoff', type=float, required=True, help='Cutoff value for clustering.')
    return parser.parse_args()

def main(args):
    matrix, seq_names = load_dissimilarity_matrix(args.input_file)
    linkage_matrix = single_linkage_clustering(matrix)
    labels = assign_cluster_labels(linkage_matrix, args.cutoff)
    save_clusters(args.output_file, labels, seq_names)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
