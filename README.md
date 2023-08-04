**************************************************************************************************************
# Using DCT to Compress Protein Embeddings
**************************************************************************************************************

Miscellaneous scripts for benchmarking, clustering, embedding, and parsing different DCT transformed datasets.

**************************************************************************************************************
# Benchmarking
**************************************************************************************************************

A replication of the results from https://academic.oup.com/bioinformatics/article/32/17/2636/2450749 where
multiple state-of-the-art protein database searching programs are ran on the pfam-max50 dataset. conda.sh
sets up the conda environment with all of the necessary dependancies. The benchmarking script is then run
to compare each protein pair in the dataset to each other to get a bit score and E-value. The results are
then parsed and the AUROC is calculated.

Instructions on how to use each tool for the purposes of this project can be found in guide.txt in the
benchmarking directory.

**************************************************************************************************************
# Clustering
**************************************************************************************************************

Attemping to find similarities between the clustering of the DCT vectors of pfam and scop domains and the
grouping of the sequences as defined by the pfam and scop databases. Several different clustering algorithms
are used.

Results:

The best results of clustering pfam dct vectors was given by average linkage hierarchical clustering with
a manhattan distance metric with a precision and recall of 0.559 and 0.269, respectively. The number of
clusters was set to 2000, which is much lower than the number of pfam families (19,632).

The best results of clustering scop dct vectors was given by single linkage hierarchical clustering with
a manhatten distance metric with a precision and recall of 0.892 and 0.762, respectively. The number of
of clusters was set to 9500, which is much higher than the number of scop families (5936).

While the scop dct clustering performed much better than pfam dct clustering, the number of clusters does
not match the actual number of families so it is unclear what the clusters are actually representing.

**************************************************************************************************************
# Embedding
**************************************************************************************************************

The pfam-max50 dataset is embedded using various layers of ESM2 and transformed with different DCT parameters.
The manhattan distance (transformed into a similarity metric) is then calculated between each pair of
transformed embeddings and the AUROC is calculated from the actual homology classifications to determine which
layer and DCT parameters are best for determining homology. 

Results:

For ESM2_t33_650M, the best AUROC was given by layers 16 and 22 with DCT dimensions of 3x80 and 5x60,
respectively. The AUROC for each layer can be found in in embedding/graphs/auroc.png, and heat maps for the DCT
parameters can be found in the same directory.

**************************************************************************************************************
# Parsing
**************************************************************************************************************

The FUpred predictions for the pfam-max50 dataset are parsed to examine the predictions.