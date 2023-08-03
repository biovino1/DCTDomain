#!/bin/bash
# this script is used to setup conda environment for all benchmarking tools
# benchmarks.py actually calls each tool for each pair in pfam_max50 dataset

# channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# setting up environment
conda create -n benchmarking
source activate base
conda activate benchmarking
conda install python=3.10.6
python -m pip install Bio
python -m pip install scikit-learn
conda install -c bioconda hhsuite  # lots of incompatible dependencies between hhsuite and latest version of blast (2.14.0)
conda install -c bioconda blast
conda install -c bioconda blast-legacy
conda install -c bioconda csblast
conda install -c bioconda fasta3
conda install -c bioconda hmmer

# manually install usearch and add to path
wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
mv usearch11.0.667_i86linux32.gz $CONDA_PREFIX/bin/usearch.gz
cd $CONDA_PREIX
gunzip bin/usearch.gz
mkdir -p ./etc/conda/activate.d
touch ./etc/conda/activate.d/env_vars.sh
echo '#!'"/bin/sh" > ./etc/conda/activate.d/env_vars.sh
echo -en '\n' >> ./etc/conda/activate.d/env_vars.sh
echo 'export PATH=$PATH:$CONDA_PREFIX/bin/usearch' >> ./etc/conda/activate.d/env_vars.sh
chmod -R 777 bin/usearch

### TESTING ###

## CSBLAST COMMANDS
# formatdb -t test -i test_db.seq - p T
# csblast -i test_query.seq -d test_db.seq -D /home/ben/anaconda3/data/K4000.lib --blast-path /home/ben/anaconda3/envs/benchmarking/bin
