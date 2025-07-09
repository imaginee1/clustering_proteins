#!/bin/bash
module purge
module load python/3.10.4

python -m venv ~/myenvs/proti_clust_env
source ~/myenvs/proti_clust_env/bin/activate
pip install --upgrade pip
pip install networkx
