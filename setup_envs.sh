#!/bin/bash

# ======== Settings ========
CONDA_ENV_NAME="proti_clust"
PYTHON_VERSION="3.10"

# ======== Load Miniconda module (adjust to your system) ========
module load miniconda3/23.11.0

# ======== Set up Conda environment ========
echo "Creating Conda environment: $CONDA_ENV_NAME"
conda create -y -n "$CONDA_ENV_NAME" python=$PYTHON_VERSION networkx
echo "Conda environment '$CONDA_ENV_NAME' created."

# Activate Conda environment and optionally install more packages
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_NAME"
# conda install -y biopython   # Add if needed
conda deactivate

echo "Conda environment setup complete."
echo "To use it: conda activate $CONDA_ENV_NAME"
