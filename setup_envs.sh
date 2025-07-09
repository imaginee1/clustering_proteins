#!/bin/bash

# ======== Settings ========
CONDA_ENV_NAME="proti_clust"
VENV_DIR="./proti_clust_env"
PYTHON_VERSION="3.10"

# ======== Load Miniconda module (adjust to your system) ========
module load miniconda3/23.11.0

# ======== Set up Conda environment ========
echo "Creating Conda environment: $CONDA_ENV_NAME"
conda create -y -n "$CONDA_ENV_NAME" python=$PYTHON_VERSION networkx
echo "Conda environment '$CONDA_ENV_NAME' created."

# Optionally activate and install additional packages
source activate "$CONDA_ENV_NAME"
# conda install -y biopython   # Add if needed
conda deactivate

# ======== Set up venv environment ========
echo "Creating venv at: $VENV_DIR"
python -m venv "$VENV_DIR"

# Activate venv
source "$VENV_DIR/bin/activate"

# Upgrade pip and install packages
pip install --upgrade pip
pip install networkx
# pip install biopython   # Add if needed

# Deactivate venv
deactivate

echo "Both Conda and venv environments created."
echo "To use Conda: source activate $CONDA_ENV_NAME"
echo "To use venv: source $VENV_DIR/bin/activate"
