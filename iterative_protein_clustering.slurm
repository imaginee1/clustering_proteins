#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=lonepeak  
#SBATCH --account=jacksonr
#SBATCH --qos=lonepeak
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=64000 # memory given in MB
#SBATCH --nodes=1   # number of nodes
# #SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="itercust"
# #SBATCH --mail-user=A02353507@usu.edu   # insert email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# load modules
module load mmseqs2/oct24
module load muscle/5.3
module load perl/5.36.0

# Load your environment manager
module load miniconda3/23.11.0

# Activate the environment
source activate proti_clust

# Run the script
python3 /path/to/iterative_clustering.py \
  --input-file "$INPUT_FASTA" \
  --output-dir "$OUTDIR" \
  --result-name "cluster_output"
