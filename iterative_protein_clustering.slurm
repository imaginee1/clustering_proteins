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

# Check if arguments were provided
if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch iqtree2_infile_outdir.slurm <input_fasta_file> <output_directory> <result_name>"
    exit 1
fi

input_file=$(realpath "$1")
output_dir=$(realpath "$2")
result_name="$3"

# Load modules
module purge
module load miniconda3/23.11.0
module load muscle/5.3
module load perl/5.36.0

# Source Conda
source $EBROOTMINICONDA3/etc/profile.d/conda.sh

# Create env (if not already created)
# bash setup_envs.sh  # Uncomment only if setup_envs.sh is idempotent

# Activate environment
conda activate proti_clust  # Replace with your env name

#load mmseqs2
module load mmseqs2/oct24

# Go to working directory
cd $SLURM_SUBMIT_DIR

# Run clustering.py with arguments
python /path/to/clustering.py \
  --input-file "${input_file}" \
  --output-dir "${output_dir}" \
  --result-name "${result_name}" \
  --cluster-sensitivity 7.5 \
  --search-sensitivity 7.5 \
  --max-iter-strict 10
