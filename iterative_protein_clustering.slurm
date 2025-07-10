#!/bin/bash

#SBATCH --partition=kingspeak
#SBATCH --account=jacksonr
# #SBATCH --qos=lonepeak
#SBATCH --time=24:00:00         # Walltime limit (HH:MM:SS)
#SBATCH --mem=32000            # Memory in MB (32 GB)
#SBATCH --nodes=1              # One node
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH --ntasks=1             # One task (single process)
#SBATCH --cpus-per-task=16     # 24 CPUs (threads) for the task
#SBATCH --job-name="itercust"  # Job name
# #SBATCH --mail-user=A02353507@usu.edu
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

cluster=/uufs/chpc.utah.edu/common/home/u6054205/JacksonLab/scripts-git/clustering_proteins/clustering.py

module purge
module load python/3.12.4
module load mmseqs2/oct24
module load muscle/5.3
module load perl/5.36.0

export MMSEQS_NUM_THREADS=$SLURM_CPUS_PER_TASK

source ~/myenvs/proti_clust_env/bin/activate
# pip install networkx

# Run clustering.py with arguments
python3 "$cluster" \
  --input-file "${input_file}" \
  --output-dir "${output_dir}" \
  --result-name "${result_name}" \
  --cluster-sensitivity 7.5 \
  --search-sensitivity 7.5 \
  --max-iter-strict 10

#deactivate
