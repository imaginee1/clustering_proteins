#!/bin/bash

#SBATCH --partition=notchpeak-guest          # Use the full node (not shared)
#SBATCH --account=owner-guest                # Your general Notchpeak account
#SBATCH --qos=notchpeak-guest                     # QoS for full node usage
#SBATCH -C rom                              # Request AMD Rome architecture
#SBATCH --time=24:00:00                     # Max allowed walltime (3 days)
#SBATCH --mem=32000                        # Full 256 GB memory in MB
#SBATCH --nodes=1                           # One physical node
#SBATCH --ntasks=1                          # One task (job step)
#SBATCH --cpus-per-task=64                  # Use all 64 physical cores
#SBATCH -o slurm-%j.out-%N                  # Output log file
#SBATCH --job-name="itercust"              # Job name
# #SBATCH --mail-user=YOUR_EMAIL@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# If you use Intel MKL, set this to optimize on AMD CPUs:
export MKL_DEBUG_CPU_TYPE=5



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
