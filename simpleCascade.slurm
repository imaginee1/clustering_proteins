#!/bin/bash

#SBATCH --partition=notchpeak-shared-short         
#SBATCH --account=notchpeak-shared-short              
#SBATCH --qos=notchpeak-shared-short                   
# #SBATCH -C rom                              
#SBATCH --time=8:00:00                     
#SBATCH --mem=64000                        
#SBATCH --nodes=1                           
#SBATCH --ntasks=1                          
#SBATCH --cpus-per-task=16               
#SBATCH -o slurm-%j.out-%N                  
#SBATCH --job-name="cascade"              
# #SBATCH --mail-user=YOUR_EMAIL@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Check if arguments were provided
if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch cascade.slurm <input_fasta_file> <output_directory> <result_name>"
    exit 1
fi

input_file=$(realpath "$1")
output_dir=$(realpath "$2")
result_name="$3"

cluster=/uufs/chpc.utah.edu/common/home/u6054205/JacksonLab/scripts-git/clustering_proteins/simpleCascade.py

# Load required modules
module purge
module load python/3.12.4
module load mmseqs2/oct24
module load muscle/5.3

# Set MMseqs to use all available cores (if --cpus-per-task is used)
export MMSEQS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export $SLURM_CPUS_PER_TASK

# Activate your Python virtual environment
source ~/myenvs/proti_clust_env/bin/activate

# Set and create SCRATCH directory
export SCRDIR="/scratch/general/vast/$USER/$SLURM_JOB_ID"
mkdir -p "$SCRDIR"

# Make sure Python sees it as $SCRATCH
export SCRATCH="$SCRDIR"

# Optionally set TMPDIR (for any other temp tools)
export TMPDIR="$SCRDIR/tmp"
mkdir -p "$TMPDIR"

echo "Using scratch directory: $SCRDIR"
echo "Using TMPDIR: $TMPDIR"

# Run clustering pipeline
python3 "$cluster" \
  --input-file "${input_file}" \
  --output-dir "${output_dir}" \
  --result-name "${result_name}" \
  --cascade-sensitivity 7.5 \
  --search-sensitivity 7.5 \
  --cascade-eval 1e-6 \
  --search-eval 1e-6

# Optional: clean up after the run
# rm -rf "$SCRDIR"

# deactivate  # optional in batch jobs
