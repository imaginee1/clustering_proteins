#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=notchpeak-shared-short 
#SBATCH --account=notchpeak-shared-short
#SBATCH --qos=notchpeak-shared-short
#SBATCH --time=8:00:00   # walltime limit (HH:MM:SS)
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

cluster=/uufs/chpc.utah.edu/common/home/u6054205/JacksonLab/scripts-git/clustering_proteins/clustering.py

module purge
module load python/3.12.4
module load mmseqs2/oct24
module load muscle/5.3
module load perl/5.36.0

source ~/myenvs/proti_clust_env/bin/activate

# Run clustering.py with arguments
python3 "$cluster" \
  --input-file "${input_file}" \
  --output-dir "${output_dir}" \
  --result-name "${result_name}" \
  --cluster-sensitivity 7.5 \
  --search-sensitivity 7.5 \
  --max-iter-strict 10
