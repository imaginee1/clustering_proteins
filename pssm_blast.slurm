#!/bin/bash

#SBATCH --partition=notchpeak-shared-short          # Use the full node (not shared)
#SBATCH --account=notchpeak-shared-short                # Your general Notchpeak account
#SBATCH --qos=notchpeak-shared-short                   # QoS for full node usage
# #SBATCH -C rom                              # Request AMD Rome architecture
#SBATCH --time=8:00:00                     # Max allowed walltime (3 days)
#SBATCH --mem=64000                        # Full 256 GB memory in MB
#SBATCH --nodes=1                           # One physical node
#SBATCH --ntasks=1                          # One task (job step)
# #SBATCH --cpus-per-task=64                  # Use all 64 physical cores
#SBATCH -o slurm-%j.out-%N                  # Output log file
#SBATCH --job-name="pssm"              # Job name
# #SBATCH --mail-user=YOUR_EMAIL@usu.edu
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

module load blast/2.11.0

psiblast -in_pssm cas12a2_3.pssm -query NiCas12a3.fasta -db nr \
  -out results.out -evalue 1e-5 -num_threads 8 -outfmt 6




