# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 14:10:23 2025

@author: imagi
"""

import argparse as ap
# import numpy as np
import subprocess as subp
from pathlib import Path
import sys
import networkx as nx
from collections import defaultdict
from typing import List, Set, Tuple
import logging
import os


parser = ap.ArgumentParser()

parser.add_argument("--input-file",
                    action="store",
                    type=str,
                    dest="input_file",
                    help="input filepath for a multi-fasta") #, 
                   # default="input")

parser.add_argument("--output-dir",
                    action="store",
                    type=str,
                    dest="output_dir",
                    help="output path for result folder", 
                    default="output")

parser.add_argument("--result-name",
                    action="store",
                    type=str,
                    dest="result_name",
                    help="name of result files", 
                    default="result-files")

# optional arguments

parser.add_argument("--min-clust-id",
                    action="store",
                    type=float,
                    dest="min_id",
                    help="decimal identity threshold for permissive clustering", 
                    default=float(0.30))

parser.add_argument("--min-clust-cov",
                    action="store",
                    type=float,
                    dest="min_cov",
                    help="decimal cov threshold for permissive clustering", 
                    default=float(0.8))

parser.add_argument("--max-clust-id",
                    action="store",
                    type=float,
                    dest="max_id",
                    help="decimal identity threshold for strict clustering", 
                    default=float(0.9))

parser.add_argument("--max-clust-cov",
                    action="store",
                    type=float,
                    dest="max_cov",
                    help="decimal coverage threshold for strict clustering", 
                    default=float(0.9))

parser.add_argument("--search-cov",
                    action="store",
                    type=float,
                    dest="search_cov",
                    help="decimal coverage threshold for strict clustering", 
                    default=float(0.7))

parser.add_argument("--cascade-iter",
                    action="store",
                    type=int,
                    dest="cascade_iter",
                    help="max number of strict clustering/search iterations till convergence", 
                    default=int(10))
"""
arser.add_argument("--final-clust-id",
                    action="store",
                    type=float,
                    dest="strict_id",
                    help="decimal identity threshold for strict clustering", 
                    default=float(0.9))

parser.add_argument("--final-clust-cov",
                    action="store",
                    type=float,
                    dest="strict_cov",
                    help="decimal coverage threshold for strict clustering", 
                    default=float(0.9))
"""

# optional optional arguments
"""
parser.add_argument("--ref-protein-fasta",
                    action="store",
                    dest="ref_fasta",
                    help="fasta alignment of the reference protein", 
                    default=None)

parser.add_argument("--align-ref",
                    action="store",
                    type=bool,
                    dest="store",
                    help="set as True if your reference protein fasta is unaligned",
                    default=False)
"""

parser.add_argument("--cascade-sensitivity",
                    action="store",
                    type=float,
                    dest="cascade_sens",
                    help="desired mmseqs sequence clustering sensitivity (1.0 ≤ s ≤ 9.5)", 
                    default=float(7.5))

parser.add_argument("--search-sensitivity",
                    action="store",
                    type=float,
                    dest="search_sens",
                    help="desired mmseqs profile-sequence search sensitivity (1.0 ≤ s ≤ 9.5)", 
                    default=float(7.5))




"""
parser.add_argument("--require-reciprocal",
                    action="store",
                    dest="search_sens",
                    help="when identifying connected components, if search hits must be reciprocal to be included", 
                    default=False)
"""

inarg = parser.parse_args()





############################################

## load modules (CHPC)
#subp.run("""
#         module load mmseqs2/oct24
#         module load muscle/5.3
#         module load perl/5.36.0
#         """, shell=True, check=True)
# add commands to make sure reformat.pl and join.awk are executable


SCRIPT_DIR = Path(__file__).resolve().parent

REFORMAT_PL = SCRIPT_DIR / "reformat.pl"


subp.run(f"chmod +x {REFORMAT_PL}", shell=True, check=True)
         
############################################



# check that input file and output directory exist; exit code if no

input_file = Path(inarg.input_file).resolve()

if input_file.is_file():
    print(f"Exists at: {input_file}")
else:
    print(f"Error: output file not found at {input_file}", file=sys.stderr)
    
output_dir = Path(inarg.output_dir).resolve()

if output_dir.is_dir():
    logging.info(f"output directory exists at: {output_dir}")
    print(f"output directory exists at: {output_dir}")
else:
    logging.error(f"error: ouput directory not found at {output_dir}")
    print(f"error: ouput directory not found at {output_dir}", file=sys.stderr)


if not(input_file.is_file()) or not(output_dir.is_dir()):
    logging.error("input file or output directory don't exist")
    sys.exit(1)

## create initial folders in output dir with tmp directory
output_path = output_dir / inarg.result_name
output_path.mkdir(parents=True, exist_ok=True)

# Paths
# inputDbDir = output_path / "inputDbDir"
refDbDir = output_path / "refDbDir"
refDbDir.mkdir(parents=True, exist_ok=True)

######################################
user = os.environ.get("USER", "nouser")
job_id = os.environ.get("SLURM_JOB_ID", "nojobid")

# Safely handle unset SCRATCH environment variable
scratch_env = os.environ.get("SCRATCH")
if scratch_env:
    scratch_base = Path(scratch_env)
else:
    scratch_base = Path(f"/scratch/general/vast/{user}/{job_id}")

# Use scratch_base, fallback to /tmp if needed
if scratch_base:
    tmpDir = scratch_base / "tmpDir"
else:
    tmpDir = Path(f"/tmp/{user}/test_scratch/tmpDir")

tmpDir.mkdir(parents=True, exist_ok=True)
os.environ["TMPDIR"] = str(tmpDir)

print(f"Scratch directory set up at: {tmpDir}")
######################################


#inputDbDir.mkdir(parents=True, exist_ok=True)




finalPath = output_path / "results"
finalPath.mkdir(parents=True, exist_ok=True)


mmseqs = "mmseqs"

def linear_steps(start, end, steps):
    return [round(start + i * (end - start) / (steps - 1), 4) for i in range(steps)]

id_steps = linear_steps(inarg.max_id, inarg.min_id, inarg.cascade_iter)
cov_steps = linear_steps(inarg.max_cov, inarg.min_cov, inarg.cascade_iter)

cascade_dir = output_path / "cascade_dir"
cascade_dir.mkdir(parents=True, exist_ok=True)

## create initial db
input_db = cascade_dir / "cascade_step0"
subp.run([mmseqs, "createdb", input_file, str(input_db)], check=True)

## cascaded clustering

cluster_dbs = []

prev_db = input_db
for i, (identity, coverage) in enumerate(zip(id_steps, cov_steps)):
    step = i + 1
    cluster_db = cascade_dir / f"cascade_clust{step}"
    rep_db = cascade_dir / f"cascade_rep{step}"
    
    print(f"Step {step}: Clustering at identity ≥ {identity:.2f}")
    
    # cluster
    subp.run([
        mmseqs, "cluster",
        str(prev_db), str(cluster_db), str(tmpDir),
        "--min-seq-id", str(identity),
        "-c", str(coverage),
        "--cov-mode", "0",
        "--cluster-mode", "0",
        "-s", str(inarg.cascade_sens),
        "-e", "1e-4"
    ], check=True)
    
    # save list of cluster_dbs
    cluster_dbs.append(cluster_db)
    
    # create subdb of representatives for next round
    subp.run([
        mmseqs, "createsubdb",
        str(cluster_db), str(cluster_db), str(prev_db), str(rep_db)
    ], check=True)
    
    prev_db = rep_db


print("\nMerging clusters from all cascade steps...")

merged_db = cluster_dbs[0]
for i, curr_db in enumerate(cluster_dbs[1:], start=1):
    new_merged_db = cascade_dir / f"merged_clust_step{i}"
    print(f"Merging step {i}: {merged_db} + {curr_db} → {new_merged_db}")

    subp.run([
        mmseqs, "mergeclusters",
        str(input_db),
        str(new_merged_db),
        str(merged_db),
        str(curr_db)
    ], check=True)

    merged_db = new_merged_db

print(f"\nFinal merged cluster DB: {merged_db}")


## convert cluster result to profiles

merged_rep_db = cascade_dir / "merged_rep_db"

subp.run([
    mmseqs, "createsubdb",
    str(merged_db),
    str(input_db),
    str(merged_rep_db)
], check=True)

input_db_h = input_db.with_name(input_db.name + "_h")
merged_rep_db_h = merged_rep_db.with_name(merged_rep_db.name + "_h")

subp.run([
    mmseqs, "createsubdb",
    str(merged_db),
    str(input_db_h),
    str(merged_rep_db_h)
], check=True)

merged_profile_db = cascade_dir / "merged_profile_db"

subp.run([
    mmseqs, "result2profile",
    str(merged_rep_db),
    str(input_db),
    str(merged_db),
    str(merged_profile_db),
    "--match-mode", "1",
    "--match-ratio", "1"
], check=True)


## profile-profile search
merged_profile_db_cons = cascade_dir / "merged_profile_db_cons"

subp.run([
    mmseqs, "profile2consensus",
    str(merged_profile_db),
    str(merged_profile_db_cons)
], check=True)

merged_profile_result = cascade_dir / "merged_profile_result"

subp.run([
    mmseqs, "search",
    str(merged_profile_db),
    str(merged_profile_db_cons),
    str(merged_profile_result),
    str(tmpDir),
    "--add-self-matches", "-a",
    "-s", str(inarg.search_sens),
    "-e", "1e-4",
    "-c", str(inarg.search_cov),
    "--cov-mode", "0",
    "--min-seq-id", "0.15"
], check=True)

merged_profile_clust = cascade_dir / "merged_profile_clust"

subp.run([
    mmseqs, "clust",
    str(merged_profile_db),
    str(merged_profile_result),
    str(merged_profile_clust),
    "--cluster-mode", "0"
], check=True)

subp.run([
    mmseqs, "createtsv",
    str(merged_profile_clust),
    str(merged_profile_clust),
    str(merged_profile_clust),
    "--cluster-mode", "0"
], check=True)





























