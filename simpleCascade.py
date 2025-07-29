# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 16:45:48 2025

@author: imagi
"""

import argparse as ap
# import numpy as np
import subprocess as subp
from pathlib import Path
import sys
import networkx as nx
from collections import defaultdict
#from typing import List, Set, Tuple
import logging
import os
#import multiprocessing

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

parser.add_argument("--search-cov",
                    action="store",
                    type=float,
                    dest="search_cov",
                    help="decimal coverage threshold for strict clustering", 
                    default=float(0.7))

parser.add_argument("--search-id",
                    action="store",
                    type=float,
                    dest="search_id",
                    help="decimal coverage threshold for strict clustering", 
                    default=float(0.1))

parser.add_argument("--cascade-eval",
                    action="store",
                    type=str,
                    dest="cascade_eval",
                    help="decimal coverage threshold for strict clustering", 
                    default="1e-5")

parser.add_argument("--search-eval",
                    action="store",
                    type=str,
                    dest="search_eval",
                    help="decimal coverage threshold for strict clustering", 
                    default="1e-5")

# optional optional arguments


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

# If available, use CHPC's SCRATCH (which is usually /scratch/general/vast/$USER)
scratch_base = Path(os.environ.get("SCRATCH", f"/scratch/general/vast/{user}"))

# Use job-specific tmp directory
if job_id != "nojobid":
    tmpDir = scratch_base / job_id / "tmpDir"
else:
    tmpDir = scratch_base / "manual_job" / "tmpDir"

tmpDir.mkdir(parents=True, exist_ok=True)
os.environ["TMPDIR"] = str(tmpDir)

print(f"Temporary directory created at: {tmpDir}")

num_threads = os.environ.get("SLURM_CPUS_PER_TASK")
# num_threads = str(multiprocessing.cpu_count())
######################################


#inputDbDir.mkdir(parents=True, exist_ok=True)


finalPath = output_path / "results"
finalPath.mkdir(parents=True, exist_ok=True)

mmseqs = "mmseqs"
muscle = "muscle"

cascade_dir = output_path / "cascade_dir"
cascade_dir.mkdir(parents=True, exist_ok=True)

## create initial db
input_db = cascade_dir / "cascade_step0"
subp.run([mmseqs, "createdb", input_file, str(input_db)], check=True)


## cascaded clustering
cluster_db = cascade_dir / "cascade_cluster"
subp.run([
    mmseqs, "cluster",
    str(input_db), str(cluster_db), str(tmpDir),
    "--min-seq-id", str(inarg.min_id),
    "-c", str(inarg.min_cov),
    "--cov-mode", "0",
    "--cluster-mode", "0",
    "-s", str(inarg.cascade_sens),
    "-e", str(inarg.cascade_eval),
    "--cluster-steps", "3",
    "--cluster-reassign"
], check=True)

## print tsv

cascade_tsv = finalPath / "cascade-sequence_tsv.tsv"

subp.run([
    mmseqs, "createtsv",
    str(input_db),
    str(input_db),
    str(cluster_db),
    str(cascade_tsv)
], check=True)



## convert cluster result to profiles

rep_db = cascade_dir / "rep_db"

subp.run([
    mmseqs, "createsubdb",
    str(cluster_db),
    str(input_db),
    str(rep_db)
], check=True)

input_db_h = input_db.with_name(input_db.name + "_h")
rep_db_h = rep_db.with_name(rep_db.name + "_h")

subp.run([
    mmseqs, "createsubdb",
    str(cluster_db),
    str(input_db_h),
    str(rep_db_h)
], check=True)

profile_dir = output_path / "profile_dir"
profile_dir.mkdir(parents=True, exist_ok=True)

profile_db = profile_dir / "profile_db"

subp.run([
    mmseqs, "result2profile",
    str(rep_db),
    str(input_db),
    str(cluster_db),
    str(profile_db),
], check=True)


## profile-profile search
profile_db_consensus = profile_dir / "profile_db_consensus"

subp.run([
    mmseqs, "profile2consensus",
    str(profile_db),
    str(profile_db_consensus)
], check=True)

profile_result = profile_dir / "profile_result"

subp.run([
    mmseqs, "search",
    str(profile_db),
    str(profile_db_consensus),
    str(profile_result),
    str(tmpDir),
    "--add-self-matches", "-a",
    "-s", str(inarg.search_sens),
    "-e", str(inarg.search_eval),
    "-c", str(inarg.search_cov),
    "--cov-mode", "0",
    "--min-seq-id", str(inarg.search_id)
], check=True)

#################################
        # write m8 file #
#################################

profile_mate = profile_dir / "profile_mate_result"

subp.run([
    mmseqs, "convertalis",
    str(profile_db),
    str(profile_db_consensus),
    str(profile_result),
    str(profile_mate)
], check=True)

hit_pairs = set()
all_seqs = set()
searchGraph = nx.Graph()

with open(profile_mate, "r") as f:
    for line in f:
        query, target = line.strip().split("\t")[0:2]
        all_seqs.update([query, target])
        if query != target:
            hit_pairs.add((query, target))

for query, target in hit_pairs:
    if (target, query) in hit_pairs:
        searchGraph.add_edge(query, target)

searchGraph.add_nodes_from(all_seqs)

con_comp = sorted(nx.connected_components(searchGraph), key=len, reverse=True)

# write a tsv

profile_tsv = finalPath / "profileclust-cascade_tsv.tsv"

with open(profile_tsv, "w") as out:
    for i, components in enumerate(con_comp):
        proclust_id = f"proclust_{i:06d}"
        for component in components:
            out.write(f"{proclust_id}\t{component}\n")


# double tsv index 


cascade_dict = defaultdict(list)

with open(cascade_tsv, "r") as f:
    for line in f:
        rep, member = line.strip().split("\t")
        cascade_dict[rep].append(member)

profile_dict = defaultdict(list)

with open(profile_tsv, "r") as f:
    for line in f:
        rep, member = line.strip().split("\t")
        profile_dict[rep].append(member)

proseq_tsv = finalPath / "connectedprofile-sequence_tsv.tsv"

with open(proseq_tsv, "w") as out:
    for rep, intereps in profile_dict.items():
        for interep in intereps:
            final_members = cascade_dict.get(interep, [])
            for member in final_members:
                out.write(f"{rep}\t{member}\n")

input_db_lookup = cascade_dir / "cascade_step0.lookup"



# loop through tsv, to write in idx using already loaded lookup

lookup_dict = {}
with open(input_db_lookup, "r") as lookup_file:
    for line in lookup_file:
        idx, match_name, col3 = line.strip().split("\t")
        lookup_dict[match_name] = idx

rep_to_indices = defaultdict(list)

with open(proseq_tsv, "r") as tsv:
    for line in tsv:
        rep, seq = line.strip().split("\t")
        idx = lookup_dict.get(seq)
        if idx:
            rep_to_indices[rep].append(idx)

sorted_reps = sorted(rep_to_indices.items(), key=lambda item: len(item[1]), reverse=True)

index_path = output_path / "profile_indices"
index_path.mkdir(parents=True, exist_ok=True)


for i, (rep, indices) in enumerate(sorted_reps):
    cluster_id = f"cluster_{i:06d}.txt"  # e.g. cluster_000000.txt
    filepath = index_path / cluster_id
    with open(filepath, "w") as out:
        out.write("\n".join(indices) + "\n")

# now loop through reps and subdb and print fastas
index_file_list = [f for f in index_path.glob("*.txt")]

cluster_subdb = output_path / "cluster_subdb"
cluster_subdb.mkdir(parents=True, exist_ok=True)

cluster_fastas = finalPath / "cluster_fastas"
cluster_fastas.mkdir(parents=True, exist_ok=True)

for file_name in index_file_list:
    subdb_path = cluster_subdb / file_name.stem
    subp.run([
        mmseqs, "createsubdb",
        str(file_name),
        str(input_db),
        str(subdb_path)
    ])
    clust_fast_path = cluster_fastas / f"{file_name.stem}.fasta"
    subp.run([
        mmseqs, "convert2fasta",
        str(subdb_path),
        str(clust_fast_path)
    ])

fasta_file_list = [f for f in cluster_fastas.glob("*.fasta")]

cluster_MSAs = finalPath / "final_cluster_MSAs"
cluster_MSAs.mkdir(parents=True, exist_ok=True)

# now align all fastas with muscle

for i, file_name in enumerate(fasta_file_list):
    cluster_ids = f"{file_name.stem}.afa"
    outfilepath =  cluster_MSAs / cluster_ids
    subp.run([
        muscle, "-super5", str(file_name),
        "-output", str(outfilepath),
        "-threads", str(num_threads)
    ])
    
cmd = f"""
find {cluster_MSAs} -type f -name "*.afa" -exec bash -c 'echo -e "$(basename "{{}}")\t$(grep -c "^>" "{{}}")"' \; | sort -k2,2nr > {finalPath}/cluster_member_counts.tsv
"""
subp.run(cmd, shell=True, executable="/bin/bash")
