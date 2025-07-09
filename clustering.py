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

parser.add_argument("--loose-clust-id",
                    action="store",
                    type=float,
                    dest="loose_id",
                    help="decimal identity threshold for permissive clustering", 
                    default=float(0.35))

parser.add_argument("--loose-clust-cov",
                    action="store",
                    type=float,
                    dest="loose_cov",
                    help="decimal cov threshold for permissive clustering", 
                    default=float(0.7))

parser.add_argument("--strict-clust-id",
                    action="store",
                    type=float,
                    dest="strict_id",
                    help="decimal identity threshold for strict clustering", 
                    default=float(0.9))

parser.add_argument("--strict-clust-cov",
                    action="store",
                    type=float,
                    dest="strict_cov",
                    help="decimal coverage threshold for strict clustering", 
                    default=float(0.9))

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

parser.add_argument("--cluster-sensitivity",
                    action="store",
                    type=float,
                    dest="clust_sens",
                    help="desired mmseqs sequence clustering sensitivity (1.0 ≤ s ≤ 9.5)", 
                    default=float(7.5))

parser.add_argument("--search-sensitivity",
                    action="store",
                    type=float,
                    dest="search_sens",
                    help="desired mmseqs profile-sequence search sensitivity (1.0 ≤ s ≤ 9.5)", 
                    default=float(7.5))


parser.add_argument("--max-iter-strict",
                    action="store",
                    type=int,
                    dest="max_iter_strict",
                    help="max number of strict clustering/search iterations till convergence", 
                    default=int(10))

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
    
output_resolved = Path(inarg.output_dir).resolve()

if output_resolved.is_dir():
    logging.info(f"output directory exists at: {output_resolved}")
    print(f"output directory exists at: {output_resolved}")
else:
    logging.error(f"error: ouput directory not found at {output_resolved}")
    print(f"error: ouput directory not found at {output_resolved}", file=sys.stderr)


if not(input_file.is_file()) or not(output_resolved.is_dir()):
    logging.error("input file or output directory don't exist")
    sys.exit(1)

## create initial folders in output dir with tmp directory
output_path = output_resolved / inarg.result_name
output_path.mkdir(parents=True, exist_ok=True)

# Paths
inputDbDir = output_path / "inputDbDir"
tmpDir = output_path / "tmpDir"
refDbDir = output_path / "refDbDir"


#inputDbDir.mkdir(parents=True, exist_ok=True)
tmpDir.mkdir(parents=True, exist_ok=True)
refDbDir.mkdir(parents=True, exist_ok=True)



finalPath = output_path / "results"
finalPath.mkdir(parents=True, exist_ok=True)

log_path = finalPath / f"{inarg.result_name}.log"
logging.basicConfig(
    filename=log_path,
    level=logging.INFO,        # log messages INFO and above
    format="%(asctime)s - %(levelname)s - %(message)s"
)

###############################################################################

## define necessary functions


def mmCluster(input_name: Path, fasta_or_db, iter_name: str, out_direct: Path = output_path, inarg=inarg):
    stems = input_name.stem
    genPath = out_direct / iter_name / stems
    for fold in ["Db", "Clu", "SubDb", "Seq", "Out", "RepDb"]:
        dirPath = genPath / fold
        dirPath.mkdir(parents=True, exist_ok=True)
    
    if iter_name == "perm":
        ident = inarg.loose_id
        cov = inarg.loose_cov
    else:
        ident = inarg.strict_id
        cov = inarg.strict_cov
        
    if fasta_or_db == "fasta":
        inDataPath = genPath / "Db" / stems
        cmdOne = f"mmseqs createdb {input_name} {inDataPath}"
        subp.run(cmdOne, shell=True, check=True, executable="/bin/bash")
    else:
        inDataPath = input_name
    if iter_name == "perm":
        global fullInputDatabase
        fullInputDatabase = inDataPath
    
    cmdTwo = f"""
    mmseqs cluster {inDataPath} {genPath}/Clu/{stems} {tmpDir} --min-seq-id {ident} -c {cov} -s {inarg.clust_sens}
    mmseqs createtsv {inDataPath} {inDataPath} {genPath}/Clu/{stems} {genPath}/Out/{stems}.tsv
    mmseqs createsubdb {genPath}/Clu/{stems} {inDataPath} {genPath}/RepDb/{stems}
    
    cut -f1 {genPath}/Out/{stems}.tsv | sort | uniq > {genPath}/Out/{stems}_uniqueID.tsv
    
    i=1
    while read -r rep; do
        cluster_id=$(printf "Clust_%06d" "$i")

        awk -v r="$rep" '$1 == r {{ print $2 }}' {genPath}/Out/{stems}.tsv > {tmpDir}/${{cluster_id}}_ids.txt
        
        awk 'NR==FNR {{a[$1]; next}} $2 in a {{print $1}}' {tmpDir}/${{cluster_id}}_ids.txt {inDataPath}.lookup > {tmpDir}/${{cluster_id}}_indices.txt

        
        mmseqs createsubdb {tmpDir}/${{cluster_id}}_indices.txt {inDataPath} {genPath}/SubDb/${{cluster_id}}
        mmseqs convert2fasta {genPath}/SubDb/${{cluster_id}} {genPath}/Seq/${{cluster_id}}.fasta

        i=$((i + 1))
    done < {genPath}/Out/{stems}_uniqueID.tsv
    """
    subp.run(cmdTwo, shell=True, check=True, executable="/bin/bash")


def muscleAlign(input_name: Path, iter_name: str, out_direct: Path): #will need to align each 90percent SubDb
    stems = input_name.stem
    genPath = out_direct / iter_name / stems # / "Seq" / "90clusters"
    
    dirPath = genPath / "Align"
    dirPath.mkdir(parents=True, exist_ok=True)
    
    ## gather the filenames for aligning
    
    
    databaseDir = genPath / "Seq"  # define database
    clusterList = [f.resolve().stem for f in databaseDir.glob("*.fasta")]
    for clusters in clusterList:
        cmdOne = f"muscle -super5 {genPath}/Seq/{clusters}.fasta -output {genPath}/Align/{clusters}.afa"
        cmdTwo = f"perl {REFORMAT_PL} fas sto {genPath}/Align/{clusters}.afa {genPath}/Align/{clusters}.sto"
        subp.run(cmdOne, shell=True, check=True, executable="/bin/bash")
        subp.run(cmdTwo, shell=True, check=True, executable="/bin/bash")


def stoicProfile(input_name: Path, iter_name: str, out_direct: Path):
    stems = input_name.stem
    genPath = out_direct / iter_name / stems # / "Seq" / "90clusters"
    
    dirPath = genPath / "AlSubDb"
    dirTwoPath = genPath / "AlProfile"
    dirPath.mkdir(parents=True, exist_ok=True)
    dirTwoPath.mkdir(parents=True, exist_ok=True)
    # just calling 90percent clusters subdbs
    
    
    databaseDir = genPath / "Align"  # define database
    clusterList = [f.resolve().stem for f in databaseDir.glob("*.sto")]
    
    for clusters in clusterList:
        cmdOne = f""" 
        mmseqs convertmsa {genPath}/Align/{clusters}.sto {genPath}/AlSubDb/{clusters}
        mmseqs msa2profile {genPath}/AlSubDb/{clusters} {genPath}/AlProfile/{clusters} --match-mode 1
        """
        subp.run(cmdOne, shell=True, check=True, executable="/bin/bash")
        
# subject of search in {genPath}/RepDb/{stems}
# query of search in {genPath}/AlProfile/{clusters}
def mmSearch(input_name: Path, iter_name: str, out_direct: Path, inarg=inarg):
    stems = input_name.stem
    genPath = out_direct / iter_name / stems
    
    dirPath = genPath / "ProRepSeqSearchMate"
    dirPath.mkdir(parents=True, exist_ok=True)
    
    dirDbPath = genPath / "ProRepSeqSearchDb"
    dirDbPath.mkdir(parents=True, exist_ok=True)
    
    sens = inarg.search_sens
    
    # collect profile names to loop through -_h
    databaseDir = genPath / "AlProfile" # define database
    clusterList = [
        f.stem  # get filenames
        for f in databaseDir.glob("*.dbtype")
        if not f.name.endswith("_h.dbtype")
    ]
    
    
    # command to be looped
    
    for clusters in clusterList:
        cmdOne = f"""
        mmseqs search {genPath}/AlProfile/{clusters} {genPath}/RepDb/{stems} {genPath}/ProRepSeqSearchDb/{clusters} {tmpDir} -s {sens} --num-iterations 4 -e 1e-4 -c 0.7
        mmseqs convertalis {genPath}/AlProfile/{clusters} {genPath}/RepDb/{stems} {genPath}/ProRepSeqSearchDb/{clusters} {genPath}/ProRepSeqSearchMate/{clusters}.m8
        """
        subp.run(cmdOne, shell=True, check=True, executable="/bin/bash")
    
    # combines all .m8 files across iteration
    cmdTwo = f"cat {out_direct}/{iter_name}/*/ProRepSeqSearchMate/*.m8 > {out_direct}/{iter_name}/all_iter_hits.m8"
    subp.run(cmdTwo, shell=True, check=True, executable="/bin/bash")
    
    

def conComp(iter_name: str, out_direct: Path) -> Tuple[List[Set[str]], List[Path], Path]:
    
    dirPath = out_direct / iter_name / "conComp"
    dirPath.mkdir(parents=True, exist_ok=True)
    
    dirPathDb = out_direct / iter_name / "conCompDb"
    dirPathDb.mkdir(parents=True, exist_ok=True)
    searchGraph = nx.Graph()
    all_hits = out_direct / iter_name / "all_iter_hits.m8"
    with open(all_hits) as f:
        for line in f:
            parts = line.strip().split("\t")
            query = parts[0]
            target = parts[1]
            if query != target:
                searchGraph.add_edge(query,target)
    con_comp = sorted(nx.connected_components(searchGraph), key=len, reverse=True)
    
    
    # iterate over all sets in con_comp
    for i, component in enumerate(con_comp):
        cluster_id = f"cluster_{i:07d}"
        id_paths = dirPath / f"{cluster_id}_ids.txt"
        index_paths = dirPath / f"{cluster_id}_indices.txt"
        subdb_paths = dirPathDb / cluster_id
        
        # write ids to file
        with open(id_paths, "w") as f:
            for node in component:
                f.write(f"{node}\n")
        
        # bash awk command
        cmdOne = f"""
        awk 'NR==FNR {{a[$1]; next}} $2 in a {{print $1}}' {id_paths} {fullInputDatabase}.lookup > {index_paths}
        mmseqs createsubdb {index_paths} {fullInputDatabase} {subdb_paths}
        """
        subp.run(cmdOne, shell=True, check=True, executable="/bin/bash")
        
        # redefine datab_names
    datab_names = [
        f.with_suffix('')
        for f in dirPathDb.glob("*.dbtype")
        if not f.name.endswith("_h.dbtype")]
        
    return con_comp, datab_names, dirPathDb


def mmTsvToSets(tsv_path: Path):
    cluster_dict = defaultdict(set)
    with open(tsv_path) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            cluster_dict[rep].add(member)
            
    return list(cluster_dict.values())

def jaccard(set1, set2):
    return len(set1 & set2) / len(set1 | set2)


# search output concatenation commands will include input for previous and next iteration so that the next iteration is compatible with previous functions

###############################################################################


## permissive clustering
perm_iter = "perm"

logging.info("beginning permissive clustering")
mmCluster(input_file, "fasta", perm_iter, out_direct=output_path)
logging.info("permissive clustering finished")

initial_tsv_clusters = output_path / perm_iter / input_file.stem / "Out" / f"{input_file.stem}.tsv"

cluster_set1 = mmTsvToSets(initial_tsv_clusters)

## iterative strict clustering and search

#define paths
stemor = input_file.stem
databaseDir = output_path / "perm" / stemor / "SubDb"
#datab_names = [f.with_suffix('') for f in databaseDir.glob("*.dbtype")]
datab_names = [
    f.with_suffix('')  # strip suffix, keep full path
    for f in databaseDir.glob("*.dbtype")
    if not f.name.endswith("_h.dbtype")
]


iter_var = 0
weighted_jaccard = 0

while weighted_jaccard < 1 and iter_var < inarg.max_iter_strict:
    logging.info(f"strict clustering iteration {iter_var}")
    strict_iter = f"strict{iter_var}"
    for databayes in datab_names:
        mmCluster(databayes, "db", strict_iter, out_direct=output_path)
        muscleAlign(databayes, strict_iter, out_direct=output_path)
        stoicProfile(databayes, strict_iter, out_direct=output_path)
        mmSearch(databayes, strict_iter, out_direct=output_path)
        
        
    cluster_set2, datab_names, finalDbDir = conComp(strict_iter, out_direct=output_path)
    
    total_score = 0
    total_weights = 0
    
    logging.info(f"calculating jaccard values for iteration {iter_var}")
    for ca in cluster_set2:
        best_match = max(jaccard(ca, cb) for cb in cluster_set1)
        total_score += best_match * len(ca)
        total_weights += len(ca)
        
    weighted_jaccard = total_score / total_weights if total_weights else 0
    logging.info(f"weighted jaccard value: {weighted_jaccard}")
    cluster_set1 = [set(s) for s in cluster_set2]
    
    iter_var += 1


logging.info(f"total strict iterations: {iter_var + 1}")


## final output folders
    ## in finalPath folder
finalFasta = finalPath / "cluster_fasta"
finalFasta.mkdir(parents=True, exist_ok=True)

finalMSA = finalPath / "cluster_MSAs"
finalSto = finalMSA / "sto"
stoMsaDb = finalSto / "msaDb"
finalPro = finalPath / "cluster_mmProfiles"

finalMSA.mkdir(parents=True, exist_ok=True)
finalSto.mkdir(parents=True, exist_ok=True)
stoMsaDb.mkdir(parents=True, exist_ok=True)
finalPro.mkdir(parents=True, exist_ok=True)

# unaligned fasta files
    # use finalDbDir or datab_names?
for databayes in datab_names:
    cmdFasta = f"""
    mmseqs convert2fasta {databayes} {finalFasta}/{databayes.stem}.fasta
    """
    subp.run(cmdFasta, shell=True, check=True, executable="/bin/bash")

# muscle aligned clusters
    # align the previous

fasta_files = [
    f for f in finalFasta.glob("*.fasta")]

for file in fasta_files:
    cmdMuscle = f"""
    muscle -super5 {file} -output {finalMSA}/{file.stem}.afa
    perl reformat.pl fas sto {finalMSA}/{file.stem}.afa {finalSto}/{file.stem}.sto
    """
    subp.run(cmdMuscle, shell=True, check=True, executable="/bin/bash")


# mmseqs (profile?) database files
    # profile it

stoAwayList = [f.resolve().stem for f in finalSto.glob("*.sto")]


for stowaways in stoAwayList:
    cmdOne = f""" 
    mmseqs convertmsa {finalSto}/{stowaways}.sto {stoMsaDb}/{stowaways}
    mmseqs msa2profile {stoMsaDb}/{stowaways} {finalPro}/{stowaways} --match-mode 1
    """
    subp.run(cmdOne, shell=True, check=True, executable="/bin/bash")


# index files
    # cluster NAME and members .tsv file
    # make from cluster_set2

finaltsv = finalPath / "final_sorted_cluster_index.tsv"
with open(finaltsv, "w") as f:
    for i, cluster in enumerate(cluster_set2):
        label = f"cluster_{i:07d}"
        f.writelines(map(lambda member: f"{label}\t{member}\n", cluster))
    
logging.info("done")


###############################################################################




































