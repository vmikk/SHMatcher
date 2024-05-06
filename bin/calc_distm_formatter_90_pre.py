#!/usr/bin/env python

# Script to run usearch complete-linkage clustering for 90 percent clusters

# Main output is from `usearch -cluster_aggd -clusterout`,
# It has an `cluster output` file format = tabbed text with two fields: 1. cluster number and 2. label

import argparse
import csv
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description="Script to run usearch complete-linkage clustering for 90 percent clusters")

parser.add_argument("--cluster",    help="File for clustering")
parser.add_argument("--uclust_dir", default="clusters",       help="Directory with clusters")
parser.add_argument("--out_dir",    default="calc_distm_out", help="Output directory")
parser.add_argument("--cl_tmp",     default="tmp.txt",        help="Temp file with clustering results")
parser.add_argument("--threads",    default=8,                help="CPU threads")

args = parser.parse_args()

# read in args
name = args.cluster
name_folder = name + "_folder"

# Infiles
# user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
uclust_dir  = args.uclust_dir # user_dir / "clusters_pre" / "clusters"
out_dir     = args.out_dir    # uclust_dir / name_folder / "calc_distm_out"
cl_tmp_file = args.cl_tmp     # uclust_dir / name_folder / "tmp.txt"
threads     = args.threads    # 8

usearch_program = "usearch"

# Ensure the output directory exists
os.makedirs(out_dir, exist_ok=True)

# Function to handle subprocesses
def run_subprocess(command):
    result = subprocess.run(command, stdout=subprocess.DEVNULL)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)

# Get cluster codes
with open(cl_tmp_file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        code = row[0]
        # code_url = uclust_dir / name_folder / "clusters" / code
        code_url = os.path.join("clusters", code)
        mx_code = code + "_mx_03"
        mx_code_url = os.path.join(out_dir, mx_code)
        out_code_005 = code + "_out_005"
        out_code_url_005 = os.path.join(out_dir, out_code_005)

        ## Calculate distance matrix
        # usearch -calc_distmx ClusterX -tabbedout mx_005.txt -maxdist 0.005 -threads 8
        run_subprocess([
            usearch_program,
            "-calc_distmx", code_url,
            "-tabbedout",   mx_code_url,
            "-maxdist",     "0.005",
            "-threads",     threads,
        ])

        ## Complete linkage clustering
        # usearch -cluster_aggd mx_005.txt -clusterout clusters.txt -id 0.995 -linkage max
        run_subprocess([
            usearch_program,
            "-cluster_aggd", mx_code_url,
            "-clusterout",   out_code_url_005,
            "-id",           "0.995",
            "-linkage",      "max"
        ])

        ## Clean up intermediate files
        run_subprocess(["rm", str(mx_code_url)])
