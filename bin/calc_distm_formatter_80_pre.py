#!/usr/bin/env python

# Script to run usearch complete-linkage clustering for 80 percent clusters

import argparse
import csv
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description="Script to run usearch complete-linkage clustering for 80 percent clusters")

parser.add_argument("--uclust", default="clusters",       help="Uclust direcotry")
parser.add_argument("--out",    default="calc_distm_out", help="Output directory")
parser.add_argument("--name",    help="Cluster name")
parser.add_argument("--threads", default=8, help="Number of CPUs")

args = parser.parse_args()

## Infiles
# user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
uclust_dir = args.uclust    # user_dir / "clusters_pre" / "clusters"
out_dir = args.out          # uclust_dir / "calc_distm_out"
threads = args.threads      # 8

usearch_program = "usearch"

code = args.name
code_url = uclust_dir / code
mx_code = code + "_mx_005"
mx_code_url = out_dir / mx_code
out_code_005 = code + "_out_005"
out_code_url_005 = out_dir / out_code_005

# Function to handle subprocesses
def run_subprocess(command):
    result = subprocess.run(command, stdout=subprocess.DEVNULL)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)

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
run_subprocess([
    usearch_program,
    "-cluster_aggd", mx_code_url,
    "-clusterout",   out_code_url_005,
    "-id",           "0.995",
    "-linkage",      "max"
])

## Clean up intermediate files
run_subprocess(["rm", str(mx_code_url)])
