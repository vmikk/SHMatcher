#!/usr/bin/env python

# Script to run usearch single-linkage clustering for 90 percent clusters

import argparse
import csv
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description="Script to run usearch single-linkage clustering for 90 percent clusters")

parser.add_argument("--cluster", help="File for clustering")
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

# get cluster codes
with open(cl_tmp_file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        code = row[0]
        code_url = uclust_dir / name_folder / "clusters" / code
        mx_code = code + "_mx_03"
        mx_code_url = out_dir / mx_code
        out_code_005 = code + "_out_005"
        out_code_url_005 = out_dir / out_code_005

        # usearch -calc_distmx ClusterX -tabbedout mx_005.txt -maxdist 0.005 -threads 8
        usearch_cmd_1 = subprocess.run(
            [
                usearch_program,
                "-calc_distmx",
                code_url,
                "-tabbedout",
                mx_code_url,
                "-maxdist",
                "0.005",
                "-threads",
                threads,
            ],
            stdout=subprocess.DEVNULL,
        )

        # usearch -cluster_aggd mx_005.txt -clusterout clusters.txt -id 0.995 -linkage min
        usearch_cmd_2 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_005,
                "-id",
                "0.995",
                "-linkage",
                "max",
            ],
            stdout=subprocess.DEVNULL,
        )

        rm_cmd_1 = subprocess.run(["rm", str(mx_code_url)], stdout=subprocess.DEVNULL)
