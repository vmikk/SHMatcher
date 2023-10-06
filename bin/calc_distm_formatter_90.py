#!/usr/bin/env python

import argparse
import csv
# import logging
import os
from pathlib import Path
import subprocess

## Script to run usearch single-linkage clustering for 90 percent clusters

parser = argparse.ArgumentParser(description="Script to run usearch single-linkage clustering for 90 percent clusters")
parser.add_argument("--cluster",    help="File for clustering")
parser.add_argument("--uclust_dir", default="compounds",       help="Directory with clusters")
parser.add_argument("--out_dir",    default="calc_distm_out", help="Output directory")
parser.add_argument("--cl_tmp",     default="tmp.txt",        help="Temp file with clustering results")
parser.add_argument("--threads",    default=8,                help="CPU threads")

args = parser.parse_args()

# read in args
name = args.cluster
name_folder = name + "_folder"

# Infiles
uclust_dir  = args.uclust_dir # user_dir / "compounds"
out_dir     = args.out_dir    # uclust_dir / name_folder / "calc_distm_out"
cl_tmp_file = args.cl_tmp     # uclust_dir / name_folder / "tmp.txt"
threads     = args.threads    # 8

usearch_program = "usearch"

# get cluster codes
with open(cl_tmp_file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        code = row[0]
        # code_url = uclust_dir / name_folder / "clusters" / code
        code_url = os.path.join("clusters", code)
        mx_code = code + "_mx_03"
        mx_code_url = os.path.join(out_dir, mx_code)
        out_code_03 = code + "_out_03"
        out_code_url_03 = os.path.join(out_dir, out_code_03)
        out_code_025 = code + "_out_025"
        out_code_url_025 = os.path.join(out_dir, out_code_025)
        out_code_02 = code + "_out_02"
        out_code_url_02 = os.path.join(out_dir, out_code_02)
        out_code_015 = code + "_out_015"
        out_code_url_015 = os.path.join(out_dir, out_code_015)
        out_code_01 = code + "_out_01"
        out_code_url_01 = os.path.join(out_dir, out_code_01)
        out_code_005 = code + "_out_005"
        out_code_url_005 = os.path.join(out_dir, out_code_005)

        # usearch -calc_distmx ClusterX -tabbedout mx_03.txt -maxdist 0.03 -threads 8
        usearch_cmd_1 = subprocess.run(
            [usearch_program, "-calc_distmx", code_url, "-tabbedout", mx_code_url, "-maxdist", "0.03", "-threads", "8"],
            stdout=subprocess.DEVNULL,
        )

        # usearch -cluster_aggd mx_03.txt -clusterout clusters.txt -id 0.97 -linkage min
        usearch_cmd_2 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_03,
                "-id",
                "0.97",
                "-linkage",
                "min",
            ],
            stdout=subprocess.DEVNULL,
        )
        usearch_cmd_3 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_025,
                "-id",
                "0.975",
                "-linkage",
                "min",
            ],
            stdout=subprocess.DEVNULL,
        )
        usearch_cmd_4 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_02,
                "-id",
                "0.98",
                "-linkage",
                "min",
            ],
            stdout=subprocess.DEVNULL,
        )
        usearch_cmd_5 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_015,
                "-id",
                "0.985",
                "-linkage",
                "min",
            ],
            stdout=subprocess.DEVNULL,
        )
        usearch_cmd_6 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_01,
                "-id",
                "0.99",
                "-linkage",
                "min",
            ],
            stdout=subprocess.DEVNULL,
        )
        usearch_cmd_7 = subprocess.run(
            [
                usearch_program,
                "-cluster_aggd",
                mx_code_url,
                "-clusterout",
                out_code_url_005,
                "-id",
                "0.995",
                "-linkage",
                "min",
            ],
            stdout=subprocess.DEVNULL,
        )

        rm_cmd_1 = subprocess.run(["rm", str(mx_code_url)], stdout=subprocess.DEVNULL)
