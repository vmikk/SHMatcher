#!/usr/bin/python

## Script to format vsearch --fastx_uniques output to restore original sequence.

import argparse
import csv
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Script to format vsearch --fastx_uniques output to restore original sequence."
)

parser.add_argument("--infile",  help="Input file")
parser.add_argument("--outfile", help="Output file")

args = parser.parse_args()

## Read in args
# user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile  = args.infile   # user_dir / f"source_{run_id}_fastauc"
outfile = args.outfile  # user_dir / f"source_{run_id}_fastanames"

matches_dict = {}
with open(infile, "r") as i:
    dataReader = csv.reader(i, delimiter="\t")
    for row in dataReader:
        if row[0] == "S":
            matches_dict[row[8]] = row[8]
        elif row[0] == "H":
            matches_dict[row[9]] = matches_dict[row[9]] + "," + row[8]

with open(outfile, "w") as o:
    for key in matches_dict:
        o.write(f"{key}\t{matches_dict[key]}\n")
