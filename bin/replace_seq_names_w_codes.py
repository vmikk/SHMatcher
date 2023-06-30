#!/usr/bin/env python

## Script to replace original seq. ids with unique internal ids.

import argparse
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to replace original seq. ids with unique internal ids.")

# Read in args
parser.add_argument("--infile",     help="Input file")
parser.add_argument("--fasta_file", default="source_fasta", help="Output sequences")
parser.add_argument("--names_file", default="source_names", help="Output names")

args = parser.parse_args()

# user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile     = args.infile       # user_dir / f"source_{run_id}"
fasta_file = args.fasta_file   # user_dir / f"source_{run_id}_fasta"
names_file = args.names_file   # user_dir / f"source_{run_id}_names"

with open(infile, "r") as handle, open(fasta_file, "w") as f, open(names_file, "w") as n:
    for counter, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
        new_id = f"id_{counter}"

        f.write(f">i{new_id}i\n")
        f.write(f"{record.seq}\n")

        n.write(f"{record.id}\ti{new_id}i\t1\n")
