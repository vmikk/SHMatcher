#!/usr/bin/env python

## Script to select from max 5 best usearch hits the one with the highest identity)

import argparse
import csv
import logging
import os
from pathlib import Path

from Bio import SeqIO
from decimal import Decimal

parser = argparse.ArgumentParser(
    description="Script to select from max 5 best usearch hits the one with the highest identity)"
)

parser.add_argument("--infile",       default="iupac_out_full.fasta", help="Input fasta")
parser.add_argument("--hits_fasta",   default="hits.fasta",   help="Output - hits, sequences")
parser.add_argument("--nohits_fasta", default="nohits.fasta", help="Output - no hits, sequences")
parser.add_argument("--hits",      default="hits.txt",        help="Output - hits")
parser.add_argument("--map_file",  default="closedref.80.map.uc",           help="Mapping file")
parser.add_argument("--best_hits", default="closedref.80-best-hits.map.uc", help="Best hits")
parser.add_argument("--log",       default="err.log", help="Log file")


args = parser.parse_args()

# Read in args
## user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile   = args.infile       # user_dir / "iupac_out_full.fasta"
outfile1 = args.hits_fasta   #  user_dir / "hits.fasta"
outfile2 = args.nohits_fasta #  user_dir / "nohits.fasta"
outfile3 = args.hits         #  user_dir / "hits.txt"
map_file = args.map_file     #  user_dir / "closedref.80.map.uc"
## print out one best hit for each sequence to create compound clusters in later stage
best_hits_file = args.best_hits #  user_dir / "closedref.80-best-hits.map.uc"
log_file = args.log_file #   user_dir / f"err_{run_id}.log"

# Logging conf

logging.basicConfig(
    filename=log_file,
    filemode="a",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    level="INFO",
)

full_dict = {}
with open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        full_dict[name] = str(record.seq)

ident_dict = {}
best_match_dict = {}
best_match_dict_full = {}

with open(outfile2, "w") as o2, open(map_file) as map_f:
    # TODO - possibility for csv.DictReader
    dataReader_map = csv.reader(map_f, delimiter="\t")
    for row in dataReader_map:
        if row[0] == "N":
            o2.write(f">{row[8]}\n")
            o2.write(f"{full_dict[row[8]]}\n")
        elif row[0] == "H":
            if row[8] in best_match_dict:
                # at least one of sequences best matches already processed
                if Decimal(row[3]) > Decimal(ident_dict[row[8]]):
                    best_match_dict[row[8]] = f"{row[8]}\t{row[2]}\t{row[3]}\t{row[9]}\n"
                    best_match_dict_full[row[8]] = row
                    ident_dict[row[8]] = row[3]
            else:
                best_match_dict[row[8]] = f"{row[8]}\t{row[2]}\t{row[3]}\t{row[9]}\n"
                best_match_dict_full[row[8]] = row
                ident_dict[row[8]] = row[3]
        else:
            logging.info(f"USEARCH\tERR: {row[0]}")

with open(outfile1, "w") as o1, open(outfile3, "w") as o3, open(best_hits_file, "w") as bh:
    for best_match in best_match_dict:
        o3.write(best_match_dict[best_match])
        o1.write(f">{best_match}\n")
        o1.write(f"{full_dict[best_match]}\n")
        bh_string = "\t".join(best_match_dict_full[best_match])
        bh.write(bh_string + "\n")
