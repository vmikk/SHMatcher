#!/usr/bin/env python

import argparse
import csv
import os
from pathlib import Path

## Script to merge parsed matches output into one file

parser = argparse.ArgumentParser(description="Script to merge parsed matches output into one file")

parser.add_argument("--matches", default="matches", help="Matches dir")
parser.add_argument("--outfile", default="matches_out_all.csv", help="Output file")
parser.add_argument("--besthits", default="closedref.80-best-hits.map.uc", help="UC file with best hits")

## Matches dir must contain  `matches_out_{thld}.csv`

args = parser.parse_args()

# user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
matches_dir    = args.matches    # user_dir / "matches"
outfile        = args.outfile    # matches_dir / f"matches_out_all.csv"
best_hits_file = args.besthits   # user_dir / "closedref.80-best-hits.map.uc"

threshold_list = ["03", "025", "02", "015", "01", "005"]
one_line_str_dict = {}

print("Starting script execution...")

# get best hits into mapping and similarity tables
print(f"Reading best hits from {best_hits_file}...")
bh_mapping_dict = {}
bh_sim_dict = {}
with open(best_hits_file) as bh:
    dataReader = csv.reader(bh, delimiter="\t")
    for row in dataReader:
        if row[0] == "H":
            query_id = row[8].replace("i", "")
            subject_fields = row[9].split("_")
            plutof_seq_url = "https://app.plutof.ut.ee/sequence/view/" + str(subject_fields[0].replace("i", ""))
            bh_mapping_dict[query_id] = plutof_seq_url
            bh_sim_dict[query_id] = str(row[3])
print("Best hits loaded successfully.")

# merge threshold-based files
print("Merging threshold-based files...")
with open(outfile, "w") as o:
    for thld in threshold_list:
        matches_file = Path(matches_dir) / f"matches_out_{thld}.csv"
        if matches_file.is_file():
            print(f"Processing file: {matches_file}")
            with open(matches_file) as m:
                dataReader_m = csv.reader(m, delimiter="\t")
                row_ct = 0
                for row in dataReader_m:
                    row_ct += 1
                    if thld == "03":
                        one_line_str_dict[row[0]] = ""
                        ## seq_id_tmp   seq_accno   status (3.0)    SH code (3.0)   SH/compound taxonomy (3.0)
                        for i in range(5):
                            one_line_str_dict[row[0]] += row[i] + "\t"
                    elif thld == "005":
                        ## status (1.0) SH code (1.0)   SH/compound taxonomy (1.0)    compound_cl_code (1.0)    Compound taxonomy (1.0)
                        if row_ct == 1:
                            # header line, add columns for best hits
                            one_line_str_dict[row[0]] += (
                                row[2]
                                + "\t"
                                + row[3]
                                + "\t"
                                + row[4]
                                + "\t"
                                + row[5]
                                + "\t"
                                + row[6]
                                + "\tMatched sequence\tSimilarity percentage"
                                + "\n"
                            )
                        else:
                            one_line_str_dict[row[0]] += (
                                row[2] + "\t" + row[3] + "\t" + row[4] + "\t" + row[5] + "\t" + row[6]
                            )
                            if row[0] in bh_mapping_dict:
                                one_line_str_dict[row[0]] += "\t" + bh_mapping_dict[row[0]] + "\t" + bh_sim_dict[row[0]]
                            else:
                                if row[7] and row[7] in bh_mapping_dict:
                                    one_line_str_dict[row[0]] += (
                                        "\t" + bh_mapping_dict[row[7]] + "\t" + bh_sim_dict[row[7]]
                                    )
                                else:
                                    one_line_str_dict[row[0]] += "\t" + "\t"
                            one_line_str_dict[row[0]] += "\n"
                    else:
                        ## status (2.0) SH code (2.0)   SH/compound taxonomy (2.0)
                        one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + row[4] + "\t"
    print("Finished merging data from threshold-based files.")

    print("Processing additional threshold-based files...")
    for thld in threshold_list:
        matches_1_file = Path(matches_dir) / f"matches_1_out_{thld}.csv"
        if matches_1_file.is_file():
            with open(matches_1_file) as m1:
                dataReader_m1 = csv.reader(m1, delimiter="\t")
                row_ct = 0
                for row in dataReader_m1:
                    row_ct += 1
                    if row_ct > 1:
                        if thld == "03":
                            one_line_str_dict[row[0]] = ""
                            ## seq_id_tmp  seq_accno   status (3.0)    SH code (3.0)
                            for i in range(4):
                                one_line_str_dict[row[0]] += row[i] + "\t"
                            one_line_str_dict[row[0]] += "\t"
                        elif thld == "005":
                            ## status (0.05)    SH code (0.05)   compound_cl_code (0.05)
                            one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + "\t" + row[4]
                            if row[0] in bh_mapping_dict:
                                one_line_str_dict[row[0]] += "\t" + bh_mapping_dict[row[0]] + "\t" + bh_sim_dict[row[0]]
                            else:
                                if row[5] and row[5] in bh_mapping_dict:
                                    one_line_str_dict[row[0]] += (
                                        "\t" + bh_mapping_dict[row[5]] + "\t" + bh_sim_dict[row[5]]
                                    )
                                else:
                                    one_line_str_dict[row[0]] += "\t" + "\t"
                            one_line_str_dict[row[0]] += "\n"
                        else:
                            ## status (2.0)    SH code (2.0)
                            one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + "\t"

    print("All data merged successfully. Writing to output file...")
    for line in one_line_str_dict:
        o.write(one_line_str_dict[line])

    print(f"Script execution completed. Data written to {outfile}.")
