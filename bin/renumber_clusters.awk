
## Cluster ID Renumbering Script

### Description:
# This script is designed to concatenate multiple tab-delimited text files,
# where the first column of each file is a numeric cluster ID.
# The script ensures that cluster IDs remain unique across all concatenated files
# by renumbering them as needed.

### Usage:
# awk -f renumber_clusters.awk file1.txt file2.txt file3.txt
## or
# find . -name "*.txt" | sort --version-sort | xargs awk -f renumber_clusters.awk

### Input Format:
# Each input file should be a tab-delimited text file.
# The first column should contain numeric cluster IDs, starting from zero.
# The second column contains sequence ID (will not be modified).

# Example:
# ```
# 0	seq1
# 0	seq2
# 1	seq3
# 2	seq4
# ```

### Output:
# The script will print the concatenated files to the standard output
# with renumbered cluster IDs to ensure uniqueness.

### How it Works:
# 1. The script initializes an `offset` to 0.
# 2. For each line in the current file, the script checks if
#    the current cluster ID is greater than the previously observed maximum cluster ID (`max_id`).
#    If it is, `max_id` is updated.
# 3. The script prints the current line with the cluster ID offset by the value of `offset`.
# 4. At the end of each file, the `offset` is updated based on the maximum cluster ID observed in the file.
# 5. The process is repeated for each file.

### Notes:
# - The script assumes that cluster IDs in each file start from zero and are consecutive.
# - The script can handle any number of input files.


BEGIN {
	FS  = "\t";    # Set input delimiter to tab
    OFS = "\t";    # Set output delimiter to tab
    offset = 0
  }

# Process each line of the file
{
    # If the current cluster ID is greater than max_id, update max_id
    if ($1 > max_id) {
        max_id = $1
    }

    # Print the current line with the updated cluster ID
    print offset + $1, $2
}

# At the end of each file, update the offset for the next file
ENDFILE {
    offset += max_id + 1
    max_id = -1
}

