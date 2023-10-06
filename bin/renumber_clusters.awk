
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

