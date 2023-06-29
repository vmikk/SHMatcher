// Species hypothesis (SH) matching analysis tool
//  Based on the SH MATCHING analysis tool: https://github.com/TU-NHM/sh_matching_pub


params.input = false

params.its_region = "itsfull"  // alternatively, "its2"

// Input sequence preparion:
// Replace sequence identifiers with unique codes,
// Remove duplicate sequences
process seq_prep {

    label "main_container"

    // cpus 1

    input:
      path input

    output:
      path "SeqQualities.txt.gz", emit: quals

    script:
    """
    echo -e "Aggregating sequence qualities"

    ## replace sequence identifiers with unique codes for the analysis
    python3 "$script_dir"/replace_seq_names_w_codes.py "$run_id"


    pushd "$user_dir"
    "$program_dir/vsearch/bin/vsearch" --fastx_uniques $infile_new_w_dir --fastaout "$infile_new""unique" --uc "$infile_new""uc"
    popd
    
    python3 "$script_dir"/reformat_fastx_uniques_output.py "$run_id"


    ## Additional quality controls - Remove low quality sequences (too short or with too many non-IUPAC symbols)
    python3 "$script_dir/exclude_non_iupac.py" "$run_id" 6



    echo -e "..Done"

    """
}

}


workflow {
}

