// Species hypothesis (SH) matching analysis tool
//  Based on the SH MATCHING analysis tool: https://github.com/TU-NHM/sh_matching_pub


params.input = false

params.its_region = "itsfull"  // alternatively, "its2"

params.db     = "sanger_refs_sh.udb"
params.dbfull = "sanger_refs_sh_full.udb"



// Input sequence preparion:
// Replace sequence identifiers with unique codes,
// Remove duplicate sequences
process seq_prep {

    label "main_container"

    // cpus 1

    input:
      path input

    output:
      path "source_names",               emit: namesorig
      path "source_fasta_unique",        emit: unique
      path "source_fasta_uc",            emit: uc
      path "source_fasta_names",         emit: namesuniq

    script:
    """
    echo -e "Input sequence preparion\n"

    ## Replace sequence identifiers with unique codes for the analysis
    echo -e "\n..Replacing sequence names with codes\n"
    replace_seq_names_w_codes.py \
      --infile ${input} \
      --fasta_file source_fasta \
      --names_file source_names

    ## Find the set of unique sequences in an input file
    echo -e "\n..Dereplicating sequences\n"
    vsearch \
      --fastx_uniques source_fasta \
      --fastaout      source_fasta_unique \
      --uc            source_fasta_uc

    ## Format output to restore original sequence
    echo -e "\n..Reformatting sequences\n"
    reformat_fastx_uniques_output.py \
      --infile  source_fasta_uc \
      --outfile source_fasta_names

    echo -e "..Done"
    """
}

// Chimera filtering
process chimera_filtering {

    label "main_container"

    cpus 8

    input:
      path input
      path db        // sanger_refs_sh.udb

    output:
      path "seqs_out_chim.fasta",                  emit: fasta
      path "usearch_global.full.75.map.uc",        emit: uc
      path "usearch_global.full.75.blast6out.txt", emit: blast6out
      path "excluded.txt",                         emit: excluded, optional: true

    script:
    """
    echo -e "Chimera filtering\n"

    ## vsearch usearch_global
    echo -e "..Running VSEARCH\n"
    vsearch \
      --usearch_global ${input} \
      --db             ${db} \
      --strand         plus \
      --id             0.75 \
      --threads        ${task.cpus} \
      --uc             usearch_global.full.75.map.uc \
      --blast6out      usearch_global.full.75.blast6out.txt \
      --output_no_hits

    ## Handle all potentially chimeric sequences from usearch_global
    echo -e "\n..Processing VSEARCH output\n"
    exclude_chims.py \
      --global_infile usearch_global.full.75.blast6out.txt \
      --infile        ${input} \
      --outfile       seqs_out_chim.fasta \
      --log_file      err.log \
      --ex_file       excluded.txt \
      --region        ${params.its_region} \

    echo -e "..Done"

    """
}


// Additional quality controls
process exclude_non_iupac {

    label "main_container"

    // cpus 1

    input:
      path input
      path inputuniq

    output:
      path "iupac_out_full.fasta",       emit: iupac
      path "iupac_out_vsearch_96.fasta", emit: iupac96
      path "excluded.txt",               emit: excluded

    script:
    """
    echo -e "Additional quality controls\n"

    ## Additional quality controls
    ## Remove low quality sequences (too short or with too many non-IUPAC symbols)
    exclude_non_iupac.py \
      --infile              ${input} \
      --infile_orig         ${inputuniq} \
      --outfile             iupac_out_full.fasta \
      --outfile_vsearch_96  iupac_out_vsearch_96.fasta \
      --log_file            err.log \
      --ex_file             excluded.txt \
      --allowed_number      6

    echo -e "..Done"
    """
}


// Allow query sequences vary in length at 100% similarity
process seqlen_variation {

    label "main_container"

    cpus 8

    input:
      path input   // iupac_out_vsearch_96.fasta

    output:
      path "centroids_100.fasta", emit: seqs
      path "clusters_100.uc",     emit: uc

    script:
    """
    echo "Running vsearch 100% clustering"
    
    vsearch \
      --cluster_fast ${input} \
      --id           1 \
      --iddef        0 \
      --threads      ${task.cpus} \
      --uc           clusters_100.uc \
      --centroids    centroids_100.fasta \
      --query_cov    0.96 \
      --target_cov   0.96

    echo -e "..Done"

    """
}


// Prepare representative sequences
process select_representatives {

    label "main_container"

    // cpus 1

    input:
      path centroids   // centroids_100.fasta
      path iupac       // iupac_out_full.fasta

    output:
      path "iupac_out_vsearch.fasta", emit: fasta

    script:
    """
    echo "Printing out vsearch representatives"
    
    ## vsearch representatives (the sequence count diff. is 9.5% for vsearch 4%)
    select_vsearch_reps.py \
      --infile_centroids ${centroids} \
      --infile_iupac     ${iupac} \
      --outfile          iupac_out_vsearch.fasta \
      --log_file         err.log

    echo -e "..Done"

    """

}


    ## Additional quality controls - Remove low quality sequences (too short or with too many non-IUPAC symbols)
    python3 "$script_dir/exclude_non_iupac.py" "$run_id" 6



    echo -e "..Done"

    """
}

}


workflow {
}

