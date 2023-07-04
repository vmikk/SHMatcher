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

// Create module aliases
// (to reuse the same process code, but chain the input-output of these processes)
include { clustering as clustering_97 } from './modules/clustering.nf'
include { clustering as clustering_95 } from './modules/clustering.nf'
include { clustering as clustering_90 } from './modules/clustering.nf'
include { clustering as clustering_80 } from './modules/clustering.nf'

/*
// Sequence clustering (97-95-90-80% clustering)
process clustering {
    tag "$threshold"
    label "main_container"

    // cpus 10

    input:
      path input
      val  threshold
      path iupac       // iupac_out_vsearch.fasta

    output:
      path "output_*.fasta",      emit: fasta
      path "clusters_*.uc",       emit: uc
      path "clusters_out_*.txt",  emit: txt

    script:
    def iupacc = iupac.name != 'NO_FILE' ? "--fasta ${iupac}" : "--fasta ${input}"

    """
    echo -e "Sequence clustering"
    echo -e "Similarity threshold:" ${threshold}

    ## Clustering
    echo -e "\n..Running USEARCH\n"
    usearch \
      -cluster_fast ${input} \
      -id           ${threshold} \
      -gapopen      0.0/0.0E \
      -gapext       1.0/0.5E \
      -sort         other \
      -uc           clusters_${threshold}.uc

    ## Parsing cluster information
    echo -e "\n..Parsing clusters\n"
    clusterparser_preclust_pre.py \
      --uc           clusters_${threshold}.uc \
      ${iupacc} \
      --clusters_txt clusters_out_${threshold}.txt \
      --output       output_${threshold}.fasta \
      --log_file     err.log

    echo -e "..Done"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch: \$(usearch -version | sed 's/usearch //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
*/


// Final round of sequence clustering (80%),
// Aggregation of clustering results
process clustering_final {

    label "main_container"

    // cpus 10

    input:
      path input       // in_80_pre.fasta
      path iupac       // iupac_out_vsearch.fasta
      path clust97     // clusters_out_97_pre.txt
      path clust95     // clusters_out_95_pre.txt
      path clust90     // clusters_out_90_pre.txt

    output:
      path "clusters/*",          emit: clusters
      path "singletons/*",        emit: singletons
      path "clusters_80.uc",      emit: uc80
      path "clusters_out_80.txt", emit: txt80

    script:

    """
    echo -e "Sequence clustering"
    echo -e "Similarity threshold: 80"

    ## Clustering
    echo -e "\n..Running USEARCH\n"
    usearch \
      -cluster_fast ${input} \
      -id           0.80 \
      -gapopen      0.0/0.0E \
      -gapext       1.0/0.5E \
      -sort         other \
      -uc           clusters_80.uc

    mkdir -p singletons
    mkdir -p clusters

    ## Parsing and aggregating clustering information
    echo -e "\n..Parsing and aggregating clustering info\n"
    clusterparser_preclust_final_pre.py \
      --cl80uc   clusters_80.uc  \
      --cl97     ${clust97} \
      --cl95     ${clust95} \
      --cl90     ${clust90} \
      --cl80     clusters_out_80.txt \
      --allseqs  ${iupac} \
      --log      err.log

    echo -e "..Done"
    """
}

}


//  Workflow
workflow {
      
  // Input file sequences (FASTA)
  ch_input = Channel.value(params.input)

  // Sequence preparaion
  seq_prep(ch_input)

}

}

