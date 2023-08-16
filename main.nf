// Species hypothesis (SH) matching analysis tool
//  Based on the SH MATCHING analysis tool: https://github.com/TU-NHM/sh_matching_pub


params.input = false

params.its_region = "itsfull"  // alternatively, "its2"

params.db     = "sanger_refs_sh.udb"
params.dbfull = "sanger_refs_sh_full.udb"
params.shdata = "sh_matching/data"



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
      --iddef        2 \
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
      path "excluded.txt", emit: excluded, optional: true

    script:
    """
    echo "Printing out vsearch representatives"
    
    ## vsearch representatives (the sequence count diff. is 9.5% for vsearch 4%)
    select_vsearch_reps.py \
      --infile_centroids ${centroids} \
      --infile_iupac     ${iupac} \
      --outfile          iupac_out_vsearch.fasta \
      --excluded         excluded.txt \
      --log_file         err.log

    echo -e "..Done"

    """

}

// Create clustering module aliases
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
      path fastanames  // source_runID_fastanames  // source_fasta_names
      path clust100uc  // clusters_100.uc

    output:
      path "clusters/*",          emit: clusters,   optional:true
      path "singletons/*",        emit: singletons, optional:true
      path "clusters_80.uc",      emit: uc80
      path "clusters_out_80.txt", emit: txt80
      path "clusters.txt",        emit: clusters_list
      path "singletons.txt",      emit: singletons_list
      path "duplic_seqs.txt",     emit: duplicates

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

    ## List clusters and singletons
    echo -e "\n..Listing clusters and singletons\n"
    
    find clusters -type f -name "Cluster*" \
      | sed 's/clusters\\///' \
      | sort --version-sort \
      > clusters.txt

    find singletons -type f -name "Singleton*" \
      | sed 's/singletons\\///' \
      | sort --version-sort \
      > singletons.txt

    ## Write vsearch clustering duplicates into `duplic_seqs.txt` file
    echo -e "\n..Parsing results\n"
    usearch_parser.py \
      --clusters    clusters.txt \
      --singletons  singletons.txt \
      --cov100_uniq ${fastanames} \
      --cov96_uniq  ${clust100uc} \
      --duplicates  duplic_seqs.txt \
      --log err.log

    echo -e "..Done"
    """
}



// Go through 80% uclust clusters and run 97% usearch clustering if needed (if >16000 in cluster size)
// Calculate 0.5% clusters (USEARCH calc_distmx & cluster_aggd)
process agglomerative_clustering {

    label "main_container"
    tag "$input"
    cpus 1

    input:
      path input    // Cluster0
      path iupac    // iupac_out_vsearch.fasta

    output:
      path "*_out_005",   emit: clusters
      // path "*_mx_005", emit: dist

    script:
    """
    echo -e "Usearch single-linkage clustering"
    echo -e "Input: " ${input}

    ## Count number of sequences in a cluster
    NUMSEQS="\$(grep -c '>' ${input})"
    echo -e "Number of sequences in the cluster: " "\$NUMSEQS"\n

    if (( "\$NUMSEQS" > "16000" ))
    then
       
       echo "to be split:"${input}":""\$NUMSEQS"\n

       ## Clustering
       echo -e "\n..97% clustering\n"
       usearch \
         -cluster_fast ${input} \
         -id 0.97 \
         -gapopen 0.0/0.0E \
         -gapext 1.0/0.5E \
         -sort other \
         -uc ${input}_clusters_2_90.uc \
         -threads ${task.cpus}
       
       mkdir -p "clusters"
       mkdir -p "singletons"
       mkdir -p "calc_distm_out"

       ## Parse usearch clusters
       echo -e "\n..Parsing USEARCH clustering output\n"
       clusterparser_usearch_90_pre.py \
         --name ${input} \
         --file ${input}_clusters_2_90.uc \
         --tmp_file1       clusters_out_2_90_pre.txt \
         --tmp_file_nohits ${iupac} \
         --tmp_cl_file     tmp.txt \
         --tmp_singl_file  singletons.txt \
         --log_file        err.log

       ## Calculate usearch distance matrix
       echo -e "\n..Single-linkage clustering\n"
       calc_distm_formatter_90_pre.py \
         --cluster    ${input} \
         --uclust_dir clusters \
         --out_dir    calc_distm_out \
         --cl_tmp     tmp.txt \
         --threads    ${task.cpus}

      ## Concatenate clusters
      find calc_distm_out -name "Cluster*" -exec cat {} + \
        > ${input}_out_005

    else

      echo "sm:"${input}":""\$NUMSEQS"\n

      ## Calculate usearch distance matrix
      # calc_distm_formatter_80_pre.py

      ## Generate a distance matrix
      echo -e "\n..Generating distance matrix\n"
      usearch \
        -calc_distmx ${input} \
        -tabbedout   ${input}_mx_005 \
        -maxdist     0.005 \
        -threads     ${task.cpus}

      ## Agglomerative clustering
      echo -e "\n..Single-linkage clustering\n"
      usearch \
        -cluster_aggd ${input}_mx_005 \
        -clusterout ${input}_out_005 \
        -id 0.995 \
        -linkage min

    fi

    echo -e "..Done"
    """
}

// Take 0.5% representatives as RepS, add USEARCH singletons
process select_core_reps {

    label "main_container"
    // cpus 8

    input:
      path clusters_list   // clusters.txt     // tmp.txt
      path singletons_list // singletons.txt
      path iupac           // iupac_out_vsearch.fasta
      path (distm, stageAs: "calc_distm_out/*")    // calc_distm_out/*_out_005
      path (singletons, stageAs: "singletons/*")   // "singletons/Singleton*"

    output:
      path "core_reps_pre.fasta", emit: corereps
      path "seq_mappings.txt",    emit: seqmappings

    script:
    """
    echo -e "Selecting core representative sequences"

    select_core_reps_usearch.py \
      --cl_list       ${clusters_list} \
      --singl_list    ${singletons_list} \
      --infile_iupac  ${iupac} \
      --reps_outfile  core_reps_pre.fasta \
      --seq_mappings  seq_mappings.txt \
      --log_file      err.log

    echo -e "..Done"
    """
}


// Find best matches to user’s sequences
// in the existing SH sequence dataset using `usearch_global` algorithm
process sh_matching {

    label "main_container"
    cpus 8

    input:
      path input  // core_reps_pre.fasta
      path db     // sanger_refs_sh_full.udb

    output:
      path "closedref.80.map.uc", emit: matches

    script:
    """
    echo -e "Matching core reps to existing SH sequence dataset\n"

    echo -e "..Running VSEARCH\n"
    vsearch \
      --usearch_global ${input} \
      --db ${db} \
      --strand plus \
      --id 0.8 \
      --threads ${task.cpus} \
      --iddef 0 \
      --gapopen 0I/0E \
      --gapext 2I/1E \
      --uc "closedref.80.map.uc" \
      --maxaccepts 3 \
      --maxrejects 0

    echo -e "..Done"
    """
}

}


//  Workflow
workflow {
      
  // Input file sequences (FASTA)
  ch_input = Channel.value(params.input)

  // Databases
  ch_db  = Channel.value(params.db)
  // Sequence preparaion
  seq_prep(ch_input)

  // Chimera filtering
  chimera_filtering(
    seq_prep.out.unique,
    ch_db)

  // Additional quality controls
  exclude_non_iupac(
    chimera_filtering.out.fasta,
    seq_prep.out.unique)

  // Allow query sequences vary 4% in length at 100% similarity
  seqlen_variation(exclude_non_iupac.out.iupac96)

  // Selecting representative sequences
  select_representatives(
    seqlen_variation.out.seqs,
    exclude_non_iupac.out.iupac)


  // 97% pre-clustering
  // NB. `iupac_out_vsearch.fasta` is used as input twice,
  //  but we can not use the same channel twice,
  //  therefore pass "NO_FILE" as the third input,
  //  and it should be handled by the process
  no_file = file('NO_FILE')
  clustering_97(
    select_representatives.out.fasta,
    0.97,
    no_file)

  // 95% pre-clustering
  clustering_95(
    clustering_97.out.fasta,
    0.95,
    select_representatives.out.fasta)

  // 90% pre-clustering
  clustering_90(
    clustering_95.out.fasta,
    0.90,
    select_representatives.out.fasta)

  // 80% pre-clustering
  clustering_final(
    clustering_90.out.fasta,
    select_representatives.out.fasta,
    clustering_97.out.txt,
    clustering_95.out.txt,
    clustering_90.out.txt,
    seq_prep.out.namesuniq,
    seqlen_variation.out.uc
    )

  // Channel with clusters
  ch_clusters = clustering_final.out.clusters.flatten()
// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

}

