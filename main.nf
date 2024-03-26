// Species hypothesis (SH) matching analysis tool
//  Based on the SH MATCHING analysis tool: https://github.com/TU-NHM/sh_matching_pub


params.input = false

params.its_region = "itsfull"  // alternatively, "its2"

params.db     = "sanger_refs_sh.udb"
params.dbfull = "sanger_refs_sh_full.udb"
params.shdata = "sh_matching/data"

// Sequence filtering by length
if(params.its_region == "itsfull" & params.minlen1 == null){
    params.minlen1 = 400
}
if(params.its_region == "its2" & params.minlen1 == null){
    params.minlen1 = 100
}

if(params.its_region == "itsfull" & params.minlen2 == null){
    params.minlen2 = 350
}
if(params.its_region == "its2" & params.minlen2 == null){
    params.minlen2 = 100
}

// Allow query sequences vary 4% in length at 100% similarity
params.seqlenvariation = false  // == `include_vsearch_step` param of sh_matching_pub

// Sequence filtering by alignment coverage
params.mincoverage = 85


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

// Sequence filtering
// putative `chimera` removal based on low similarity to the database
process chimera_filtering {

    label "main_container"
    // cpus 8

    publishDir "${params.outdir}/Excluded_sequences",
        pattern: "excluded.txt",
        saveAs: { filename -> "1_Chimerae.txt" }

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
        --minlen1       ${params.minlen1} \
        --minlen2       ${params.minlen2} \
        --mincoverage   ${params.mincoverage}

    echo -e "..Done"

    """
}


// Additional quality controls
process exclude_non_iupac {

    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/Excluded_sequences",
        pattern: "excluded.txt",
        saveAs: { filename -> "2_NonIUPAC.txt" }

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
    // cpus 8

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


// Skipping the vsearch 100% clustering step with 96% length coverage
process no_seqlen_variation {

    label "main_container"
    // cpus 8

    input:
        path input   // iupac_out_vsearch_96.fasta

    output:
        path "centroids_100.fasta", emit: seqs
        path "clusters_100.uc",     emit: uc

    script:
    """
    echo "Running vsearch fastx_uniques"

    vsearch \
        --fastx_uniques ${input} \
        --threads       ${task.cpus} \
        --uc            clusters_100.uc \
        --fastaout      centroids_100.fasta

    echo -e "..Done"

    """
}

// Prepare representative sequences
process select_representatives {

    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/Excluded_sequences",
        pattern: "excluded.txt",
        saveAs: { filename -> "3_NoRepresentative.txt" }

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

    // cpus 8

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
    // cpus 8

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
// Calculate 0.5% clusters (USEARCH calc_distmx & cluster_aggd, complete-linkage)
process agglomerative_clustering {

    label "main_container"
    tag "$input"
    // cpus 4

    input:
        path input    // Cluster0
        path iupac    // iupac_out_vsearch.fasta

    output:
        path "*_out_005",   emit: clusters
        // path "*_mx_005", emit: dist

    script:
    """
    echo -e "Usearch complete-linkage clustering"
    echo -e "Input: " ${input}

    ## Count number of sequences in a cluster
    NUMSEQS="\$( grep -c '^>' ${input} || : )"
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
        mkdir -p "calc_distm_out"

        ## Parse usearch clusters
        echo -e "\n..Parsing USEARCH clustering output\n"
        clusterparser_usearch_90_pre.py \
            --name            ${input} \
            --file            ${input}_clusters_2_90.uc \
            --tmp_file1       clusters_out_2_90_pre.txt \
            --tmp_file_nohits ${iupac} \
            --tmp_cl_file     tmp.txt \
            --tmp_singl_file  singletons.txt \
            --log_file        err.log

        ## Calculate usearch distance matrix
        echo -e "\n..Complete-linkage clustering of sub-clusters\n"
        calc_distm_formatter_90_pre.py \
            --cluster    ${input} \
            --uclust_dir clusters \
            --out_dir    calc_distm_out \
            --cl_tmp     tmp.txt \
            --threads    ${task.cpus}

        ## Concatenate clusters
        echo -e "\n..Concatenating sub-clusters\n"
        find calc_distm_out -name "Cluster*" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
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
        echo -e "\n..Complete-linkage clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_005 \
            -clusterout   ${input}_out_005 \
            -id 0.995 \
            -linkage max

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
    // cpus 8

    publishDir "${params.outdir}/sh_matching"

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


// Parse SH matching results, select best hits
process parse_sh {

    label "main_container"
    // cpus 1

    input:
        path iupac  // iupac_out_full.fasta
        path uc     // closedref.80.map.uc

    output:
        path "closedref.80-best-hits.map.uc", emit: matches
        path "hits.txt",     emit: hits
        path "hits.fasta",   emit: hits_seq
        path "nohits.fasta", emit: nohits_seq, optional: true

    script:
    """
    echo -e "Parsing results of SH matching\n"

    parse_usearch_results.py \
        --infile         ${iupac} \
        --map_file       ${uc} \
        --hits_fasta     hits.fasta \
        --nohits_fasta   nohits.fasta \
        --hits           hits.txt \
        --best_hits      closedref.80-best-hits.map.uc \
        --log            err.log

    echo -e "..Done"
    """
}


// HITS: Create compound clusters (!of core dataset only!)
// with user's sequences added to the existing data
process compound_clusters {

    label "main_container"
    // cpus 1

    input:
        path datadir // sh_matching/data
        path iupac   // iupac_out_full.fasta
        path uc      // closedref.80-best-hits.map.uc

    output:
        path "compounds/*.fas", emit: compounds
        path "compounds.txt",   emit: compounds_list   // compounds/tmp.txt

    script:
    """
    echo -e "Creating compound clusters of core dataset\n"

    mkdir -p compounds

    echo -e "..Preparing compounds\n"
    create_compound_clusters.py \
        --datadir  ${datadir} \
        --infile   ${iupac} \
        --uc       ${uc} \
        --log      err.log \

    echo -e "..Listing clusters\n"
    find compounds -type f -name "UCL10_*.fas" \
        | sed 's/compounds\\///' \
        | sort --version-sort \
        > compounds.txt

    echo -e "..Done"
    """
}


// Go through compound clusters and
// run 97% usearch clustering if needed (if >16000 in cluster size)
// -> calc 3.0% distance matrix to form SHs based on these
process clustering_compounds {

    label "main_container"
    // cpus 3

    // tag "$input"
    tag "${input.getName().replaceFirst(/.fas$/, "")}"

    input:
        path input    // UCL10_*.fas
        path iupac    // iupac_out_vsearch.fasta

    output:
        path "calc_distm_out/*_out_*", emit: clusters

    script:
    """
    echo -e "Clustering of compound clusters\n"
    echo -e "Input: " ${input}

    ## Count number of sequences in a cluster
    NUMSEQS="\$( grep -c '^>' ${input} || : )"
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
        clusterparser_usearch_90.py \
            --name            ${input} \
            --file            ${input}_clusters_2_90.uc \
            --tmp_file1       clusters_out_2_90_pre.txt \
            --tmp_file_nohits ${iupac} \
            --tmp_cl_file     tmp.txt \
            --tmp_singl_file  singletons.txt \
            --log_file        err.log

        ## Remove clusters without user sequences
        # echo -e "..Removing redundant clusters"
        # grep -L '^>iid_' -r clusters   | parallel -j1 "rm {}"
        # grep -L '^>iid_' -r singletons | parallel -j1 "rm {}"
        #
        # ## Fix the cluster and singleton lists
        # find clusters -type f -name "Cluster*" \
        #   | sed 's/clusters\\///' \
        #   | sort --version-sort \
        #   > tmp.txt
        #
        # find singletons -type f -name "Singleton*" \
        #   | sed 's/singletons\\///' \
        #   | sort --version-sort \
        #   > singletons.txt

        ## Calculate SHs (max 3.0% distance)
        ## Calculate usearch distance matrix and generate (SH) clusters
        echo -e "\n..Calculating SHs, single-linkage clustering\n"
        calc_distm_formatter_90.py \
            --cluster    ${input} \
            --uclust_dir clusters \
            --out_dir    calc_distm_out \
            --cl_tmp     tmp.txt \
            --threads    ${task.cpus}

        ## Merge clusters
        echo -e "\n..Merging sub-clusters\n"
        mv calc_distm_out calc_distm_out_subcl
        mkdir -p calc_distm_out

        echo -e "...0.5%\n"
        find calc_distm_out_subcl -regex "calc_distm_out_subcl/Cluster[0-9]+_out_005\$" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
            > calc_distm_out/${input}_out_005

        echo -e "...1.0%\n"
        find calc_distm_out_subcl -regex "calc_distm_out_subcl/Cluster[0-9]+_out_01\$" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
            > calc_distm_out/${input}_out_01

        echo -e "...1.5%\n"
        find calc_distm_out_subcl -regex "calc_distm_out_subcl/Cluster[0-9]+_out_015\$" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
            > calc_distm_out/${input}_out_015

        echo -e "...2.0%\n"
        find calc_distm_out_subcl -regex "calc_distm_out_subcl/Cluster[0-9]+_out_02\$" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
            > calc_distm_out/${input}_out_02

        echo -e "...2.5%\n"
        find calc_distm_out_subcl -regex "calc_distm_out_subcl/Cluster[0-9]+_out_025\$" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
            > calc_distm_out/${input}_out_025

        echo -e "...3.0%\n"
        find calc_distm_out_subcl -regex "calc_distm_out_subcl/Cluster[0-9]+_out_03\$" \
            | sort --version-sort \
            | xargs awk -f \$(which "renumber_clusters.awk") \
            > calc_distm_out/${input}_out_03


    else

        echo "sm:"${input}":""\$NUMSEQS"\n

        mkdir -p "calc_distm_out"

        ## Calculate usearch distance matrix and generate (SH) clusters
        # calc_distm_formatter_80.py

        ## Generate a distance matrix
        echo -e "\n..Generating distance matrix\n"
        usearch \
            -calc_distmx ${input} \
            -tabbedout   ${input}_mx_03 \
            -maxdist     0.03 \
            -threads     ${task.cpus}

        ## Agglomerative clustering (single-linkage)
        ## 97.0%
        echo -e "\n..97.0% clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_03 \
            -clusterout   calc_distm_out/${input}_out_03 \
            -id      0.97 \
            -linkage "min"

        ## 97.5%
        echo -e "\n..97.5% clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_03 \
            -clusterout   calc_distm_out/${input}_out_025 \
            -id      0.975 \
            -linkage "min"

        ## 98.0%
        echo -e "\n..98.0% clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_03 \
            -clusterout   calc_distm_out/${input}_out_02 \
            -id      0.98 \
            -linkage "min"

        ## 98.5%
        echo -e "\n..98.5% clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_03 \
            -clusterout   calc_distm_out/${input}_out_015 \
            -id      0.985 \
            -linkage "min"

        ## 99.0%
        echo -e "\n..99.0% clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_03 \
            -clusterout   calc_distm_out/${input}_out_01 \
            -id      0.99 \
            -linkage "min"

        ## 99.5%
        echo -e "\n..99.5% clustering\n"
        usearch \
            -cluster_aggd ${input}_mx_03 \
            -clusterout   calc_distm_out/${input}_out_005 \
            -id      0.995 \
            -linkage "min"

    fi

    echo -e "..Done"
    """
}





// Create parsing module aliases
// (to reuse the same process code, but chain the input-output of these processes)
include { analyse_usearch_output as analyse_usearch_output_030 } from './modules/analyse_usearch_output.nf'
include { analyse_usearch_output as analyse_usearch_output_025 } from './modules/analyse_usearch_output.nf'
include { analyse_usearch_output as analyse_usearch_output_020 } from './modules/analyse_usearch_output.nf'
include { analyse_usearch_output as analyse_usearch_output_015 } from './modules/analyse_usearch_output.nf'
include { analyse_usearch_output as analyse_usearch_output_010 } from './modules/analyse_usearch_output.nf'
include { analyse_usearch_output as analyse_usearch_output_005 } from './modules/analyse_usearch_output.nf'

/*
// Parse usearch output
process analyse_usearch_output {

    label "main_container"
    // cpus 1

    input:
      val  threshold     // 03
      path compound_map  // data/sh2compound_mapping.txt
      path compounds     // compounds/tmp.txt
      path matches_prev
      path (small_clusters, stageAs: "compounds/calc_distm_out/*")   // compounds/calc_distm_out/*.fas_out_threshold
      path (large_clusters, stageAs: "compounds/*")                  // compounds/UCL10_000035.fas_folder/...
      path best_hits_uc  // closedref.80-best-hits.map.uc

    output:
      path "matches_${threshold}.txt", emit: matches
      path "compounds/calc_distm_out/*_out_${threshold}_bc", emit: bc

    script:
    """
    echo -e "Parsing usearch output"

    analyse_usearch_output.py \
      --threshold     ${threshold} \
      --sh2compound   ${compound_map} \
      --matches       matches_${threshold}.txt \
      --matches_prev  ${matches_prev} \
      --clusters      ${compounds} \
      --besthitsuc    ${best_hits_uc} \
      --log err.log

    echo -e "..Done"
    """
}
*/


// Parse matches files to output information about input sequences
// and their belonging to SHs on different thresholds
process parse_matches {

    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/Matches"

    input:
        path (matches, stageAs: "matches/*")   // "matches/matches_*.txt"
        path sh2compound_mapping  // "/sh_matching/data/sh2compound_mapping.txt"
        path shs_out              // "/sh_matching/data/shs_out.txt"
        path compounds_out        // "/sh_matching/data/compounds_out.txt"
        path centroid2sh_mappings // "/sh_matching/data/centroid2sh_mappings.txt"
        path source_names         // "source_names"
        path duplic_seqs          // "duplic_seqs.txt"
        path seq_mappings         // "seq_mappings.txt"

    output:
        path "matches/matches_out_*.csv", emit: shmatches

    script:
    """
    echo -e "Parsing SH matches\n"

    parse_matches.pl \
        --matches_dir          "matches" \
        --sh2compound_file     ${sh2compound_mapping} \
        --shs_file             ${shs_out} \
        --compound_file        ${compounds_out} \
        --centroid2sh_file     ${centroid2sh_mappings} \
        --accno_seqs_file      ${source_names} \
        --duplicate_seqs_file1 ${duplic_seqs} \
        --duplicate_seqs_file2 ${seq_mappings}

    echo -e "..Done"
    """
}


}



// Merge `parse_matches*.pl` output into one CSV file
process merge_matches {

    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/Matches"

    input:
        path (matches, stageAs: "matches/*")   // "matches/matches_out_*.csv"
        path uc                                // "closedref.80-best-hits.map.uc"

    output:
        path "matches_out_all.csv", emit: matchesall

    script:
    """
    echo -e "Meging matches\n"

    merge_matches.py \
        --matches  "matches" \
        --besthits ${uc} \
        --outfile  "matches_out_all.csv"

    echo -e "..Done"
    """
}

    echo -e "..Done"
    """
}


// Convert matches into HTML report
process parse_matches_html {
    tag "$threshold"
    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/Matches"

    input:
        tuple val(threshold), path (matches, stageAs: "matches/*")
        // val(threshold)                        // 005
        // path (matches, stageAs: "matches/*")  // matches_out_*.csv & matches_1_out_*.csv

    output:
        path "matches_out_*.html", emit: html

    script:
    """
    echo -e "Parsing matches for html output\n"

    parse_matches_html.py \
        --threshold  ${threshold} \
        --matchesdir matches \
        --outfile    matches_out_${threshold}.html

    echo -e "..Done"
    """
}


// Create Krona chart
process krona {
    tag "$threshold"
    label "main_container"
    // cpus 1

    publishDir "${params.outdir}/Krona", pattern: "*.html"

    input:
        tuple val(threshold), path (matches, stageAs: "matches/*")
        // val(threshold)                        // 005
        // path (matches, stageAs: "matches/*")  // matches_out_*.csv

    output:
        path "krona_*.html",   emit: kronahtml
        // path "krona_*.txt", emit: kronatxt

    script:
    """
    echo -e "Creating Krona charts\n"

    ## Create input for Krona chart
    shmatches2kronatext.py \
        --threshold  ${threshold} \
        --matchesdir matches \
        --outfile    krona_${threshold}.txt

    ## Export Krona charts
    ktImportText \
        -o krona_${threshold}.html \
        krona_${threshold}.txt

    echo -e "..Done"
    """
}



//  Workflow
workflow {

  // Input file sequences (FASTA)
  ch_input = Channel.value(params.input)

  // Databases
  ch_db  = Channel.value(params.db)       // sanger_refs_sh.udb
  ch_dbf = Channel.value(params.dbfull)   // sanger_refs_sh_full.udb
  ch_shd = Channel.value(params.shdata)   // sh_matching/data
  ch_map = Channel.fromPath( params.shdata + '/sh2compound_mapping.txt' ).first()

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

  if(params.seqlenvariation){
    // Allow query sequences vary 4% in length at 100% similarity  (`include_vsearch_step` == "yes")
    seqlen_variation(exclude_non_iupac.out.iupac96)
    ch_seqlen_seqs = seqlen_variation.out.seqs
    ch_seqlen_uc   = seqlen_variation.out.uc
  } else {
    // Skipping the vsearch 100% clustering step with 96% length coverage (`include_vsearch_step` == "no")
    no_seqlen_variation(exclude_non_iupac.out.iupac96)
    ch_seqlen_seqs = no_seqlen_variation.out.seqs
    ch_seqlen_uc   = no_seqlen_variation.out.uc
  }


  // Selecting representative sequences
  select_representatives(
    ch_seqlen_seqs,
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
    ch_seqlen_uc
    )

  // Channel with clusters
  ch_clusters = clustering_final.out.clusters.flatten()

  // Distance matrix estimation, agglomerative clustering
  agglomerative_clustering(
    ch_clusters,
    select_representatives.out.fasta
    )

  // Channel with clustering info
  ch_distm = agglomerative_clustering.out.clusters.collect()

  // Select representative sequences
  select_core_reps(
    clustering_final.out.clusters_list,
    clustering_final.out.singletons_list,
    select_representatives.out.fasta,
    ch_distm,
    clustering_final.out.singletons
    )

  // SH matching
  sh_matching(
    select_core_reps.out.corereps,
    ch_dbf
    )

  // Parse results, select best hits
  parse_sh(
    exclude_non_iupac.out.iupac,
    sh_matching.out.matches
    )

  // Create compound clusters of core dataset
  compound_clusters(
    ch_shd,
    exclude_non_iupac.out.iupac,
    parse_sh.out.matches
    )

  // Channel with compound clusters
  ch_compounds = compound_clusters.out.compounds.flatten()

  // Clustering compounds and distance estimation
  clustering_compounds(
    ch_compounds,
    select_representatives.out.fasta
    )

  // Channel with cluster memberships
  ch_compclusters = clustering_compounds.out.clusters.collect()


  //// Parse usearch output
  // 3.0%
  no_prev_match = file('None')
  analyse_usearch_output_030(
    "03",
    ch_map,
    compound_clusters.out.compounds_list,
    no_prev_match,
    ch_compclusters,
    parse_sh.out.matches
    )

  // 2.5%
  analyse_usearch_output_025(
    "025",
    ch_map,
    compound_clusters.out.compounds_list,
    analyse_usearch_output_030.out.matches,
    ch_compclusters,
    parse_sh.out.matches
    )

  // 2.0%
  analyse_usearch_output_020(
    "02",
    ch_map,
    compound_clusters.out.compounds_list,
    analyse_usearch_output_025.out.matches,
    ch_compclusters,
    parse_sh.out.matches
    )

  // 1.5%
  analyse_usearch_output_015(
    "015",
    ch_map,
    compound_clusters.out.compounds_list,
    analyse_usearch_output_020.out.matches,
    ch_compclusters,
    parse_sh.out.matches
    )

  // 1.0%
  analyse_usearch_output_010(
    "01",
    ch_map,
    compound_clusters.out.compounds_list,
    analyse_usearch_output_015.out.matches,
    ch_compclusters,
    parse_sh.out.matches
    )

  // 0.5%
  analyse_usearch_output_005(
    "005",
    ch_map,
    compound_clusters.out.compounds_list,
    analyse_usearch_output_010.out.matches,
    ch_compclusters,
    parse_sh.out.matches
    )


  // Parse matches
  ch_all_matches = analyse_usearch_output_030.out.matches.concat(
    analyse_usearch_output_025.out.matches,
    analyse_usearch_output_020.out.matches,
    analyse_usearch_output_015.out.matches,
    analyse_usearch_output_010.out.matches,
    analyse_usearch_output_005.out.matches
    ).collect()

  ch_map = Channel.fromPath( params.shdata + '/sh2compound_mapping.txt' ).first()
  ch_shs = Channel.fromPath( params.shdata + '/shs_out.txt' ).first()
  ch_cmp = Channel.fromPath( params.shdata + '/compounds_out.txt' ).first()
  ch_cnt = Channel.fromPath( params.shdata + '/centroid2sh_mappings.txt' ).first()

  parse_matches(
    ch_all_matches,
    ch_map,                           // "data/sh2compound_mapping.txt"
    ch_shs,                           // "data/shs_out.txt"
    ch_cmp,                           // "data/compounds_out.txt"
    ch_cnt,                           // "data/centroid2sh_mappings.txt"
    seq_prep.out.namesorig,           // "source_names"
    clustering_final.out.duplicates,  // "duplic_seqs.txt"
    select_core_reps.out.seqmappings  // "seq_mappings.txt"
    )


  // Channel with matches
  ch_all_matches = parse_matches.out.shmatches
  // TODO - add `matches_1_out_*` from parse_matches_1 [ NOHITS part ]
  // ch_all_matches = parse_matches.out.shmatches + nohit_clusters...

  // Merge matches
  merge_matches(
    ch_all_matches,         // "matches/matches_out_{thld}.csv"
    parse_sh.out.matches    // "closedref.80-best-hits.map.uc"
    )

  // Channel with threshold values
  ch_thresholds = Channel.fromList( ['005', '01', '015', '02', '025', '03'] )

  // Combine thresholds and matches into tuples
  // (to reuse value of `ch_all_matches`)
  ch_threshold_match = ch_thresholds.combine(ch_all_matches.toList())

    // ch_threshold_match.view()
    //  [005, [matches_out_005.csv, matches_out_01.csv, matches_out_015.csv, matches_out_02.csv, matches_out_025.csv, matches_out_03.csv]]
    //  [01,  [matches_out_005.csv, matches_out_01.csv, matches_out_015.csv, matches_out_02.csv, matches_out_025.csv, matches_out_03.csv]]
    //  [015, [matches_out_005.csv, matches_out_01.csv, matches_out_015.csv, matches_out_02.csv, matches_out_025.csv, matches_out_03.csv]]
    //  [02,  [matches_out_005.csv, matches_out_01.csv, matches_out_015.csv, matches_out_02.csv, matches_out_025.csv, matches_out_03.csv]]
    //  [025, [matches_out_005.csv, matches_out_01.csv, matches_out_015.csv, matches_out_02.csv, matches_out_025.csv, matches_out_03.csv]]
    //  [03,  [matches_out_005.csv, matches_out_01.csv, matches_out_015.csv, matches_out_02.csv, matches_out_025.csv, matches_out_03.csv]]

  // Convert matches into HTML report
  parse_matches_html(ch_threshold_match)

  // Create Krona charts
  krona(ch_threshold_match)

}


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

