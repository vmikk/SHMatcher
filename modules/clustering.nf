
// NB. for 97% clustering, two inputs are identical:
//   input = iupac_out_vsearch.fasta
//   iupac = iupac_out_vsearch.fasta
// But the same channel could not be used twice,
//  therefore "NO_FILE" will be passed as `iupac`,
//  and it is handled by `def iupacc = ...`

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
    ## NB. vsearch scores = 2 * usearch scores
    echo -e "\n..Running VSEARCH\n"
    vsearch \
        --cluster_smallmem ${input} \
        --id               ${threshold} \
        --gapopen          0I/0E \
        --gapext           2I/1E \
        --usersort \
        --uc               clusters_${threshold}.uc \
        --threads          ${task.cpus}

    # echo -e "..Running USEARCH"
    # usearch \
    #     -cluster_fast ${input} \
    #     -id           ${threshold} \
    #     -gapopen      0.0/0.0E
    #     -gapext       1.0/0.5E
    #     -sort         other \
    #     -uc           clusters_${threshold}.uc

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
