
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

