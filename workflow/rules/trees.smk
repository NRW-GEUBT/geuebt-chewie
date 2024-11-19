checkpoint merge_profiles:
    input:
        new_profiles="qcfilter/allele_profiles.tsv",
        clusters="join_clusters/clusters.json",
    output:
        profile_dir=directory("trees/merged_profiles")
    params:
        ext_profiles=config["external_profiles"],
    message:
        "[Trees] Joining allele profiles per clusters"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_profiles.log"
    script:
        "../scripts/merge_profiles.py"


rule grapetree:
    input:
        profile="trees/merged_profiles/{cluster_id}.tsv"
    output:
        tree="trees/trees/{cluster_id}.tre"
    message:
        "[Trees] Growing tree for {wildcards.cluster_id}"
    conda:
        "../envs/grapetree.yaml"
    log:
        "logs/{cluster_id}_grapetree.log"
    shell:
        """
        exec 2> {log}
        # grapetree crashes on 0 matrix, need to genereate single node tree manually
        # basically just "(id1,id2,id3);" NB: omitting distance is the same as distance 0
        grapetree -p {input.profile} -m MSTreeV2 > {output.tree} || echo \($(cut -f1 {input.profile} | tail -n+2 | paste -sd "," -)\)\; > {output.tree}
        """


rule aggregate_trees:
    input: 
        trees=aggregate_trees,
    output:
        merged="trees/trees.json",
    message:
        "[Trees] Aggregating trees"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/aggregate_trees.log"
    script:
        "../scripts/aggregate_trees.py"