checkpoint merge_profiles:
    input:
        new_profiles="qcfilter/allele_profiles.tsv",
        clusters="join_clusters/clusters.json",
        ext_profiles="dbdata/profiles_ext.tsv",
    output:
        profile_dir=directory("trees/merged_profiles")
    params:
        
    message:
        "[Trees] Joining allele profiles per clusters"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/merge_profiles.log"
    script:
        "../scripts/merge_profiles.py"


# Grapetree crashes on null matrix so catch exit 0 and manually make the tree
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