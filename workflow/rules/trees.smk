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
        if grapetree -p {input.profile} -m MSTreeV2 > {output.tree}; then
            echo "GrapeTree succeeded."
        else
            echo "GrapeTree failed. Generating fallback Newick tree."
            # Extract sample names and construct a simple Newick tree
            samples=$(cut -f1  {input.profile} | tail -n +2 | awk '{{printf "%s:0,", $1}}' | sed 's/,$//')
            echo "($samples);" > {output.tree}
        fi
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