gh# GET data


# POST data
rule post_profiles:
    input:
        qc="qc/qc_alleles.json",
        jsons=aggregate_json_call,
    output:
        qc="qc/qc_post_alleles.json",
        isolate_sheet_dir=directory("staging/isolate_sheets"),
        merged="staging/isolate_sheets.json",
    params:
        url=config["API_url"],
    message:
        "[Staging] POST profiles"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/post_profiles.log",
    script:
        "../scripts/post_profiles.py"


rule post_clusters:
    input:
        qc="qc/qc_post_alleles.json",
        trees="trees/trees.json",
        clusters="join_clusters/clusters.json",
    output:
        qc="staging/qc_status.json",
        custer_sheet_dir=directory("staging/clusters"),
        merged="staging/clusters.json"
    params:
        url=config["API_url"],
    message:
        "[Staging] POST clusters"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/post_clusters.log",
    script:
        "../scripts/post_clusters.py"
