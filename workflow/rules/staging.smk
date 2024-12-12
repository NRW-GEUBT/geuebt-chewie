# GET data
rule get_settings:
    output:
        settings="dbdata/settings.json",
    params:
        url=config["API_url"],
        organism=config["organism"],
    message:
        "[BakCharak] Getting analysis settings"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/get_settings.log",
    script:
        "../scripts/get_settings.py"


rule get_clusters_repr:
    input:
        settings="dbdata/settings.json",
    output:
        main="dbdata/clusters_ext.tsv",
        sub="dbdata/subclusters_ext.tsv",
    params:
        url=config["API_url"],
        organism=lambda w, input: get_setting_value(input.settings, "organism"),
    message:
        "[Staging] Getting cluster information from databse"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/get_clusters_repr.log",
    script:
        "../scripts/get_clusters.py"


rule get_isolate_info:
    input:
        settings="dbdata/settings.json",
        # Ensure that profile fetching happens after profile posting
        # and that for new samples, profiles are fetched locally
        sample_list="common/sample_list.txt",
    output:
        profiles="dbdata/profiles_ext.tsv",
        statistics="dbdata/statistics_ext.tsv",
        timestamps="dbdata/timestamps_ext.tsv",
    params:
        url=config["API_url"],
        organism=lambda w, input: get_setting_value(input.settings, "organism"),
        scheme=lambda w, input: get_setting_value(input.settings, "scheme_path"),
    message:
        "[Staging] Getting isolates information from database"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/get_isolate_info.log",
    script:
        "../scripts/get_isolates.py"


# POST data
rule post_profiles:
    input:
        jsons=aggregate_json_call,
    output:
        qc="staging/qc_status.json",
        isolate_sheet_dir=directory("staging/isolate_sheets"),
        merged="staging/isolate_sheets.json",
        sample_list="common/sample_list.txt",
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
        trees="trees/trees.json",
        clusters="join_clusters/clusters.json",
    output:
        qc="staging/cluster_response_status.json",
        cluster_sheet_dir=directory("staging/clusters"),
        merged="staging/clusters.json",
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
