# Use chewie to perform alle calling
# Produce QC status stage allele profiles


checkpoint chewie_call:
    output:
        outdir=directory("allele_calling"),
        jsons=directory("allele_calling/cgmlst/json"),
        profiles="allele_calling/cgmlst/allele_profiles.tsv",
        statistics="allele_calling/cgmlst/allele_statistics.tsv",
        timestamps="allele_calling/cgmlst/timestamps.tsv",
    params:
        chewie=os.path.join(config["chewie_path"], "chewieSnake.py"),
        samples_sheet=config["sample_sheet"],
        max_threads_per_job=config["max_threads_per_job"],
        cgmlst_scheme=config["cgmlst_scheme"],
        prodigal=config["prodigal"],
        max_missing_loci=config["max_missing_loci"],
        distance_method=config["distance_method"],
        clustering_method=config["clustering_method"],
        conda_prefix={workflow.conda_prefix},
    message:
        "[Allele calling] Calling alleles using ChewieSnake"
    threads: workflow.cores
    conda:
        "../envs/chewie.yaml"
    log:
        "logs/chewie_call.log",
    shell:
        """
        exec 2> {log}
        
        python {params.chewie} \
            --sample_list {params.samples_sheet} \
            --working_directory {output.outdir} \
            --scheme {params.cgmlst_scheme} \
            --prodigal {params.prodigal} \
            --max_fraction_missing_loci {params.max_missing_loci} \
            --distance_method {params.distance_method} \
            --clustering_method {params.clustering_method} \
            --noreport \
            --use_conda \
            --condaprefix {params.conda_prefix} \
            --threads {threads} \
            --threads_sample {params.max_threads_per_job}
        """


rule qc_status:
    input:
        jsons=aggregate_json_call,
    output:
        qc_status="staging/qc_status.json",
    params:
        max_missing_loci=config["max_missing_loci"],
    message:
        "[Allele calling] Checking allele profiles quality"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/qc_status.log",
    script:
        "../scripts/qc_status.py"


checkpoint stage_profiles:
    input:
        jsons=aggregate_json_call,
        qc="staging/qc_status.json",
    output:
        isolate_sheet_dir=directory("staging/isolate_sheets"),
    message:
        "[Allele calling] Writting isolate sheets"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/stage_profiles.log",
    script:
        "../scripts/stage_profiles.py"
