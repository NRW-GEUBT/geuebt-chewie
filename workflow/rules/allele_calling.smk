# Use chewie to perform alle calling
# Produce QC status stage allele profiles


import os


checkpoint chewie_call:
    input:
        settings="dbdata/settings.json",
    output:
        outdir=directory("allele_calling"),
        jsons=directory("allele_calling/cgmlst/json"),
        profiles="allele_calling/cgmlst/allele_profiles.tsv",
        statistics="allele_calling/cgmlst/allele_statistics.tsv",
        timestamps="allele_calling/cgmlst/timestamps.tsv",
    params:
        chewie=os.path.expanduser(f"~/.nrw-geuebt/geuebt-chewie-{version}/chewieSnake/chewieSnake.py"),
        samples_sheet=config["sample_sheet"],
        max_threads_per_job=config["max_threads_per_job"],
        cgmlst_scheme=lambda w, input: get_setting_value(input.settings, "scheme_path"),
        prodigal=lambda w, input: get_setting_value(input.settings, "prodigal_path"),
        max_missing_loci=1,  # checks for max_missing_loci happen at the API level
        distance_method=lambda w, input: get_setting_value(input.settings, "distance_method"),
        clustering_method=lambda w, input: get_setting_value(input.settings, "clustering_method"),
        conda_prefix=get_conda_prefix,
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
