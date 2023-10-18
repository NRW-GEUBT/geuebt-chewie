# Join new samples to existing clusters


rule filter_samples:
    input:
        sample_list="common/sample_list.txt",
        profiles="allele_calling/cgmlst/allele_profiles.tsv",
        statistics="allele_calling/cgmlst/allele_statistics.tsv",
        timestamps="allele_calling/cgmlst/timestamps.tsv",
    output:
        new_profiles="qcfilter/allele_profiles.tsv",
        new_statistics="qcfilter/allele_statistics.tsv",
        new_timestamps="qcfilter/timestamps.tsv",
    message:
        "[Join clusters] Collecting samples passing QC"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_samples.log",
    shell:
        """
        exec 2> {log}
        # Profiles
        head -n1 {input.profiles} > {output.new_profiles}
        grep -w -f {input.sample_list} {input.profiles} >> {output.new_profiles}
        # stats
        head -n1 {input.statistics} > {output.new_statistics}
        grep -w -f {input.sample_list} {input.statistics} >> {output.new_statistics}
        # timestamps
        head -n1 {input.timestamps} > {output.new_timestamps}
        grep -w -f {input.sample_list} {input.timestamps} >> {output.new_timestamps}
        """


rule make_join_list:
    input:
        new_profiles="qcfilter/allele_profiles.tsv",
        new_statistics="qcfilter/allele_statistics.tsv",
        new_timestamps="qcfilter/timestamps.tsv",
    output:
        serovars="dummy/serovar_info.tsv",
        samples="qcfilter/sample_list.tsv",
    params:
        ext_profiles=config["external_profiles"],
        ext_timestamps=config["external_timestamps"],
        ext_statistics=config["external_statistics"],
    message:
        "[Join clusters] Collecting external clusters information"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/make_join_list.log",
    shell:
        """
        exec 2> {log}
        # Make realpaths
        ext_profiles=$(realpath {params.ext_profiles})
        ext_statistics=$(realpath {params.ext_statistics})
        ext_timestamps=$(realpath {params.ext_timestamps})
        new_profiles=$(realpath {input.new_profiles})
        new_statistics=$(realpath {input.new_statistics})
        new_timestamps=$(realpath {input.new_timestamps})
        # Write table
        echo "sender\tprofiles\tstatistics\ttimestamps" > {output.samples}
        echo "external\t$ext_profiles\t$ext_statistics\t$ext_timestamps" >> {output.samples}
        echo "new\t$new_profiles\t$new_statistics\t$new_timestamps" >> {output.samples}
        echo "sample\tcluster_name" > {output.serovars}
        """


# First join at main cluters level
rule chewie_join_main:
    input:
        samplelist="qcfilter/sample_list.tsv",
        serovars="dummy/serovar_info.tsv",
    output:
        outdir=directory("join_clusters/main/"),
        sample_cluster="join_clusters/main/merged_db/sample_cluster_information.tsv",
        orphans="join_clusters/main/merged_db/orphan_samples.tsv",
        distances="join_clusters/main/merged_db/distance_matrix.tsv",
    params:
        chewie=os.path.join(config["chewie_path"], "chewieSnake_join.py"),
        clustering_method=config["clustering_method"],
        distance_threshold=config["cluster_distance"],
        external_main_clusters=config["external_main_clusters"],
        distance_method=config["distance_method"],
        species_shortname=config["cluster_prefix"],
        conda_prefix={workflow.conda_prefix},
    message:
        "[Join clusters] Joining samples to precomputed subclusters with ChewieSnake-join"
    threads: workflow.cores
    conda:
        "../envs/chewie.yaml"
    log:
        "logs/chewie_join_main.log",
    shell:
        """
        exec 2> {log}

        python {params.chewie} \
            --sample_list {input.samplelist} \
            --working_directory {output.outdir} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --serovar_info {input.serovars} \
            --external_cluster_names {params.external_main_clusters} \
            --cluster \
            --distance_method {params.distance_method} \
            --species_shortname {params.species_shortname} \
            --use_conda \
            --condaprefix {params.conda_prefix} \
            --threads {threads}
        """


# Then join at subcluster level
rule chewie_join_sub:
    input:
        samplelist="qcfilter/sample_list.tsv",
        serovars="dummy/serovar_info.tsv",
    output:
        outdir=directory("join_clusters/sub/"),
        sample_cluster="join_clusters/sub/merged_db/sample_cluster_information.tsv",
        orphans="join_clusters/sub/merged_db/orphan_samples.tsv",
    params:
        chewie=os.path.join(config["chewie_path"], "chewieSnake_join.py"),
        clustering_method=config["clustering_method"],
        distance_threshold=config["subcluster_distance"],
        external_sub_clusters=config["external_sub_clusters"],
        distance_method=config["distance_method"],
        species_shortname=config["cluster_prefix"],
        organism=config["organism"],
        conda_prefix={workflow.conda_prefix},
    message:
        "[Join clusters] Joining samples to precomputed clusters with ChewieSnake-join"
    threads: workflow.cores
    conda:
        "../envs/chewie.yaml"
    log:
        "logs/chewie_join_sub.log",
    shell:
        """
        exec 2> {log}

        python {params.chewie} \
            --sample_list {input.samplelist} \
            --working_directory {output.outdir} \
            --clustering_method {params.clustering_method} \
            --distance_threshold {params.distance_threshold} \
            --serovar_info {input.serovars} \
            --external_cluster_names {params.external_sub_clusters} \
            --cluster \
            --distance_method {params.distance_method} \
            --species_shortname {params.species_shortname} \
            --use_conda \
            --condaprefix {params.conda_prefix} \
            --threads {threads}
        """


checkpoint stage_clusters:
    input:
        clusters_main="join_clusters/main/merged_db/sample_cluster_information.tsv",
        orphans_main="join_clusters/main/merged_db/orphan_samples.tsv",
        clusters_sub="join_clusters/sub/merged_db/sample_cluster_information.tsv",
        distances="join_clusters/main/merged_db/distance_matrix.tsv",
    output:
        dirout=directory("staging/clusters/"),
        merged="staging/clusters.json",
    params:
        prefix=config["cluster_prefix"],
        main_threshold=config["cluster_distance"],
        sub_threshold=config["subcluster_distance"],
        organism=config["organism"],
    message:
        "[Join clusters] Consolidating and staging clusters"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/stage_clusters.log",
    script:
        "../scripts/stage_clusters.py"
