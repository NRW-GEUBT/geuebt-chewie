import os
import time
import json
import pandas as pd
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


def validate_input_param(path, schema):
    try:
        df = pd.read_csv(path, index_col="sample", sep="\t", engine="python")
        validate(df, schema=schema)
    except FileNotFoundError:
        path = os.path.join(workflow.basedir, "..", ".tests", "integration", path)
        df = pd.read_csv(path, index_col="sample", sep="\t", engine="python")
        validate(df, schema=schema)


def get_conda_prefix(wildcards):
    try:
        # snakemake < 8.0
        return workflow.conda_prefix
    except:
        # snakemake > 8
        return workflow.deployment_settings.conda_prefix


# Input functions ------------------------------------
def aggregate_json_call(wildcards):
    "Aggregate JSON outputs of chewie_call"
    checkpoint_output = checkpoints.chewie_call.get(**wildcards).output["jsons"]
    ids_map = glob_wildcards(
        os.path.join(checkpoint_output, "{sample}.chewiesnake.json")
    ).sample
    return expand(
        "allele_calling/cgmlst/json/{sample}.chewiesnake.json", sample=ids_map
    )


def aggregate_qc_pass(wildcards):
    "Aggregate sample sheets of samples passing QC"
    checkpoint_output = checkpoints.stage_profiles.get(**wildcards).output[
        "isolate_sheet_dir"
    ]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{qc_pass}.json")).qc_pass
    return expand("staging/isolate_sheets/{qc_pass}.json", qc_pass=ids_map)


def aggregate_trees(wildcards):
    "Aggregate Newick trees"
    checkpoint_output = checkpoints.merge_profiles.get(**wildcards).output[
        "profile_dir"
    ]
    ids_map = glob_wildcards(os.path.join(checkpoint_output, "{cluster_id}.tsv")).cluster_id
    return expand("trees/trees/{cluster_id}.tre", cluster_id=ids_map)


def get_setting_value(file, value):
    with open(file, "r") as fi:
        d = json.load(fi)
    return d[value]


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")
validate_input_param(config["sample_sheet"], schema="../schema/samples.schema.yaml")
