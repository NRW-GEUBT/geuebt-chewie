import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
validate(config, schema="../schema/config.schema.yaml")


# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


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
    "Aggregte sample sheets of samples passing QC"
    checkpoint_output = checkpoints.stage_profiles.get(**wildcards).output[
        "isolate_sheet_dir"
    ]
    ids_map = glob_wildcards(
        os.path.join(checkpoint_output, "{qc_pass}.json")).qc_pass
    return expand("staging/isolate_sheets/{qc_pass}.json", qc_pass=ids_map)
