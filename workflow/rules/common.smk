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
def aggregate_samples(wildcards):
    checkpoint_output = checkpoints.checkpoint_name.get(**wildcards).output[0]
    ids_map = glob_wildcards(
        os.path.join(checkpoint_output, "{wildcard_name}.fa")
    ).wildcard_name
    return expand("file/path/pattern/{wildcard_name}.ext", wildcard_name=ids_map)

