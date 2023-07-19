#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import shutil


def is_empty(filepath):
    "basically checks if first line contians any info"
    with open(filepath, "r") as fi:
        if not fi.readline().strip():
            return True
    return False


def main(
    profiles_template,
    ext_main,
    ext_sub,
    ext_profiles,
    ext_timestamps,
    ext_statistics,
    ext_main_out,
    ext_sub_out,
    ext_profiles_out,
    ext_timestamps_out,
    ext_statistics_out
):
    # Check if input files are empty, if they are create a file with correct heading
    # correct heading may have to be inferred from chewie output eg profiles)
    for fpi, fpo, header in zip(
        [ext_main, ext_sub, ext_timestamps, ext_statistics],
        [ext_main_out, ext_sub_out, ext_timestamps_out, ext_statistics_out],
        [
            "sample\tcluster_name\n",
            "sample\tcluster_name\n",
            "sample\tdate\n",
            "sample\tEXC\tINF\tLNF\tPLOT\tNIPH\tALM\tASM\n",
        ]
    ):
        if is_empty(fpi):
            with open(fpo, "w") as fo:
                fo.write(header)
        else:
            shutil.copy(fpi, fpo)

    if is_empty(ext_profiles):
        with open(profiles_template, "r") as fi:
            header = fi.readline()
        with open(ext_profiles_out, "w") as fo:
            fo.write(header)
    else:
        shutil.copy(ext_profiles, ext_profiles_out)


if __name__ == '__main__':
    main(
        snakemake.input['new_profiles'],
        snakemake.params['external_main_clusters'],
        snakemake.params['external_sub_clusters'],
        snakemake.params['ext_profiles'],
        snakemake.params['ext_timestamps'],
        snakemake.params['ext_statistics'],
        snakemake.output['external_main_clusters'],
        snakemake.output['external_sub_clusters'],
        snakemake.output['ext_profiles'],
        snakemake.output['ext_timestamps'],
        snakemake.output['ext_statistics']
    )
