#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
from pathlib import Path


def main(tree_list, merged_path):
    json_out = {}
    
    for tree_file in tree_list:
        cluster_id = Path(tree_file).stem
        with open(tree_file, "r") as fi:
            json_out[cluster_id] = fi.readline().strip()
    
    with open(merged_path, "w") as fo:
        json.dump(json_out, fo, indent=4)        


if __name__ == '__main__':
    main(
        snakemake.input['trees'],
        snakemake.output['merged'],
    )
