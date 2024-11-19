#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import json
import pandas as pd


def main(new_profiles, clusters, profile_dir, ext_profiles):
    os.makedirs(profile_dir, exist_ok=True)

    # load clusters
    with open(clusters, "r") as fi:
        cluslist = json.load(fi)

    # merge profiles
    new_prof = pd.read_csv(new_profiles, index_col=False, sep="\t", dtype=str)
    ext_prof = pd.read_csv(ext_profiles, index_col=False, sep="\t", dtype=str)
    profiles = pd.concat([new_prof, ext_prof], ignore_index=True)

    # map representative id to cluster id for each cluster
    id_mapping = {
        record.get("representative", "orphan_dummy"): record["cluster_id"] for record in cluslist
        }

    for cluster in cluslist:
        if cluster["cluster_number"] > 0:
            # get member profiles
            mem_prof = profiles.loc[profiles["#FILE"].isin(cluster["members"])]
            # output to file
            mem_prof.to_csv(
                os.path.join(profile_dir,  f"{cluster['cluster_id']}.tsv"),
                sep ="\t", index=False, na_rep="-",
            )
        else:
            # Special case orphan cluster
            all_members = cluster["members"] + list(id_mapping.values())
            # replace repr name in profile by cluster name
            profiles["#FILE"] = profiles["#FILE"].replace(id_mapping)
            # get member profiles
            mem_prof = profiles.loc[profiles["#FILE"].isin(all_members)]
            # output to file
            mem_prof.to_csv(
                os.path.join(profile_dir, f"{cluster['cluster_id']}.tsv"),
                sep ="\t", index=False, na_rep="-",
            )


if __name__ == '__main__':
    main(
        snakemake.input['new_profiles'],
        snakemake.input['clusters'],
        snakemake.output['profile_dir'],
        snakemake.params['ext_profiles'],
    )
