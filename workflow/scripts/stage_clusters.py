#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Assumes Cluster naming to be in the form:
PREFIX-<int main number>.<int sub number>
"""

import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os
import json
import numpy as np
import pandas as pd
from json import loads


def _load_matrix(filepath):
    df = pd.read_csv(filepath, header=None, index_col=0, sep="\t")
    return (df.index.to_list(), df.to_numpy())


def _get_cluster_numbers(x, sep):
    if not pd.isna(x["external_cluster_name"]):
        return x["external_cluster_name"].split(sep)[1]
    else:
        return pd.NA


def _assign_subcluster_names(x, mapping):
    if pd.isna(x['external_cluster_name']):
        return mapping[x['cluster_name']]
    else:
        return x['external_cluster_name']


def _get_subclusters_list(merged_df, prefix, dist):
    listout = []
    for clsname in merged_df['sub_name'].dropna().unique():
        current = merged_df.loc[merged_df["sub_name"] == clsname].reset_index()
        d = {
            "subcluster_id": f"{prefix}-{current.at[0, 'main_name']}.{current.at[0, 'sub_name']}",
            "subcluster_number": int(current.at[0, 'sub_name']),
            "size": len(current["sample"].to_list()),
            "representative": current.at[0, 'sub_repr'],
            "AD_threshold": dist,
            "members": current['sample'].to_list()
        }
        listout.append(d)
    return listout


def _filter_matrix(ids, matrix, isolates):
    # Get indices of samples
    indices = [ids.index(id) for id in isolates]
    # filter array with indices
    distances = matrix[np.ix_(indices, indices)]
    return pd.DataFrame(
        data=np.int_(distances),
        index=isolates,
        columns=isolates
    )


def map_names(cluster_df):
    """
    Assign available names to new clusters
    For simplicity assumes that namin is incremented numbering.
    While gaps in the numbering is allowed, string naming is not.
    Use an indexed list of names if you really want strings.
    """
    # get used cluster names
    used_names = list(cluster_df['external_cluster_name'].dropna().unique())
    # extract max cluster number
    if used_names:
        max_used = max([int(clsnum) for clsnum in used_names])
    else:
        max_used = 0
    # For each new cluster, get a new name as increments of the max name used
    new_clusters_df = cluster_df.loc[cluster_df["external_cluster_name"].isna()]
    cluster_name_list = list(new_clusters_df['cluster_name'].unique())
    name_map = {
        clust: cluster_name_list.index(clust) + 1 + max_used
        for clust in cluster_name_list
    }
    return name_map


def main(
    cluster_info,
    orphans,
    subcluster_info,
    distances,
    prefix,
    main_dist,
    sub_dist,
    organism,
    dirout,
    mergedout
):
    # Make output dir
    if not os.path.isdir(dirout):
        os.mkdir(dirout)

    # Load data
    clusters = pd.read_csv(cluster_info, index_col=False, sep="\t")
    subclusters = pd.read_csv(subcluster_info, index_col=False, sep="\t")
    orphs = pd.read_csv(orphans, index_col=False, sep="\t")
    dist = _load_matrix(distances)  # tuple(index, ndarray)

    # Parse cluster names and reformat
    clusters["external_cluster_name"] = clusters.apply(
        _get_cluster_numbers,
        axis=1,
        sep="-"
    )
    subclusters["external_cluster_name"] = subclusters.apply(
        _get_cluster_numbers,
        axis=1,
        sep="."
    )

    merged_json = []

    # Main cluster name assignment
    main_map = map_names(clusters)

    # collect cluster names and repr
    cluster_tracking = []  # list of tuples (cluster_repr, cluster_name)
    # Iterate over clusters - not as efficient as groupy but more manageable
    for cluster_tmpname in clusters['cluster_name'].unique():
        cluster_df = clusters.loc[clusters['cluster_name'] == cluster_tmpname].reset_index()
        # if no new samples, just add repr and cluster name to tracker and skip
        if "new" not in cluster_df["Sender"].to_list():
            cluster_tracking.append((
                cluster_df.at[0, "representative_sample_externalcluster"],
                f"{prefix}-{cluster_df.at[0, 'external_cluster_name']}"
            ))
            continue
        # Assign cluster name if empty
        if cluster_df["external_cluster_name"].isna().all():
            cluster_df["external_cluster_name"] = main_map[cluster_tmpname]
        cluster_short = cluster_df[
            ["sample", "Sender", "representative_sample", "external_cluster_name"]
        ]
        cluster_short = cluster_short.rename(columns={
            "external_cluster_name": "main_name",
            "representative_sample": "main_repr"
        })

        # Select samples in subclusters based on current cluster
        subcluster_df = subclusters.loc[subclusters['sample'].isin(cluster_df['sample'])].reset_index()
        # Assign new names to subclusters
        # If there are no subclusters, then all samples are root level in cluster
        if subcluster_df.empty:
            merged = cluster_short
            merged["sub_name"] = pd.NA
            merged["sub_repr"] = pd.NA
        else:
            sub_map = map_names(subcluster_df)
            subcluster_df["external_cluster_name"] = subcluster_df.apply(
                _assign_subcluster_names,
                axis=1,
                mapping=sub_map
            )
            subcluster_short = subcluster_df[
                ["sample", "representative_sample", "external_cluster_name"]
            ]
            subcluster_short = subcluster_short.rename(columns={
                "external_cluster_name": "sub_name",
                "representative_sample": "sub_repr"
            })
            # Join mains with sublcusters - orphans get NA
            merged = pd.merge(
                cluster_short,
                subcluster_short,
                how='left',
                on="sample"
            )

        # Format as dict for JSON, or skip if nothing new
        if merged.empty:
            raise ValueError("merged array should not be empty")
        distance_json = _filter_matrix(
            *dist, merged['sample'].to_list()
        ).to_json(orient="records")
        d = {
            "cluster_id": f"{prefix}-{merged.at[0, 'main_name']}",
            "cluster_number": int(merged.at[0, 'main_name']),
            "organism": organism,
            "size": len(merged['sample'].to_list()),
            "representative": merged.at[0, 'main_repr'],
            "AD_threshold": main_dist,
            "root_members": merged['sample'].loc[merged['sub_name'].isna()].to_list(),
            "subclusters": _get_subclusters_list(merged, prefix, sub_dist),
            "distance_matrix": loads(distance_json)
        }
        cluster_tracking.append((
            merged.at[0, 'main_repr'],
            f"{prefix}-{merged.at[0, 'main_name']}"
        ))
        # Dump JSONs
        with open(os.path.join(dirout, f"{d['cluster_id']}.json"), 'w') as fp:
            json.dump(d, fp, indent=4)

        merged_json.append(d)

    # If there are orphans, also make an orphan JSON
    # The orphan JSON will be used as global record including the representative
    # of each main clusters, and the distances.
    if orphs.empty:
        members = []
    else:
        members = orphs["sample"].to_list()
    # getting distance matrix of orphans + representatives of mainclusters
    representatives = [tpl[0] for tpl in cluster_tracking]
    distance_global = _filter_matrix(
        *dist, [*members, *representatives]
    )
    mapper = {old: new for old, new in cluster_tracking}
    # converting name of representatives to cluster names
    distance_json = distance_global.rename(columns=mapper).to_json(orient="records")
    d = {
        "cluster_id": f"{prefix}-orphans",
        "cluster_number": 0,
        "organism": organism,
        "size": len(members),
        "AD_threshold": main_dist,
        "root_members": members,
        "distance_matrix": loads(distance_json)
    }
    with open(os.path.join(dirout, f"{d['cluster_id']}.json"), 'w') as fp:
        json.dump(d, fp, indent=4)

    merged_json.append(d)

    with open(mergedout, 'w') as fp:
        json.dump(merged_json, fp, indent=4)


if __name__ == '__main__':
    main(
        cluster_info=snakemake.input['clusters_main'],
        orphans=snakemake.input['orphans_main'],
        subcluster_info=snakemake.input['clusters_sub'],
        distances=snakemake.input['distances'],
        prefix=snakemake.params['prefix'],
        main_dist=snakemake.params['main_threshold'],
        sub_dist=snakemake.params['sub_threshold'],
        organism=snakemake.params['organism'],
        dirout=snakemake.output['dirout'],
        mergedout=snakemake.output['merged']
    )
