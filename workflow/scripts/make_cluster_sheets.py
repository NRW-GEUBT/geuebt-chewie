#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Handles cluster and subcluster name assignment and merging.
Outputs cluster sheet JSON files, with special treatment for merged clusters.
"""

import sys


# Redirect stderr to Snakemake log if running in pipeline
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
import numpy as np
import pandas as pd


def load_matrix(filepath):
    """
    Load a distance matrix from a file.
    Returns a tuple of (index list, numpy array).
    """
    df = pd.read_csv(filepath, header=None, index_col=0, sep=" ")
    return (df.index.to_list(), df.to_numpy())


def filter_matrix(ids, matrix, isolates):
    """
    Filter a distance matrix to only include specified isolates.
    Args:
        ids: List of all sample IDs in the matrix.
        matrix: Numpy array representing the distance matrix.
        isolates: List of sample IDs to filter by.
    Returns:
        A DataFrame containing the distances between the specified isolates.
    """
    # Get indices of samples
    indices = [ids.index(id) for id in isolates]
    # filter array with indices
    distances = matrix[np.ix_(indices, indices)]
    return pd.DataFrame(
        data=np.int_(distances),
        index=isolates,
        columns=isolates
    )


def annnotate_clusters(cluster_group):
    """
    Annotates a cluster row with its status and changes.
    status can be new or kept
    changes can be none, increase or merge
    """
    merged_from = None
    # Check status
    sender = set(cluster_group["Sender"])
    if sender == {"new"}:
        status, changes = "new", "increase"
    elif sender == {"kept"}:
        status, changes = "kept", "none"
    else:
        # If there are both new and kept samples, check if merged or just increase
        if ";" in cluster_group["external_cluster_name"].iloc[0]:
            # If external cluster name contains a semicolon, it means it was merged
            status, changes = "new", "merge"
            merged_from = cluster_group["external_cluster_name"].iloc[0]
        else:
            status, changes = "kept", "increase"
    # Create a new row with the status and changes
    return pd.Series({
        "representative_sample": cluster_group["representative_sample"].iloc[0],
        "external_cluster_name": cluster_group["external_cluster_name"].iloc[0],
        "status": status,
        "changes": changes,
        "merged_from": merged_from,
    })


def rename_clusters(df, prefix, sep="-"):
    """
    parse external_cluster_name to find existing numbering
    either keep the name if status is kept or assign a new name by incrementating the number
    names are in the form "prefix-number"
    """
    # Find the maximum number in the external_cluster_name field
    max_number = 0
    for name_field in set(df["external_cluster_name"]):
        if pd.isna(name_field):
            # If the name field is NaN, skip it
            continue
        for name in name_field.split(";"):
            number  = name.split(sep)[-1]
            try:
                number = int(number)
                if number > max_number:
                    max_number = number
            except ValueError:
                # If the name does not contain a number, skip it
                continue
    # Assign new names based on the status and changes
    for row in df.itertuples():
        if row.status == "new":
            # If the status is new, assign a new name with incremented number
            df.at[row.Index, "new_name"] = f"{prefix}{sep}{max_number + 1}"
            max_number += 1
        else:
            # If the status is kept, keep the external_cluster_name
            df.at[row.Index, "new_name"] = row.external_cluster_name
    return df


def main(
    cluster_info,
    orphans_info,
    subcluster_info,
    distances,
    prefix,
    main_dist,
    sub_dist,
    organism,
    cluster_json,
    merged_json
):
    # Load distance matrix
    distance_index, distance_array = load_matrix(distances)
    # Load clusters and subclusters infos
    clusters = pd.read_csv(
        cluster_info,
        sep="\t",
        index_col=0,
        usecols=[
            "sample",
            "cluster_name",
            "Sender",
            "representative_sample",
            "external_cluster_name",
            "representative_sample_externalcluster"
        ]
    )
    subclusters = pd.read_csv(
        subcluster_info,
        sep="\t",
        index_col=0,
        usecols=[
            "sample",
            "cluster_name",
            "Sender",
            "representative_sample",
            "external_cluster_name",
            "representative_sample_externalcluster"
        ]
    )
    # Generate a new df describing each main cluster with new or kept, and changes with none, increase or merge
    cluster_status = clusters.fillna(
        ""
    ).set_index(
        'cluster_name'
    ).groupby(
        "cluster_name"
    ).apply(
        annnotate_clusters
    ).reset_index(
    ).pipe(
        rename_clusters,
        prefix=prefix,
        sep="-"
    )
    cluster_json_list, merged_json_list = [], []
    # Generate a cluster sheet, getting the relevant sample IDs from the clusters df
    for cluster_row in cluster_status.itertuples():
        sample_ids = list(set(
            clusters[
                clusters["cluster_name"] == cluster_row.cluster_name
            ].index
        ))
        dist = filter_matrix(
            distance_index,
            distance_array,
            sample_ids
        )
        # Get subclusters for this cluster
        # Basically same workflow as for clusters but nested for only the cluster
        # If no subcluster subcluster_status is empty
        current_subcluster = subclusters[
            subclusters["cluster_name"] == cluster_row.cluster_name
        ]
        subcluster_status = current_subcluster.fill_na(
            ""
        ).set_index(
            'cluster_name'
        ).groupby(
            "cluster_name"
        ).apply(
            annnotate_clusters
        ).reset_index(
        ).pipe(
            rename_clusters,
            prefix=cluster_row.new_name,
            sep="."
        )
        # create a dict for the cluster sheet
        cluster_sheet = {
            "cluster_id": cluster_row.new_name,
            "cluster_number": int(cluster_row.new_name.split("-")[-1]),
            "organism": organism,
            "size": len(sample_ids),
            "representative": cluster_row.representative_sample,
            "AD_threshold": main_dist,
            "members": sample_ids,
            "root_members": list(set(sample_ids) - set(current_subcluster.index.to_list())),
            "subclusters": [
                {
                    "subcluster_id": sub_row.new_name,
                    "subcluster_number": int(sub_row.new_name.split(".")[-1]),
                    "size": len(current_subcluster[current_subcluster["cluster_name"] == sub_row.cluster_name].index.to_list()),
                    "representative": sub_row.representative_sample,
                    "AD_threshold": sub_dist,
                    "members": current_subcluster[current_subcluster["cluster_name"] == sub_row.cluster_name].index.to_list(),
                } for sub_row in subcluster_status.itertuples()
            ],
            "distance_matrix": json.loads(dist.to_json(orient="records")),
        }
        # Append to the list for JSON output
        if cluster_row.changes == "merge":
            merged_json_list.append({
                "merged_from": cluster_row.merged_from,
                "new_cluster":cluster_sheet
            })
        else:
            cluster_json_list.append(cluster_sheet)
    # Take care of orphans
    orphans = pd.read_csv(
        orphans_info,
        sep="\t",
        index_col=0
    )
    if orphans.empty:
        members = []
    else:
        members = orphans.index.to_list()
    # Get list of represensentative samples for main clusters
    repr_samples = cluster_status["representative_sample"].to_list()
    cluster_names = cluster_status["new_name"].to_list()
    # Distance_index rename repr samples to the cluster names
    for sample_id, cluster_id in zip(repr_samples, cluster_names):
        distance_index[distance_index.index(sample_id)] = cluster_id
    # If there are no orphans, then the distance matrix is empty
    dist = filter_matrix(
            distance_index,
            distance_array,
            members + cluster_names
        )
    # Create a dict for the orphans cluster sheet
    orphans_sheet = {
        "cluster_id": f"{prefix}-orphans",
        "cluster_number": 0,
        "organism": organism,
        "size": len(members),
        "representative": None,
        "AD_threshold": main_dist,
        "members": members,
        "root_members": members,
        "subclusters": [],
        "distance_matrix": json.loads(dist.to_json(orient="records")),
    }
    cluster_json_list.append(orphans_sheet)
    # Output JSON files
    with open(cluster_json, "w") as f:
        json.dump(cluster_json_list, f, indent=4)
    with open(merged_json, "w") as f:
        json.dump(merged_json_list, f, indent=4)


if __name__ == '__main__':
    main(
        cluster_info=snakemake.input["clusters_main"],
        orphans_info=snakemake.input["orphans_main"],
        subcluster_info=snakemake.input["clusters_sub"],
        distances=snakemake.input["distances"],
        prefix=snakemake.params["prefix"],
        main_dist=snakemake.params["main_threshold"],
        sub_dist=snakemake.params["sub_threshold"],
        organism= snakemake.params["organism"],
        cluster_json=snakemake.output["cluster_json"],
        merged_json=snakemake.output["merged_json"],
    )
