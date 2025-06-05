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
import requests
from urllib.parse import urljoin


USERNAME = os.getenv("GEUEBT_API_USERNAME")
PASSWORD = os.getenv("GEUEBT_API_PASSWORD")


def login(url, username, password):
    response = requests.post(
        f"{url}/users/token",
        data={"username": username, "password": password},
        headers={"Content-Type": "application/x-www-form-urlencoded"}
    )
    response.raise_for_status()
    return response.json()["access_token"]


def authenticated_request( method, endpoint, token, **kwargs):
    headers = kwargs.pop("headers", {})
    headers["Authorization"] = f"Bearer {token}"
    return requests.request(method, endpoint, headers=headers, **kwargs)


def main(trees_in, clusters_in, qc_out, cluster_dir, merged_clusters, url):
    os.makedirs(cluster_dir, exist_ok=True)

    mergedlist = []
    qc = {}

    with open(clusters_in, "r") as fi:
        clusters = json.load(fi)
    with open(trees_in, "r") as fi:
        trees = json.load(fi)
    
    if not USERNAME or not PASSWORD:
        raise RuntimeError("Missing API_USERNAME or API_PASSWORD env vars")
    token = login(url, USERNAME, PASSWORD)

    for record in clusters:
        cluster_id = record["cluster_id"]
        record["tree"] = trees[cluster_id]

        response = authenticated_request(
            "PUT",
            urljoin(url, f"clusters/{cluster_id}"),
            token,
            json=record
        )

        if response.status_code == 200:
            qc[cluster_id] = {"STATUS": "PASS", "MESSAGES": [response.json()["message"]]}

            mergedlist.append(record)

            with open(os.path.join(cluster_dir, f"{cluster_id}.json"), "w") as fo:
                json.dump(record, fo, indent=4)
        
        # On errors
        elif response.status_code == 422: 
            err_type = response.json()["detail"][0]["type"]
            err_msg = response.json()["detail"][0]["msg"]
            err_field = "/".join(response.json()["detail"][0]["loc"])
            nice_error = f"VALIDATION ERROR '{err_type}': {err_msg}; for field: '{err_field}'"
            qc[cluster_id] = {"STATUS": "FAIL", "MESSAGES": [nice_error]}
        else:  # will catch 404 - no need for special case
            qc[cluster_id] = {"STATUS": "FAIL", "MESSAGES": [
                f"An unexpected error has occured, contact a geuebt admin."
                f"Status: {response.status_code}."
                f"Body: {response.text}"]
            }

    with open(merged_clusters, "w") as fo:
        json.dump(mergedlist, fo, indent=4)
    with open(qc_out, "w") as fo:
        json.dump(qc, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['trees'],
        snakemake.input['clusters'],
        snakemake.output['qc'],
        snakemake.output['cluster_sheet_dir'],
        snakemake.output['merged'],
        snakemake.params['url'],
    )
