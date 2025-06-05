#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import requests
import pandas as pd
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


def authenticated_request(method, endpoint, token, **kwargs):
    headers = kwargs.pop("headers", {})
    headers["Authorization"] = f"Bearer {token}"
    return requests.request(method, endpoint, headers=headers, **kwargs)


def main(main_out, sub_out, url, organism):
    if not USERNAME or not PASSWORD:
        raise RuntimeError("Missing API_USERNAME or API_PASSWORD env vars")
    token = login(USERNAME, PASSWORD)
    response_list = authenticated_request(
        "GET",
        urljoin(url, "clusters"),
        token,
        params={"species": organism}
    )

    main_ids, main_repr = [], []
    sub_ids, sub_repr = [], []

    if response_list.status_code != 200:
        raise ValueError(
            f"Failed to retrieve clusters for species {organism}"
            f"Status: {response_list.status_code}."
            f"Body: {response_list.text}"
        )

    for record in response_list.json():
        cluster_id = record["cluster_id"]
        response_cluster = authenticated_request(
            "GET",
            urljoin(url, f"clusters/{cluster_id}"),
            token
        )
        jdata = response_cluster.json()
        main_ids.append(jdata["cluster_id"])
        main_repr.append(jdata["representative"])

        for subcluster in jdata.get("subclusters", []):
            sub_ids.append(subcluster["subcluster_id"])
            sub_repr.append(subcluster["representative"])
    
    d = pd.DataFrame(data={"sample": main_repr, "cluster_name": main_ids})
    d.to_csv(main_out, sep="\t", index=False)

    d = pd.DataFrame(data={"sample": sub_repr, "cluster_name": sub_ids})
    d.to_csv(sub_out, sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.output['main'],
        snakemake.output['sub'],
        snakemake.params['url'],
        snakemake.params['organism'],
    )
