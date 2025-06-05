#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import os, glob
import requests
import pandas as pd
from urllib.parse import urljoin


FASTA_EXT = [".fa", ".faa", ".fna", ".ffn", ".frn", ".fasta"]

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


def list_scheme_files(folder, exts=FASTA_EXT):
    files = []
    for ext in exts:
        pathstring = os.path.join(folder, f"*{ext}")
        files.extend(glob.glob(pathstring))
    return files


def main(sample_list, profiles_out, statistics_out, timestamps_out, url, organism, scheme):
    if not USERNAME or not PASSWORD:
        raise RuntimeError("Missing API_USERNAME or API_PASSWORD env vars")
    token = login(url, USERNAME, PASSWORD)
    response_list = authenticated_request(
        "GET",
        urljoin(url, "isolates"),
        token,
        params={"species": organism}
    )

    if response_list.status_code != 200:
        raise ValueError(
            f"Failed to retrieve isolates for species {organism}"
            f"Status: {response_list.status_code}."
            f"Body: {response_list.text}"
        )

    profile_records, statistics_records = [], []
    isolate_ids, timestamps = [], []

    with open(sample_list, "r") as fi:
        local_ids = [id.strip() for id in fi.readlines()]

    for record in response_list.json():
        isolate_id = record["isolate_id"]

        # skip local samples, nesessary to distinguish to be able to update join clusters later
        if isolate_id in local_ids:
            continue

        response_profile = authenticated_request(
            "GET",
            urljoin(url, f"isolates/{isolate_id}/allele_profile"),
            token
        )

        # No profile yet - this should actually not happens but still should be caught if it does
        if response_profile.status_code == 404:
            continue
        
        jdata = response_profile.json()

        profile = jdata["profile"]
        massaged_profile = {entry["locus"]: str(entry["allele_crc32"]) for entry in profile}
        profile_records.append({"#FILE": isolate_id, **massaged_profile})

        statistics = jdata["statistics"]
        statistics["sample"] = isolate_id
        statistics_records.append(statistics)

        timestamp = jdata["updated_at"]
        isolate_ids.append(isolate_id)
        timestamps.append(timestamp)
    
    # Need to check if database contained entries.
    # If not need to get proper headers, otherwise chewie will crash
    if profile_records:
        d = pd.DataFrame.from_records(profile_records)
        d.to_csv(profiles_out, sep="\t", index=False)
    else:
        tab = "\t"  # backlash forbidden in f-strings
        fastanames = [os.path.split(p)[1] for p in list_scheme_files(scheme)]
        header = f"#FILE\t{tab.join(fastanames)}\n"
        with open(profiles_out, "w") as fo:
            fo.write(header)

    if statistics_records:
        d = pd.DataFrame.from_records(statistics_records)
        d.to_csv(statistics_out, sep="\t", index=False)
    else:
        header = "sample\tEXC\tINF\tLNF\tPLOT\tNIPH\tALM\tASM\n"
        with open(statistics_out, "w") as fo:
                fo.write(header)

    d = pd.DataFrame(data={"sample": isolate_ids, "date": timestamps})
    d.to_csv(timestamps_out, sep="\t", index=False)


if __name__ == '__main__':
    main(
        snakemake.input["sample_list"],
        snakemake.output['profiles'],
        snakemake.output['statistics'],
        snakemake.output['timestamps'],
        snakemake.params['url'],
        snakemake.params['organism'],
        snakemake.params['scheme'],
    )
