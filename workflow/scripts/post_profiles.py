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


def main(json_files, qc_out, isolate_sheet_dir, merged, sample_list, url, ver):
    os.makedirs(isolate_sheet_dir, exist_ok=True)
    mergedlist = []
    slist = []
    qc = {}

    for jsonfile in json_files:
        with open(jsonfile,"r") as fi:
            data = json.load(fi)
        
        isolate_id = data["run_metadata"]["sample_description"]["sample"]
        allele_profile = data["allele_profile"]
        allele_stats = {k: v for k,v in data["allele_stats"].items() if k in ["ALM","ASM", "EXC", "INF", "LNF", "NIPH", "PLOT"]}
        missing_frac = data["allele_stats"]["loci_missing_fraction"]
        jdata = {
            "cgmlst": {
                    "allele_profile": allele_profile,
                    "allele_stats": allele_stats
                },
            "qc_metrics": {"cgmlst_missing_fraction": missing_frac},
            "sample_info": {"geuebt_chewie_ver": ver}
        }

        if not USERNAME or not PASSWORD:
            raise RuntimeError("Missing API_USERNAME or API_PASSWORD env vars")
        token = login(url, USERNAME, PASSWORD)
        response = authenticated_request(
            "PUT",
            urljoin(url, f"isolates/{isolate_id}/allele_profile"),
            token,
            json=jdata
        )

        if response.status_code == 200:
            qc[isolate_id] = {"STATUS": "PASS", "MESSAGES": [response.json()["message"]]}

            jdata["isolate_id"] = isolate_id        

            mergedlist.append({isolate_id: jdata})
            slist.append(isolate_id)
            with open(os.path.join(isolate_sheet_dir, f"{isolate_id}.json"), "w") as fo:
                json.dump(jdata, fo, indent=4)
        
        # On errors
        elif response.status_code == 422: 
            err_type = response.json()["detail"][0]["type"]
            err_msg = response.json()["detail"][0]["msg"]
            err_field = "/".join(response.json()["detail"][0]["loc"])
            nice_error = f"VALIDATION ERROR '{err_type}': {err_msg}; for field: '{err_field}'"
            qc[isolate_id] = {"STATUS": "FAIL", "MESSAGES": [nice_error]}
        else:  # will catch 404 - no need for special case
            qc[isolate_id] = {"STATUS": "FAIL", "MESSAGES": [
                f"An unexpected error has occured, contact a geuebt admin."
                f"Status: {response.status_code}."
                f"Body: {response.text}"]
            }

    with open(merged, "w") as fo:
        json.dump(mergedlist, fo, indent=4)
    with open(qc_out, "w") as fo:
        json.dump(qc, fo, indent=4)
    # Dump join list
    with open(sample_list, 'w') as fi:
        fi.write("\n".join(slist))


if __name__ == '__main__':
    main(
        snakemake.input['jsons'],
        snakemake.output['qc'],
        snakemake.output['isolate_sheet_dir'],
        snakemake.output['merged'],
        snakemake.output['sample_list'],
        snakemake.params['url'],
        snakemake.params['ver'],
    )
