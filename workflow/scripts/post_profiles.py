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


def main(qc_in, json_files, qc_out, isolate_sheet_dir, merged, url):
    os.makedirs(isolate_sheet_dir, exist_ok=True)
    mergedlist = []

    with open(qc_in, "r") as fi:
        qc = json.load(fi)

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
                "allele_stats": allele_stats,
                "qc_metrics": {"cgmlst_missing_fraction": missing_frac}
            }
        }

        response = requests.put(
            urljoin(url, f"isolates/{isolate_id}/allele_profile"),
            json=jdata
        )

        if response.status_code == 200:
            qc[isolate_id]["MESSAGES"].append(response.json()["message"])

            jdata["isolate_id"] = isolate_id        

            merged.append({isolate_id: jdata})
            with open(os.path.join(isolate_sheet_dir, f"{isolate_id}.json"), "w") as fo:
                json.dump(jdata, fo, indent=4)
        
        # On errors
        elif response.status_code == 422: 
            err_type = response.json()["detail"][0]["type"]
            err_msg = response.json()["detail"][0]["msg"]
            err_field = "/".join(response.json()["detail"][0]["loc"])
            nice_error = f"VALIDATION ERROR '{err_type}': {err_msg}; for field: '{err_field}'"
            qc[isolate_id]["STATUS"] = "FAIL"
            qc[isolate_id]["MESSAGES"].append(nice_error)
        else:  # will catch 404 - no need for special case
            qc[isolate_id]["STATUS"] = "FAIL"
            qc[isolate_id]["MESSAGES"].append(
                f"An unexpected error has occured, contact a geuebt admin."
                f"Status: {response.status_code}."
                f"Body: {response.text}"
                )

    with open(merged, "w") as fo:
        json.dump(mergedlist, fo, indent=4)
    with open(qc_out, "w") as fo:
        json.dump(qc, fo, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['summary'],
        snakemake.output['sheet_out'],
        snakemake.output['merged'],
        snakemake.output['qc_out'],
        snakemake.params['API_url'],
    )
