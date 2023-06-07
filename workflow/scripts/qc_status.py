#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json 


def main(jsons, max_missing_loci, json_path):
    dictout = {}
    for js in jsons:
        # Parse json
        with open(js, 'r') as fi:
            js_dict=json.load(fi)
        sample = js_dict["run_metadata"]["sample_description"]["sample"]
        missing = js_dict["allele_stats"]["loci_missing_fraction"]
        # Check fraction missing
        if missing <= max_missing_loci:
            status = "PASS"
            mssg = ['']
        else:
            status = "FAIL"
            mssg = [
                f"Too many missing core genome loci. "
                f"Fraction missing is {missing}, "
                f"maximum allowed is {max_missing_loci}."
            ]
        # Update dict
        dictout[sample] = {
            'STATUS': status,
            'MESSAGES': mssg
        }
    # Dump json
    with open(json_path, 'w') as fp:
        json.dump(dictout, fp, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['jsons'],
        snakemake.params['max_missing_loci'],
        snakemake.output['qc_status'],
    )
