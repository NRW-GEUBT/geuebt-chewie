#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import json
import os


def main(jsons, qc, outdir):
    # Make output dir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Parse QC
    with open(qc, 'r') as fi:
        qc_status = json.load(fi)

    for js in jsons:
        dictout = {}
        # Parse JSON
        with open(js, 'r') as fi:
            jsdict = json.load(fi)
        sample = jsdict['run_metadata']['sample_description']['sample']
        # Make isolate sheet only if passing QC
        if qc_status['sample']['STATUS'] != 'PASS':
            continue
        # Reformat JSON
        dictout['qc_metrics'] = {
            'cgmlst_missing_fraction': jsdict['allele_stats']['loci_missing_fraction']
        }
        dictout['cgmlst'] = {
            'allele_profile': jsdict['allele_profile'],
            'allele_stats': {
                k: jsdict['allele_stats'][k] for k in ['EXC', 'INF', 'LNF', 'PLOT', 'NIPH', 'ALM' , 'ASM']
            }
        }
        # Dump in isolate sheet
        with open(os.path.join(outdir, sample), 'w') as fi:
            json.dump(dictout, fi, indent=4)


if __name__ == '__main__':
    main(
        snakemake.input['jsons'],
        snakemake.input['qc'],
        snakemake.output['isolate_sheet_dir'],
    )
