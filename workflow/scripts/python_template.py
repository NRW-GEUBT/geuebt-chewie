#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import pandas as pd


def main(input, params, output):
    pass


if __name__ == '__main__':
    main(
        snakemake.input[''],
        snakemake.params[''],
        snakemake.output[''],
    )
