#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


# if not calling for snakemake rule
try:
    sys.stderr = open(snakemake.log[0], "w")
except NameError:
    pass


import requests
import json
from urllib.parse import urljoin


def main(settings, url, organism):
    response = requests.get(
        urljoin(url, f"settings/{organism}")
    )

    if response.status_code != 200:
        raise ValueError(
            f"Failed to retrieve clusters for species {organism}"
            f"Status: {response.status_code}."
            f"Body: {response.text}"
        )


    with open(settings, "w") as fo:
        json.dump(response.json(), fo, indent = 4)


if __name__ == "__main__":
    main(
        snakemake.output['settings'],
        snakemake.params['url'],
        snakemake.params['organism'],
    )
