# Guide for users

## Installation

### Conda

Install conda from any distribution, i.e. miniconda.
You can follow the setup guide from the [Bioconda team](https://bioconda.github.io/).

We advise installing the mamaba solver in the base environement to speed up
environments creation.

```bash
conda install mamba -n base -c conda-forge
```

### Run environment

The running environement simply required a recent python version (>= 3.9) and snakemake.
If you followed the steps above just run:

```bash
mamba create -n snakemake snakemake
```

### Install module and databases

Download the [latest realease](https://github.com/NRW-GEUBT/geuebt-chewie/releases/latest)
and unpack it.

If you're feeling brave, clone the repository form Github:

```bash
git clone https://github.com/NRW-GEUBT/geuebt-chewie
```

Most software and databases dependencies will be installed during the first run.

### cgMLST schemes

cgMLST schemes need to be provided by the user.
A good place to start for these is [pubMLST](https://pubmlst.org/).
In general it is recommended to look for schemes that are globally make a
consensus for the species of interest.

Importantly, the scheme must be formatted to be used by [CHEWBBACCA](https://github.com/B-UMMI/chewBBACA).

## Configuration

The configuaration can be defined in two ways:

- either edit and locally save the `config/config.yaml` files and provide its path
  to the snakemake command with the `--configfile` argument

- or provide the parameters directly to the snakemake command with
  `--config <ARGNAME>=<VALUE>`

### User defined parameters

Following arguments must be provided for each run:

| Parameter | Type | Description |
| --- | --- | --- |
| `workdir` | path-like string | Path to the ouptut directory |
| `sample_sheet` | path-like string | Path to the sample sheet in TSV format |
| `external_main_clusters` | path-like string | Path to the main cluster naming file |
| `external_sub_clusters` | path-like string | Path to the subcluster naming file |
| `external_profiles` | path-like string | Path to the external allele profiles file |
| `external_timestamps` | path-like string | Path to the external timestamps file |
| `external_statistics` | path-like string | Path to the external allele statistics file |
| `cgmlst_scheme` | path-like string | Path to the cgMLST scheme |

### Optional parameters

Following parameters are optional and will revert to defaults if not set:

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `max_threads_per_job` | integer | 1 | Max number of threads assigned to a single job |
| `chewie_path` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the ChewieSnake pipeline folder |
| `prodigal` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to prodigal training file for the given species |
| `max_missing_loci` | flot [0-1] | 0.05 | Maximumj fraction of loci missing for sample acceptance |
| `distance_method` | int | 3 | Distance calculation method for grapetree |
| `clustering_method` | string | "single" | The agglomeration method to be used for hierarchical <br> clustering. This should be (an unambiguous <br> abbreviation of) one of "ward.D", "ward.D2", "single", <br> "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), <br> "median" (= WPGMC) or "centroid" (= UPGMC) |
| `cluster_distance` | int | 10 | Maximum distance for samples to be attached to a cluster |
| `subcluster_distance` | int | 5 | Maximum distance for samples to be attached to a subcluster |
| `cluster_prefix` | string | "CT" | Prefix for cluster names |

## Usage

The workflow can be started with:

```bash
snakemake --use-conda --conda-prefix <PATH TO CONDA ENVS> --configfile <PATH TO CONFIG> --cores <N>
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more information.

## Input formats

See the input formatting specifications in the `config`folder.

### Sequence files

The assemblies must be provided as fasta file.
Wrapped and unwrapped fastas as well as multifastas are allowed.
There are no special requirements for sequence headers.

## Cluster naming

Cluster naming convention must be according to the following structure:

<PREFIX>-<MAIN#>.<SUB#>

whereby PREFIX is a string identifier for differenciating clusters from e.g. different
species or different sources. This is fixed by the params `cluster_prefix`, make sure you use
the same PREFIX in the external database and while running chewie to avoid colllisions.

MAIN# and SUB# are integer identifying specific clusters and subclusters repectively. They will be
simply incremented for new clusters.

## Results

Results to be used for the next steps are located in the `staging` folder in the workdir.

### Status report

A JSON report of the QC checks in the form

```json
{
    "16-LI00984-0": {
        "STATUS": "PASS",
        "MESSAGES": [
            ""
        ]
    },
    "16-LI01077-0": {
        "STATUS": "PASS",
        "MESSAGES": [
            ""
        ]
    },
    "16-LI00296-0": {
        "STATUS": "PASS",
        "MESSAGES": [
            ""
        ]
    }
}
```

### Isolate datasheets

Isolate information are summarized in single JSON files
Note that these are generated only for samples satisfying all filters.
The follow the same structure, here for a single entry:

```json
{
    "isolate_id": "16-LI00296-0",
    "qc_metrics": {
        "cgmlst_missing_fraction": 0.0035481963335304554
    },
    "cgmlst": {
        "allele_profile": [
            {
                "locus": "lmo0001.fasta",
                "allele_crc32": 3453202319
            },
            {
                "locus": "lmo0002.fasta",
                "allele_crc32": 138852938
            },
            <MANY MORE ALLELES>
        ],
        "allele_stats": {
            "EXC": 1685,
            "INF": 0,
            "LNF": 1,
            "PLOT": 0,
            "NIPH": 2,
            "ALM": 3,
            "ASM": 0
        }
    }
}
```

### Cluster datasheets

Cluster information are summarized in a JSON file per cluster.
Note that JSONs are generated only for clusters with new samples.

```json
{
    "cluster_id": "LIS-1",
    "cluster_number": 1,
    "size": 7,
    "representative": "DB-LI00339-0",
    "AD_threshold": 10,
    "root_members": [
        "16-LI00919-0",
        "16-LI00984-0"
    ],
    "subclusters": [
        {
            "subcluster_id": "LIS-1.2",
            "subcluster_number": 2,
            "size": 2,
            "representative": "DB-LI00986-0",
            "AD_threshold": 5,
            "members": [
                "16-LI00960-0",
                "DB-LI00986-0"
            ]
        },
        {
            "subcluster_id": "LIS-1.1",
            "subcluster_number": 1,
            "size": 3,
            "representative": "DB-LI00339-0",
            "AD_threshold": 5,
            "members": [
                "16-LI01076-0",
                "DB-LI00339-0",
                "DB-LI01201-0"
            ]
        }
    ]
}
```

The orphans samples are listed in a simplified version of the cluster datasheet.
This one is always generated, even if no new oprhans were found.

```json
{
    "cluster_id": "LIS-orphans",
    "cluster_number": 0,
    "size": 15,
    "AD_threshold": 10,
    "members": [
        "16-LI00944-0",
        "16-LI00963-0",
        "16-LI01134-0",
        "DB-LI00345-0",
        "DB-LI00346-0",
        "DB-LI00390-0",
        "DB-LI00572-0",
        "DB-LI00704-0",
        "DB-LI00710-0",
        "DB-LI00732-0",
        "DB-LI00733-0",
        "DB-LI00762-0",
        "DB-LI00842-0",
        "DB-LI01064-0",
        "DB-LI01068-0"
    ]
}
```