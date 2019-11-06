# CircMiMi

A toolset for investigating the interactions between circRNA - miRNA - mRNA.


# Table of Contents
- [Dependency](#dependency)
    - [Packages](#packages)
    - [Externel tools](#externel-tools)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
    - [Generate the index and references](#generate_the_index_and_references)
        - [Available species](#available_species)
    - [Run the main pipeline](#run_the_main_pipeline)
      - [Input format](#input_format)
      - [Output format](#output_format)


# Dependency

### Packages

- python3
- click
- sqlalchemy
- pandas
- xlrd

### Externel tools

- bedtools (https://github.com/arq5x/bedtools2)
- miranda (http://www.microrna.org/microrna/getDownloads.do)


# Installation

The recommended way is through `conda`, a package and environment management system. (https://docs.conda.io/en/latest/)


You may install `circmimi` by the following steps:

1. Create a new environment for `circmimi` (optional, but recommended!),
```
$ conda create -n circmimi python3
```
and activate this new environment by:
```
$ conda activate circmimi
```

2. We may install the externel tools with `conda` through the `bioconda`(https://bioconda.github.io/) channel:
```
$ conda install -c bioconda bedtools miranda
```

3. Finally, install the `circmimi` package with `pip`:
```
$ pip install circmimi
```


Now, you can try the following command,
```
$ circmimi_tools --help
```
it should print out with the help messages.



# Quick Start

1. Generate the index and references
```
$ circmimi_tools genref --species hsa --ensembl --version 98 ./refs
```


2. Run the main pipeline of CircMiMi

```
$ circmimi_tools run -r ./refs circ_events.tsv > out.tsv
```



# Usage
## Generate the index and references
```
Usage: circmimi_tools genref [OPTIONS] REF_DIR

Options:
  --species TEXT
  --source TEXT
  --gencode
  --ensembl
  --version TEXT
  --init          Create an init template ref_dir.
  --help          Show this message and exit.
```

##### Example.
```
$ circmimi_tools genref --species hsa --source ensembl --version 98 ./refs
```
or just
```
$ circmimi_tools genref --species hsa --ensembl --version 98 ./refs
```

### Available species

Key | Name                    | Ensembl | Gencode | Alternative Source
:-- | :---------------------- | :-----: | :-----: | :------------------
ath | Arabidopsis thaliana    |    △    |         | Ensembl Plants
bmo | Bombyx mori             |    △    |         | Ensembl Metazoa
bta | Bos taurus              |    ✔️   |         |
cel | Caenorhabditis elegans  |    ✔️   |         | Ensembl Metazoa
cfa | Canis familiaris        |    ✔️   |         |
dre | Danio rerio             |    ✔️   |         |
dme | Drosophila melanogaster |    ✔️   |         |
gga | Gallus gallus           |    ✔️   |         |
hsa | Homo sapiens            |    ✔️   |   ✔️    |
mmu | Mus musculus            |    ✔️   |   ✔️    |
osa | Oryza sativa            |    △    |         | Ensembl Plants
ola | Oryzias latipes         |    ✔️   |         |
oar | Ovis aries              |    ✔️   |         |
rno | Rattus norvegicus       |    ✔️   |         |
ssc | Sus scrofa              |    ✔️   |         |
tgu | Taeniopygia guttata     |    ✔️   |         |
xtr | Xenopus tropicalis      |    ✔️   |         |

△ : Only in the alternative source

##### Note:
To access the alternative source, just assign the source name to the "source" option:
```
$ circmimi_tools genref --species ath --source ensembl_plants --version 45 ./refs
```


## Run the main pipeline
```
Usage: circmimi_tools run [OPTIONS] CIRC_FILE

  Main pipeline.

Options:
  -r, --ref PATH          Assign the path of ref_dir.  [required]
  -p, --num_proc INTEGER  The number of processes.
  --no-header
  --help                  Show this message and exit.
```

##### Example.
```
$ circmimi_tools run -r ./refs -p 10 circ_events.tsv > output.tsv
```

### Input format

eg. circ_events.tsv

\#   | Column  | Description
:--: | :-----: | :----------
  1  |  chr    | Chromosome name
  2  |  pos1   | The 1st position of circRNA junction
  3  |  pos2   | The 2nd position of circRNA junction
  4  |  strand | + / -

### Output format
CircMiMi appends the following columns to the original input.

\#   | Column          | Description
:--: | :-------------: | :----------
  5  |  host_gene      | Host gene of the circRNA
  6  |  mirna          | The miRNA which may bind on the circRNA
  7  |  max_score      | The max binding score
  8  |  count          | The number of binding sites of the miRNA on the circRNA
  9  |  cross_boundary | Is there a binding site cross the junction of circRNA
 10  |  target_gene    | The miRNA-targeted gene
 11  |  ref_count      | The references count from the miRTarBase

