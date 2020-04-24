# CircMiMi

A toolset for investigating the interactions between circRNA - miRNA - mRNA.


# Table of Contents
- [Dependency](#dependency)
    - [Packages](#packages)
    - [External tools](#external-tools)
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
- numpy
- pandas
- xlrd
- networkx
- lxml

### External tools

- bedtools (https://github.com/arq5x/bedtools2)
- miranda (http://www.microrna.org/microrna/getDownloads.do)
- blat (https://genome.ucsc.edu/FAQ/FAQblat.html)


# Installation

The recommended way is via `conda`, a package and environment management system. (https://docs.conda.io/en/latest/)


You may install `circmimi` by the following steps:
```
$ conda create -n circmimi python3
$ conda activate circmimi
$ pip install circmimi
```

For the external tools, they can also be installed via `conda` with the `bioconda`(https://bioconda.github.io/) channel:
```
$ conda install -c bioconda bedtools miranda blat
```



Now, you can try the following command to test the installation,
```
$ circmimi_tools --help
```
it should print out with the help messages.



# Quick Start

1. Generate the references
```
$ circmimi_tools genref --species hsa --source ensembl --version 98 ./refs
```


2. Run the main pipeline of CircMiMi

```
$ circmimi_tools run -r ./refs -c circRNAs.tsv -o ./out/ -p 5
```


3. Create the network file for Cytoscape
```
$ circmimi_tools network create ./out/out.tsv ./out/out.xgmml
```


# Usage
## Generate the references

```
circmimi_tools genref --species SPECIES --source SOURCE [--version RELEASE_VER] REF_DIR
```

### Options
Option                | Description
:-------------------- | :------------------------------
--species SPECIES     | Assign the species for references. Use the species code for SPECIES. ***[required]***
--source SOURCE       | Available values for SOURCE: "ensembl", "ensembl_plants", "ensembl_metazoa", "gencode". ***[required]***
--version RELEASE_VER | The release version of the SOURCE. For examples,  "98" for ("hsa", "ensembl"), "M24" for ("mouse", "gencode"). If the version is not assigned, the latest will be used.




### Available species and sources

Code | Name                    |  E  |  G  |  EP |  EM | MTB | MDB |
:--  | :---------------------- | :-: | :-: | :-: | :-: | :-: | :-: |
ath  | Arabidopsis thaliana    |     |     |  V  |     |  V  |     |
bmo  | Bombyx mori             |     |     |     |  V  |  V  |     |
bta  | Bos taurus              |  V  |     |     |     |  V  |     |
cel  | Caenorhabditis elegans  |  V  |     |     |  V  |  V  |     |
cfa  | Canis familiaris        |  V  |     |     |     |  V  |  V  |
dre  | Danio rerio             |  V  |     |     |     |  V  |     |
dme  | Drosophila melanogaster |  V  |     |     |     |  V  |     |
gga  | Gallus gallus           |  V  |     |     |     |  V  |  V  |
hsa  | Homo sapiens            |  V  |  V  |     |     |  V  |  V  |
mmu  | Mus musculus            |  V  |  V  |     |     |  V  |  V  |
osa  | Oryza sativa            |     |     |  V  |     |  V  |     |
ola  | Oryzias latipes         |  V  |     |     |     |  V  |     |
oar  | Ovis aries              |  V  |     |     |     |  V  |     |
rno  | Rattus norvegicus       |  V  |     |     |     |  V  |  V  |
ssc  | Sus scrofa              |  V  |     |     |     |  V  |     |
tgu  | Taeniopygia guttata     |  V  |     |     |     |  V  |     |
xtr  | Xenopus tropicalis      |  V  |     |     |     |  V  |     |

**E**: Ensembl, **G**: Gencode, **EP**: Ensembl Plants, **EM**: Ensembl Metazoa, **MTB**: miRTarBase, **MDB**: miRDB



## Run the main pipeline

```
circmimi_tools run --ref REF_DIR --circ CIRC_FILE [-o OUT_PREFIX] [-p NUM_PROC] [--checkAA]
```

### Options
Option                      | Description
:-------------------------- | :------------------------------
-c, --circ CIRC_FILE        | The file of circRNAs. ***[required]***
-r, --ref REF_DIR           | The directory of the pre-genereated reference files. ***[required]***
-o, --out-prefix OUT_PREFIX | Assign the prefix for the output filenames. (default: "./out/")
-p, --num_proc NUM_PROC     | Assign the number of processes.
--checkAA                   | Check the circRNAs if there are ambiguous alignments.


### Input format

The input format for the CIRC_FILE.

\#   | Column  | Description
:--: | :-----: | :----------
  1  |  chr    | Chromosome name
  2  |  pos1   | One of the position of the circRNA junction site
  3  |  pos2   | Another position of the circRNA junction site
  4  |  strand | + / -

#### Note.
- The chromosome name must be the same as the name in the SOURCE.
  - For example, "1" for "ensembl", and "chr1" for "gencode".


### Output format
CircMiMi appends the following columns to the original input.

\#   | Column          | Description
:--: | :-------------: | :----------
  5  |  host_gene      | Host gene of the circRNA
  6  |  mirna          | The miRNA which may bind on the circRNA
  7  |  max_score      | The maximum binding score reported by miRanda
  8  |  count          | The number of the miRNA-binding sites on the circRNA
  9  |  cross_boundary | If there is a binding site across the junction of the circRNA
 10  |  target_gene    | The miRNA-targeted gene


And additional columns from the miRNA-target interactions database.

 \#   | Column          | Description
:--: | :-------------: | :----------
 11  |  miRTarBase      | If the miRNA-mRNA interaction is from miRTarBase
 12  |  miRDB      | If the miRNA-mRNA interaction is from miRDB
 13  |  miRTarBase__ref_count | The number of references which support the interaction
 14  |  miRDB__targeting_score | The predicted target score from miRDB


## Create the network file for Cytoscape

```
circmimi_tools network create IN_FILE OUT_FILE
```

We may use the subcommand to create the network file (in XGMML format), which can be load into the Cytoscape to visualize the network.
