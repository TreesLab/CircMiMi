# CircMiMi

A toolset for investigating the interactions between circRNA - miRNA - mRNA.


# Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
    - [Generate the references](#generate-the-references)
      - [Parameters](#parameters)
      - [Available species and sources](#available-species-and-sources)
    - [Run the main pipeline](#run-the-main-pipeline)
      - [Parameters](#parameters-1)
      - [Input file](#input-file)
      - [Output files](#output-files)
        - [summary_list.tsv](#summary_listtsv)
        - [all_interactions.tsv](#all_interactionstsv)
    - [Create the network file for Cytoscape](#create-the-network-file-for-cytoscape)
      - [Parameters](#parameters-2)


# Requirements

- Python (3.6 or above)
- External tools
  - bedtools (2.29.0) (https://github.com/arq5x/bedtools2)
  - miranda (aug2010, 3.3a) (http://www.microrna.org/microrna/getDownloads.do)
  - blat (https://genome.ucsc.edu/FAQ/FAQblat.html)


# Installation

The recommended way is via `conda`, a package and environment management system. (https://docs.conda.io/en/latest/)


You may install `circmimi` by the following steps:
```bash
$ conda create -n circmimi python3
$ conda activate circmimi
$ pip install circmimi
```

For the external tools, they can also be installed via `conda` with the `bioconda`(https://bioconda.github.io/) channel:
```bash
$ conda install -c bioconda bedtools miranda blat
```



Now, you can try the following command to test the installation,
```bash
$ circmimi_tools --help
```
it should print out with the help messages.



# Quick Start

1. Generate the references

```bash
$ circmimi_tools genref --species hsa --source ensembl --version 98 ./refs
```


2. Run the main pipeline of CircMiMi

```bash
$ circmimi_tools run -r ./refs -i circRNAs.tsv -o ./out/ -p 5 --checkAA --miranda-sc 175
```


3. Create the network file for Cytoscape

```bash
$ circmimi_tools network create ./out/out.tsv ./out/out.xgmml
```


# Usage
## Generate the references

```
circmimi_tools genref --species SPECIES --source SOURCE [--version RELEASE_VER] REF_DIR
```

### Parameters
Parameter             | Description
:-------------------- | :------------------------------
--species SPECIES     | Assign the species for references. Use the species code for SPECIES. ***[required]***
--source SOURCE       | Available values for SOURCE: "ensembl", "ensembl_plants", "ensembl_metazoa", "gencode". ***[required]***
--version RELEASE_VER | The release version of the SOURCE. For examples,  "98" for ("hsa", "ensembl"), "M24" for ("mmu", "gencode"). If the version is not assigned, the latest one will be used.




### Available species and sources

Code | Name                    |  E  |  G  |  EP |  EM |  MB | MTB | MDB |
:--  | :---------------------- | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
ath  | Arabidopsis thaliana    |     |     |  V  |     |  V  |  V  |     |
bmo  | Bombyx mori             |     |     |     |  V  |  V  |  V  |     |
bta  | Bos taurus              |  V  |     |     |     |  V  |  V  |     |
cel  | Caenorhabditis elegans  |  V  |     |     |  V  |  V  |  V  |     |
cfa  | Canis familiaris        |  V  |     |     |     |  V  |  V  |  V  |
dre  | Danio rerio             |  V  |     |     |     |  V  |  V  |     |
dme  | Drosophila melanogaster |  V  |     |     |     |  V  |  V  |     |
gga  | Gallus gallus           |  V  |     |     |     |  V  |  V  |  V  |
hsa  | Homo sapiens            |  V  |  V  |     |     |  V  |  V  |  V  |
mmu  | Mus musculus            |  V  |  V  |     |     |  V  |  V  |  V  |
osa  | Oryza sativa            |     |     |  V  |     |  V  |  V  |     |
ola  | Oryzias latipes         |  V  |     |     |     |  V  |  V  |     |
oar  | Ovis aries              |  V  |     |     |     |  V  |  V  |     |
rno  | Rattus norvegicus       |  V  |     |     |     |  V  |  V  |  V  |
ssc  | Sus scrofa              |  V  |     |     |     |  V  |  V  |     |
tgu  | Taeniopygia guttata     |  V  |     |     |     |  V  |  V  |     |
xtr  | Xenopus tropicalis      |  V  |     |     |     |  V  |  V  |     |

###### Gene annotation
   - **E**: Ensembl (https://www.ensembl.org/index.html)
   - **G**: Gencode (https://www.gencodegenes.org/)
   - **EP**: Ensembl Plants (https://plants.ensembl.org/index.html)
   - **EM**: Ensembl Metazoa (https://metazoa.ensembl.org/index.html)

###### Database for miRNAs
   - **MB**: miRBase (v22) (http://www.mirbase.org/)

###### Databases for miRNA-mRNA interactions
   - **MTB**: miRTarBase (v7.0) (http://mirtarbase.mbc.nctu.edu.tw/php/index.php)
   - **MDB**: miRDB (v6.0) (http://mirdb.org/)



## Run the main pipeline

```
circmimi_tools run -r REF_DIR -i CIRC_FILE [-o OUT_PREFIX] [-p NUM_PROC] [--checkAA] \
[--miranda-sc SCORE] [--miranda-en ENERGY] [--miranda-scale SCALE] [--miranda-strict] [--miranda-go X] [--miranda-ge Y]
```

### Parameters
Parameter                   | Description
:-------------------------- | :------------------------------
-r, --ref REF_DIR           | The directory of the pre-genereated reference files. ***[required]***
-i, --circ CIRC_FILE        | The file of circRNAs. ***[required]***
-o, --out-prefix OUT_PREFIX | Assign the prefix for the output filenames. (default: "./out/")
-p, --num_proc NUM_PROC     | Assign the number of processes.
--checkAA                   | Check the circRNAs if there are ambiguous alignments.

The miRanda parameters are also available (see [the manual of miRanda](http://cbio.mskcc.org/microrna_data/manual.html)).

Parameters | Description
:-------------------------- | :------------------------------
--miranda-sc SCORE | Set the alignment score threshold to SCORE. Only alignments with scores >= SCORE will be used for further analysis. (default: 140.0)
--miranda-en ENERGY | Set the energy threshold to ENERGY. Only alignments with energies <= ENERGY will be used for further analysis. A negative value is required for filtering to occur. (default: 1.0)
--miranda-scale SCALE | Set the scaling parameter to SCALE. This scaling is applied to match / mismatch scores in the critical 7bp region near the 5' end of the microRNA. Many known examples of miRNA:Target duplexes are highly complementary in this region. This parameter can be thought of as a contrast function to more effectively detect alignments of this type. (default: 4.0)
--miranda-strict | Require strict alignment in the seed region (offset positions 2-8). This option prevents the detection of target sites which contain gaps or non-cannonical base pairing in this region.
--miranda-go X | Set the gap-opening penalty to X for alignments. This value must be negative. (default: -4.0)
--miranda-ge Y | Set the gap-extend penalty to Y for alignments. This value must be negative. (default: -9.0)



### Input file

The input file(CIRC_FILE) is a TAB-separated file with the following columns:

\#   | Column  | Description
:--: | :-----: | :----------
  1  |  chr    | Chromosome name
  2  |  pos1   | One of the position of the circRNA junction site
  3  |  pos2   | Another position of the circRNA junction site
  4  |  strand | + / -

#### Note.
- The chromosome name must be the same as the name in the SOURCE.
  - For example, "1" for "ensembl", and "chr1" for "gencode".


### Output files
The main pipeline of CircMiMi outputs two main files: "summary_list.tsv" and "all_interactions.tsv".


#### summary_list.tsv
The summary list contains the counts of interactions and some checking results of the circRNAs.

\#   | Column          | Description
:--: | :-------------- | :----------
  1  |  chr            | Chromosome name
  2  |  pos1           | One of the position of the circRNA junction site
  3  |  pos2           | Another position of the circRNA junction site
  4  |  strand         | + / -
  5  |  pass           | 'yes' for the circRNA passing all of the checking items (column 9 to 13). Otherwise 'no'.
  6  |  #circRNA_miRNA | Count for the circRNA-miRNA interactions.
  7  |  #circRNA_mRNA  | Count for the circRNA-mRNA interactions.
  8  |  #circRNA_miRNA_mRNA | Count for the circRNA-miRNA-mRNA interactions.
  9  |  donor site not at the annotated boundary | '1' if the donor site of the circRNA is NOT at the annotated exon boundary. Otherwise '0'.
 10  |  acceptor site not at the annotated boundary | '1' if the acceptor site of the circRNA is NOT at the annotated exon boundary. Otherwise '0'.
 11  |  donor/acceptor sites not at the same transcript isoform | '1' if the donor and acceptor are not at the same annotated transcript isoform. Otherwise '0'.
 12  |  ambiguity with an co-linear explanation | '1' if the merged flanking sequence of the circRNA junction sites has an co-linear explanation. Otherwise '0'.
 13  |  ambiguity with multiple hits | '1' if the merged flanking sequence of the circRNA junction sites is with multiple hits. Otherwise '0'.


#### all_interactions.tsv

\#   | Column          | Description
:--: | :-------------- | :----------
  1  |  chr            | Chromosome name
  2  |  pos1           | One of the position of the circRNA junction site
  3  |  pos2           | Another position of the circRNA junction site
  4  |  strand         | + / -
  5  |  host_gene      | Host gene of the circRNA
  6  |  mirna          | The miRNA which may bind on the circRNA
  7  |  max_score      | The maximum binding score reported by miRanda
  8  |  count          | The number of the miRNA-binding sites on the circRNA
  9  |  cross_boundary | If there is a binding site across the junction of the circRNA
 10  |  target_gene    | The miRNA-targeted gene
 11  |  miRTarBase      | Is the miRNA-mRNA interaction reported from miRTarBase
 12  |  miRDB      | Is the miRNA-mRNA interaction reported from miRDB
 13  |  miRTarBase__ref_count | The number of references reporting the interaction
 14  |  miRDB__targeting_score | The predicted target score from miRDB


## Create the network file for Cytoscape

```
circmimi_tools network create IN_FILE OUT_FILE
```

### Parameters
Parameter     | Description
:------------ | :------------------------------
IN_FILE       | Input the file "all_interactions.tsv" produced from the CircMiMi main pipeline.
OUT_FILE      | The output filename. The file extension should be ".xgmml" or ".xml", so that the Cytoscape would recognize this file as an XGMML network file.


This command can generate a Cytoscape-executable file (.xgmml) for visualization of the input circRNA-miRNA-mRNA regulatory axes in Cytoscape.