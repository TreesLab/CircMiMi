## Examples

In this example, we use the file `circRNAs.gencode_format.tsv` as the input.


### 1. Generating the references

To use the CircMiMi package, we need to generate the reference files first.

```
$ circmimi_tools genref --species hsa --source gencode --version 34 refs/
```
The above command will automatically download and generate the needed files from the "[gencode release 34](https://www.gencodegenes.org/human/release_34.html)" and other databases.

<!-- *(explain the options here)* -->


### 2. Checking the circRNAs

```
$ circmimi_tools checking -r refs/ -i circRNAs.gencode_format.tsv -o circRNAs.gencode_format. --dist 5000
```

```
$ cat circRNAs.gencode_format.checking.results.tsv | awk -F'\t' '($9==1)&&($12==0)&&($16==1)' | cut -f '-5' > circRNAs.gencode_format.filtered.tsv
```


### 3. Predicting the interactions

```
$ circmimi_tools interactions -r refs/ -i circRNAs.gencode_format.filtered.tsv -o circRNAs.gencode_format. --miranda-sc 150
```

<!-- *(explain the options here)* -->

There are two output files:

- summary_list.tsv
- all_interactions.miRNA.tsv

For the format of these two files, please check the "[Output files](../README.md#output-files)" section.


### 4. Visualization for the interactions

```
$ circmimi_tools visualize circRNAs.gencode_format.all_interactions.miRNA.tsv circRNAs.gencode_format.all_interactions.miRNA.xgmml
```

<!-- *(explain the options here)* -->
