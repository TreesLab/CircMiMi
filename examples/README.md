## Examples

In this example, we use the file `circRNAs.gencode_format.tsv` as the input.


### 1. Generate the references

Before running the main pipeline, we need to generate the reference files first.

```
$ circmimi_tools genref --species hsa --source gencode --version 34 ./refs
```
The above command will automatically download and generate the needed files from the "[gencode release 34](https://www.gencodegenes.org/human/release_34.html)".

<!-- *(explain the options here)* -->


### 2. Run the main pipeline

Next, we run the main pipeline of CircMiMi by the following command:

```
$ circmimi_tools run -r ./refs -i circRNAs.gencode_format.tsv -o . --checkAA --miranfa-sc 170
```

<!-- *(explain the options here)* -->

There are two output files:

- summary_list.tsv
- all_interactions.tsv

For the format of these two files, please check the "[Output files](../README.md#output-files)" section.


### 3. Create the network file

```
$ circmimi_tools network create all_interactions.tsv all_interactions.xgmml
```

<!-- *(explain the options here)* -->
