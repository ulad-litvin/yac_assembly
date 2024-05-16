# yac_assembly
Reference-guided assembly of sequencing results from YAC DNA extractions.

## installation procedure

### 1. install Miniconda

Miniconda is a package manager that helps you to download and install different bioinformatic tools.
Instructions on how to install Miniconda on your computer can be found here https://docs.anaconda.com/free/miniconda/ .

### 2. clone this repository

You need to download this repository to your computer. It contains:
- `ref_assembly.py` script: python script that does the reference assembly;
- [`weeSAM` tool](https://github.com/centre-for-virus-research/weeSAM): tool that provides you with coverage statistics;
- `environment.yaml`: list of tools that will be installed by conda when you create a new environment;
- `example_dataset`: test dataset.

This folder will become a place where you run your analysis.

```bash
git clone https://github.com/ulad-litvin/yac_assembly.git
```

### 3. create a new conda environment

To perform the reference-guided assembly you need to install several bioinformatic tools and their dependencies.
This command creates a new conda environment called `yac_assembly`, downloads and installes all required tools.

```bash
cd ./yac_assembly
conda activate
conda env create -f environment.yml
```

## running procedure

### 1. activate the `yac_assembly` environment

`yac_assembly` environment need to be activated before you can run the `ref_assembly.py` script.
Otherwise your computer assumes that tools used in the script are not installed.

```bash
conda activate yac_assembly
```

### 2. make folders with reference sequence and raw reads

Inside the `yac_assembly` folder make two additional folder:
- a folder for the reference seqence (eg. `reference`). It should contain a single file of your reference sequence in the FASTA format (eg. `yac_a_ref.fasta`).
- a folder with one or more raw reads files (eg. `raw_reads`). It should contain all FASTQ files that you want to align to the reference sequence (eg. `yac_a_sample1.fast1`, `yac_a_sample2.fastq`, etc.).

### 3. run `ref_assembly.py` script

`ref_assembly.py` script can produce consensus sequences for multiple raw reads files by aligning them to the same reference sequence.
You need to provide two arguments: a folder with your reference sequence (should be .fasta) and a folder with raw reads (should be .fastq).
Optionally you can provide a folder name where consesus sequences will be stored (by default `ref_assembly.py` creates a folder called results).

```bash
python ref_assembly.py ./reference ./raw_reads ./results
```
