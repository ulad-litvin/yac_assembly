# yac_assembly
Reference-guided assembly of sequencing results from YAC DNA extractions.

## installation procedure

### 1. install Miniconda

Miniconda is a package manager that helps you to download and install different bioinformatic tools.
Instructions on how to install Miniconda on your computer can be found here https://docs.anaconda.com/free/miniconda/ .

### 2. clone this repository

Download this repository to your computer.

```bash
git clone https://github.com/ulad-litvin/yac_assembly.git
```

It contains:
- `ref_assembly.py`: python script that does the reference-guided assembly;
- `environment.yaml`: list of required tools and dependencies;
- `example_dataset`: test dataset.

This folder should be the place where you run your analysis.

### 3. create a conda environment

To perform the reference-guided assembly you need to install several bioinformatic tools and their dependencies.
This command creates a new conda environment called `yac_assembly`, downloads and installs all required tools (make sure that Miniconda is installed).

```bash
cd ./yac_assembly
```
```bash
conda activate
```
```bash
conda env create -f environment.yml
```

### 4. install [weeSAM](https://github.com/centre-for-virus-research/weeSAM)

weeSAM is a bioinformatic tool that provides you some basic statistics (number of reads mapped to the reference sequence, coverage, depth, etc.).
It creates `_html_results` folder with html and svg files. 

```bash
git clone https://github.com/centre-for-virus-research/weeSAM.git 
```


## running procedure

### 1. activate the `yac_assembly` environment

`yac_assembly` environment needs to be activated before you can run the `ref_assembly.py` script.
Otherwise your computer assumes that tools used in the script are not installed.

```bash
conda activate yac_assembly
```

### 2. reference sequence and a folder with raw reads

Inside the `yac_assembly` folder:
- have a file with the reference seqence. It should be in the FASTA format (eg. `yac-a_ref.fasta`).
- have a folder with one or more raw reads files (eg. `raw_reads`). It should contain all FASTQ files that you want to align to the reference sequence (eg. `yac-a_sample1.fast1`, `yac_a_sample2.fastq`, etc.).

Examples of the `yac-a_ref.fasta` and `raw_reads` folders with FASTA and FASTQ files can be found in the `example_dataset` folder.

### 3. run `ref_assembly.py` script

`ref_assembly.py` script can produce consensus sequences for multiple raw reads files by aligning them to the same reference sequence.

You need to provide two arguments:
- a file with your reference sequence (eg. `yac-a_ref.fasta`)
- a folder with raw reads (eg. `raw_reads`).

Optionally you can provide a folder name where consesus sequences will be stored (by default `ref_assembly.py` creates a folder called `results`).

```bash
python ref_assembly.py ./yac-a_ref.fasta ./raw_reads ./results
```


## results

`ref_assembly.py` script creates:
- a folder called `results` with consensus sequences for `raw_reads`.
- a folder called `_html_results` with html and svg files with some basic statistics on consensus sequences.
