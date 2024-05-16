# yac_assembly
Reference assembly of sequencing results from YAC DNA extraction

## installation procedure

### 1. install Miniconda on your Mac

Miniconda is a package manager that helps you to download and install different bioinformatic tools.
Instructions on how to install Miniconda can be found here https://docs.anaconda.com/free/miniconda/ .

### 2. clone this repository to your Mac

You need to download this repository to your computer. This folder will be a place where you run your analysis.

```bash
git clone https://github.com/ulad-litvin/yac_assembly.git
```

### 3. create a new conda environment

This step will download all conda packages essential for reference-guided assembly.

```bash
cd ./yac_assembly
conda activate
conda env create -f environment.yml
```

## running procedure

### 1. activate the `yac_assembly` environment

```bash
conda activate yac_assembly
```

### 2. make folders with reference sequence and raw reads

`ref_assembly.py` script can produce consensus sequences for multiple raw reads files at a time by aligning them to the same reference sequence.

### 3. 

provide folders with your reference sequence raw_reads (optionally you can provide a name for the results folder)

```bash
python ref_assembly.py ./reference ./raw_reads ./results
```
