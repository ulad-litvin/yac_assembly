"""
--------------------------------
@Title: Reference assembly of sequencing results from DNA extraction
@Description: This script is used to assemble raw reads from DNA extraction using any reference sequence
@Author: ulad-litvin with GitHub Copilot
@Date: 16/05/2024
--------------------------------
"""


"""
--------------------------------
Importing the required libraries
--------------------------------
"""

import sys
import os
import shutil
import subprocess


"""
--------------------------------
Assigning the reference sequence
--------------------------------
"""

try:
    ref_seq_folder = sys.argv[1]
    ref_seq_file = [f for f in os.listdir(ref_seq_folder) if f.endswith('.fasta')]
    
    if len(ref_seq_file) == 0:
        print("No FASTA files found in the reference sequence folder.")
        sys.exit(1)
    
    if len(ref_seq_file) > 1:
        print("Multiple FASTA files found in the reference sequence folder.\nPlease provide only one FASTA file per reference sequence.")
        sys.exit(1)
    
    ref_seq_path = os.path.join(ref_seq_folder, ref_seq_file[0])
    print(ref_seq_path)

except IndexError:
    print("Please provide the folder with your reference sequence.")
    sys.exit(1)


"""
--------------------------------
Assigning the files with raw reads
--------------------------------
"""

try:
    samples_folder = sys.argv[2]
    samples_files = [f for f in os.listdir(samples_folder) if f.endswith('.fastq')]
    print(samples_files)
    
    if len(samples_files) == 0:
        print("No FASTQ files found in the raw reads folder.")
        sys.exit(1)
    
    samples_path = [os.path.join(samples_folder, sample) for sample in samples_files]
    print(samples_path)

    samples_path_no_ext = [sample.rstrip('.fastq') for sample in samples_path]
    print(samples_path_no_ext)

    sample_names = [sample.split('/')[-1] for sample in samples_path_no_ext]
    print(sample_names)

    # make a dictionary with sample names as keys and sample paths as values
    samples_dict = dict(zip(sample_names, samples_path_no_ext))
    print(samples_dict)

except IndexError:
    print("Please provide the path to a folder with your raw reads.")
    sys.exit(1)


"""
--------------------------------
Creating the output folder if it does not exist
--------------------------------
"""

try:
    output_folder = sys.argv[3]
    os.makedirs(output_folder, exist_ok=True)
except IndexError:
    output_folder = './results'
    os.makedirs(output_folder, exist_ok=True)


"""
--------------------------------
Indexing the reference sequence
--------------------------------
"""

cmd = f"bwa index {ref_seq_path}"
subprocess.Popen([cmd], shell = True, close_fds=True).wait()


"""
--------------------------------
Produce the consensus sequences and statistics for each sample
--------------------------------
"""

for sample_name, sample_path_no_ext in samples_dict.items():

    # align trimmed reads to reference with bwa mem
    cmd = f"bwa mem {ref_seq_path} {sample_path_no_ext}.fastq > {sample_path_no_ext}.sam"
    subprocess.Popen([cmd], shell = True,close_fds=True).wait()

    # convert sam to sorted bam file
    cmd = f"samtools view -S -b {sample_path_no_ext}.sam | samtools sort > {sample_path_no_ext}.bam"
    subprocess.Popen([cmd], shell = True,close_fds=True).wait()

    # index bam file
    cmd = f"samtools index {sample_path_no_ext}.bam"
    subprocess.Popen([cmd], shell = True,close_fds=True).wait()

    # make concensus
    cmd = f"samtools consensus --ambig {sample_path_no_ext}.bam -o {sample_path_no_ext}_con.fa"
    subprocess.Popen([cmd], shell = True,close_fds=True).wait()

    # get summary statistics
    cmd = f"samtools stats {sample_path_no_ext}.bam > {sample_path_no_ext}.bam.bc"
    subprocess.Popen([cmd], shell = True,close_fds=True).wait()

    # make coverage plots using weeSAM
    cmd = f"./bin/weeSAM/weeSAM --bam {sample_path_no_ext}.bam --html {sample_path_no_ext}"
    subprocess.Popen([cmd], shell = True,close_fds=True).wait()

    # move consensus files to the output folder
    shutil.move(f"{sample_path_no_ext}_con.fa", output_folder)

