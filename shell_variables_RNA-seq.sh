#!/bin/bash

# CPU threads
CPU=7
# output directory for raw read QC
QC_dir="./QC_RNA-seq/"
# output directory for quality-trimmed read QC
QC_trim_dir="./QC_trim_RNA-seq/"
# FASTQ file directory
fastq_dir="./FASTQ_RNA-seq/"
# quality-trimmed FASTQ file directory
fastq_trim_dir="./FASTQ_trim_RNA-seq/"
# file containing the list of all FASTQ files to be processed
fastq_file_list="FASTQ_files_RNA-seq"
# file containing the list of all quality-trimmed FASTQ files to be processed
fastq_trim_file_list="FASTQ_trim_files_RNA-seq"
# FASTQ file extension
fastq_file_ext="\.fastq\.gz$"
# genome sequence and annotation folder
genome_dir="./genome/"
# file containing reference genome sequence
genome="${genome_dir}Schizosaccharomyces_pombe.fa"
# BAM file directory
bam_dir="./BAM_RNA-seq/"
# how to perform binning of genome coverage
bin_size=1
# images directory
image_dir="./images/"
# file containing the list of all BAM files to be processed
bam_file_list="BAM_files_RNA-seq"
# output of multiBamSummary
bam_summary_file="multiBamSummary.npz"
# directory for genome coverage data
coverage_dir="./coverage_RNA-seq/"

