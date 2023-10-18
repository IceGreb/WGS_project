#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp

if false 
then 
#download data
wget -P /home/nikos/WGS_project/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P /home/nikos/WGS_project/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


#download reference 
wget -P /home/nikos/WGS_project/ref_files/hg38  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /home/nikos/WGS_project/ref_files/hg38/hg38.fa.gz


#index ref - .fai file  before running haplotype caller
samtools faidx /home/nikos/WGS_project/ref_files/hg38/hg38.fa


#ref dict  - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=/home/nikos/WGS_project/ref_files/hg38/hg38.fa O=/home/nikos/WGS_project/ref_files/hg38/hg38.dict


#download known sites files for BQSR from GATK resource bundle
wget -P /home/nikos/WGS_project/ref_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /home/nikos/WGS_project/ref_files/hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

#############   VARIANT CALLING STEPS   ######################
fi

# directories
ref="/home/nikos/WGS_project/ref_files/hg38/hg38.fa"
known_sites="/home/nikos/WGS_project/ref_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/nikos/WGS_project/aligned_reads"
reads="/home/nikos/WGS_project/reads"
results="/home/nikos/WGS_project/results"
data="/home/nikos/WGS_project/data"


#   =======================
#   STEP1 : QC - Run fastqc
#   =======================

echo "STEP 1: QC    -   Run fastqc"

#fastqc  ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
#fastqc  ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# Trimming required??


#   ======================================
#   STEP 2: Map to reference using BWA-MEM  (BWA also indexes reference)
#   ======================================

echo "STEP 2: Map to reference using BWA-MEM"

#BWA index reference
#bwa index ${ref}


#BWA alignment
#bwa mem     -t  4   -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam


#       ===================================
#       STEP 3 : Mark Duplicates and Sort = GATK4
#       ===================================

echo "STEP 3: Mark Duplicates and Sort - GATK4"

#gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam   -O  ${aligned_reads}/SRR062634_sorted_dedup_reads.bam    --spark-master local[10] --tmp-dir ~/tmp  #adjusted tmp-dir for solving memory issue!


#       ============================
#       STEP 4: Base Quality Recalibration
#       ============================

echo "STEP 4: Base  Quality Recalibration"

#1. Build the model
#gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref}  --known-sites ${known_sites}  -O ${data}/recal_data.table

#2. Apply the model to adjust the base quality scores
#gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam



#           ======================================
#           STEP 5: Collect alignment & Insert  Size Metrics
#           ======================================


echo "STEP 5: Collect alignment & Insert Size Metrics"

#gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
#gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt  HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


#           ===================================
#           STEP 6: Call Variants - gatk haplotype caller
#           ===================================


echo "STEP 6: Call Variants - gatk haplotype caller"

#gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf


#  extract SNPs & INDELS

gatk SelectVariants  -R  ${ref}  -V ${results}/raw_variants.vcf  --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants  -R ${ref}   -V ${results}/raw_variants.vcf  --select-type INDEL -O ${results}/raw_indels.vcf
