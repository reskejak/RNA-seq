#!/bin/bash

# RNA-seq workflow and software commands
# April 2018
# Michigan State University
# reskejak@msu.edu

# Intended for single-end Illumina sequencing

DIR=("/mnt/home/usr/foo")
REFDIR=("/mnt/home/usr/ref_genome") # location of reference genomes
BIN=("/mnt/home/usr/bin")
samples=("foo_treat1 foo_treat2 foo_treatn foo_control1 foo_control2 foo_controln")

#########################

cd ${DIR}

# remove excess filename architecture, paired end reads
for i in $samples;
do mv ${i}_L000_R1_001.fastq.gz ${i}_R1.fastq.gz;
done

#########################
# fastqc on raw files
for i in $samples;
do fastqc ${i}_R1.fastq.gz;
done

##########################
# make trimmed output directory
mkdir ${DIR}/trimmed

# trim galore
for i in $samples;
do trim_galore ${i}_R1.fastq.gz -o trimmed/ ;
done

cd ${DIR}/trimmed

##########################
# fastqc on trimmed files
for i in $samples;
do fastqc ${i}_R1_trimmed.fq.gz;
done

##########################
# alignment via STAR
module load star

# index genome (upon first time using this genome/annotation/read-length with STAR)
cd ${REFDIR}
mkdir genome_STAR_index
STAR --runMode genomeGenerate \
--runThreadN 16 \
--genomeDir ${REFDIR}/genome_STAR_index \
--genomeFastaFiles ${REFDIR}/genome.fa \
--sjdbGTFfile ${REFDIR}/annotation.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 74
# note: --sjdbOverang should equal [read length]-1, or 75-1 = 74 (for SE 75bp reads)
# note: "--sjdbGTFtagExonParentTranscript Parent" must be used with gff3 file
# submit as job; takes <20 minutes but utilizes extensive resources

# align to indexed genome with annotation
# too long for 2 hour wall limit, so submit qsub job
cd ${DIR}/trimmed
for i in $samples;
do gunzip ${i}_R1_trimmed.fq.gz;
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir ${REFDIR}/genome_STAR_index \
--readFilesIn ${i}_R1_trimmed.fq \
--outFileNamePrefix ${i} \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts;
gzip ${i}_R1_trimmed.fq
done
# note: software not currently working with zipped files, call gunzip first...
# note: --outSAMtype can be "SAM", "BAM Unsorted", or "BAM SortedByCoordinate"
# note: --quantMode can either be "GeneCounts" or "TranscriptomeSAM"
# Gene Counts: "count reads per gene"
# TranscriptomeSAM: "output SAM/BAM alignments to transcriptome into a separate file"

############
# Generate read count matrix in R from STAR ReadsPerGene.out.tab files
R

# read STAR output gene counts
foo.treat1.counts <- read.table("foo_treat1_ReadsPerGene.out.tab", sep="\t")
foo.treat2.counts <- read.table("foo_treat2_ReadsPerGene.out.tab", sep="\t")
foo.treatn.counts <- read.table("foo_treatn_ReadsPerGene.out.tab", sep="\t")
foo.control1.counts <- read.table("foo_control1_ReadsPerGene.out.tab", sep="\t")
foo.control2.counts <- read.table("foo_control2_ReadsPerGene.out.tab", sep="\t")
foo.control3.counts <- read.table("foo_controln_ReadsPerGene.out.tab", sep="\t")

# concatenate sample gene count files for Stranded-Reverse (libtype)
counts <- cbind(foo.treat1.counts$V4, foo.treat2.counts$V4, foo.treatn.counts$V4, 
                foo.control1.counts$V4, foo.control2.counts$V4, foo.controln.counts$V4)
# note: select corresponding column for libtype; see STAR manual for more information
        
# add gene identifiers
rownames(counts) <- foo.treat1.counts$V1

# add sample identifiers
colnames(counts) <- c("foo_treat1", "foo_treat2", "foo_treatn", 
                      "foo_control1", "foo_control2", "foo_controln")			

# remove header rows
counts <- counts[-1:-4, ]

# write output table for downstream analysis via DESeq2
write.csv(counts, file="foo_RNA_raw_counts.csv")

