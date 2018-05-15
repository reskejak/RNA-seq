#!/bin/bash

###########################
# Tuxedo Pipeline
# full automated script for Illumina PE 100 bp sequencing quality control and Tuxedo analysis
#
# Jake Reske
# Michigan State University
# April 2017
# reskejak@msu.edu
###########################

# computation completes in under 72 hours with 8 cores and 32 gb RAM
# trimmomatic --> fastqc --> tophat --> cufflinks --> cuffmerge --> cuffdiff
#							      \--> cuffnorm

# REQUIRED GENOME ASSEMBLY FILES
# ensure links or original files for genome assembly are in working directory
# genes.gtf
# genome.1.bt2
# genome.2.bt2
# genome.3.bt2
# genome.4.bt2
# genome.fa
# genome.fai
# genome.rev.1.bt2
# genome.rev.2.bt2

# $HOME/raw_reads - example directory for raw reads
# FORMAT: [descriptor]_[sample]_[condition]_R[1/2].fastq.gz, i.e.
# foo_1_treat_R1.fastq.gz
# foo_1_treat_R2.fastq.gz
# foo_2_treat_R1.fastq.gz
# foo_2_treat_R2.fastq.gz
# foo_3_treat_R1.fastq.gz
# foo_3_treat_R2.fastq.gz
# foo_4_control_R1.fastq.gz
# foo_4_control_R2.fastq.gz
# foo_5_control_R1.fastq.gz
# foo_5_control_R2.fastq.gz
# foo_6_control_R1.fastq.gz
# foo_6_control_R2.fastq.gz
# R1/R2 corresponds to paired reads in PE_100 raw data

# set variable for working directory
# CHANGE AS NEEDED
workingdirectory=("/mnt/usr/dir")

# set variable for software directories and export to $PATH
# CHANGE AS NEEDED
bin=("${workingdirectory}/bin")
trimmomatic=("${bin}/Trimmomatic-0.36")
tophat=("${bin}/tophat-2.1.1.Linux_x86_64")
cufflinks=("${bin}/cufflinks-2.2.1.Linux_x86_64")
samtools=("${bin}/samtools-1.3.1")
export PATH=${trimmomatic}:$PATH
export PATH=${tophat}:$PATH
export PATH=${samtools}:$PATH
export PATH=${cufflinks}:$PATH

# create variable for experimental comparison descriptor and date, i.e. "4mo_Het-vs-WT"
# CHANGE AS NEEDED
Descriptor=("foo")
ConditionOne=("treat")
ConditionTwo=("control")
experiment=("${Descriptor}_${ConditionOne}-vs-${ConditionTwo}")

# set variable for raw reads directory
# CHANGE AS NEEDED
rawreadsfolder=("${workingdirectory}/raw_reads")

###############################################################################################
# SCRIPT START
###############################################################################################

# make variable for raw data
cd ${rawreadsfolder}
raw_data=$(ls *.fastq.gz)

# make temporary file directory
mkdir ${rawreadsfolder}/temp_files
temp=("${rawreadsfolder}/temp_files")

# move to working directory
cd ${workingdirectory}

# make file containing all raw data file names for given descriptor
for file in $raw_data; do
	file_cut_f1() {
		echo "${file}" | cut -d'_' -f1;
	}
	f1_var=$(file_cut_f1);
	if [ ${f1_var} = ${Descriptor} ]; then
		echo "${file}" >> ${temp}/${Descriptor}.txt;
	fi
done

# make separate files containing raw data file names for each 'condition' under given descriptor
for line in $(cat ${temp}/${Descriptor}.txt); do
	line_cut_f3() {
		echo "${line}" | cut -d'_' -f3;
	}
	f3_var=$(line_cut_f3);
	if [ ${f3_var} = ${ConditionOne} ]; then
		echo "${line}" >> ${temp}/${Descriptor}_${ConditionOne}.txt; 
	elif [ ${f3_var} = ${ConditionTwo} ]; then
		echo "${line}" >> ${temp}/${Descriptor}_${ConditionTwo}.txt;
	fi
done

# create variables for each 'condition'
groupAfile=("${temp}/${Descriptor}_${ConditionOne}.txt")
groupBfile=("${temp}/${Descriptor}_${ConditionTwo}.txt")
for line in $(cat $groupAfile); do
	groupA=("${groupA} ${line}");
done
for line in $(cat $groupBfile); do
	groupB=("${groupB} ${line}");
done

# strip excess information to leave only 'descriptor', 'sample,' and 'condition'
for i in $groupA; do
	short_string() {
		echo "${i}" | cut -d'_' -f1-3;
	}
	shortgroupA=("${shortgroupA} $(short_string)");
done
for i in $groupB; do
	short_string() {
		echo "${i}" | cut -d'_' -f1-3;
	}
	shortgroupB=("${shortgroupB} $(short_string)");
done

# remove duplicate values
finalgroupA=$(echo "${shortgroupA}" | tr ' ' '\n' | nl | sort -u -k2 | sort -n | cut -f2-)
finalgroupB=$(echo "${shortgroupB}" | tr ' ' '\n' | nl | sort -u -k2 | sort -n | cut -f2-)

# create array from each group
groupAarray=($finalgroupA)
groupBarray=($finalgroupB)

# create combined groupA+B variable
finalsamples=("${finalgroupA} ${finalgroupB}")

# make directory for trimmomatic output
mkdir ${workingdirectory}/trimmed_reads
mkdir ${workingdirectory}/trimmed_reads/${experiment}
#variable for experimental trimmed_reads directory
trimmed=("${workingdirectory}/trimmed_reads/${experiment}")

#########################################
# trimmomatic
# PE_100
for i in $finalsamples; do
java -jar ${trimmomatic}/trimmomatic-0.36.jar PE -phred33 \
${rawreadsfolder}/${i}_R1.fastq.gz ${rawreadsfolder}/${i}_R2.fastq.gz \
${trimmed}/trimmed_${i}_R1_paired.fastq.gz ${trimmed}/trimmed_${i}_R1_unpaired.fastq.gz \
${trimmed}/trimmed_${i}_R2_paired.fastq.gz ${trimmed}/trimmed_${i}_R2_unpaired.fastq.gz \
MINLEN:50 CROP:75 HEADCROP:15;
done

# run fastqc on newly trimmed files
for i in $finalsamples; do
fastqc ${trimmed}/trimmed_${i}_R1_paired.fastq.gz;
fastqc ${trimmed}/trimmed_${i}_R2_paired.fastq.gz;
done

# make directory for tophat output
mkdir ${workingdirectory}/thout
mkdir ${workingdirectory}/thout/${experiment}
# variable for experimental thout directory
thout=("${workingdirectory}/thout/${experiment}")

###########################################
# tophat2
for i in $finalsamples; do
${tophat}/tophat2 -p 8 -G genes.gtf -o ${thout}/${i}_thout genome \
${trimmed}/trimmed_${i}_R1_paired.fastq.gz ${trimmed}/trimmed_${i}_R2_paired.fastq.gz;
done

# sort tophat output accepted_hits.bam by name
for i in $finalsamples; do
${samtools}/samtools sort -o ${thout}/${i}_thout/accepted_hits.sorted.bam -n -@ 8 ${thout}/${i}_thout/accepted_hits.bam;
done

# make directory for HTSeq output
mkdir ${workingdirectory}/htseq_output
mkdir ${workingdirectory}/htseq_output/$experiment}
htseqout=("${workingdirectory}/htseq_output/$experiment}")

###########################################
# count reads with HTSeq
for i in $finalsamples; do
htseq-count -f bam  ${thout}/${i}_thout/accepted_hits.sorted.bam genes.gtf > ${htseqout}/${i}_counts.txt;
done

# make directory for cufflinks output
mkdir ${workingdirectory}/clout
mkdir ${workingdirectory}/clout/${experiment}
# variable for experimental clout directory
clout=("${workingdirectory}/clout/${experiment}")

###########################################
# cufflinks
for i in $finalsamples; do
${cufflinks}/cufflinks -p 8 -o ${clout}/${i}_clout -g genes.gtf ${thout}/${i}_thout/accepted_hits.bam;
done

# create assemblies text file for cuffmerge
for i in $finalsamples; do
echo "${clout}/${i}_clout/transcripts.gtf" >> assemblies_${experiment}.txt;
done

###########################################
# cuffmerge
${cufflinks}/cuffmerge -g genes.gtf -s genome.fa -p 8 assemblies_${experiment}.txt

# rename cuffmerge output folder and move to merged directory
mv merged_asm merged_asm_${experiment}
mkdir merged
mv merged_asm_${experiment} merged
# variable for merged.gtf file
merged_asm=("merged/merged_asm_${experiment}/merged.gtf")

# create variable for comma-separated locations for cuffdiff input
for sample in "${groupAarray[@]}"; do
	var=("${thout}/${sample}_thout/accepted_hits.bam");
	if [ "$groupAcuffdiffvar" = "" ]; then
		groupAcuffdiffvar=("$var");
	else
		groupAcuffdiffvar=("$groupAcuffdiffvar,$var");
	fi
done
for sample in "${groupBarray[@]}"; do
	var=("${thout}/${sample}_thout/accepted_hits.bam");
	if [ "$groupBcuffdiffvar" = "" ]; then
		groupBcuffdiffvar=("$var");
	else
		groupBcuffdiffvar=("$groupBcuffdiffvar,$var");
	fi
done

# make directory for cuffdiff output
mkdir ${workingdirectory}/diff_out
mkdir ${workingdirectory}/diff_out/diff_out_${experiment}
# variable for experimental diff_out directory
diff_out=("${workingdirectory}/diff_out/diff_out_${experiment}")

###########################################
# cuffdiff
${cufflinks}/cuffdiff -o ${diff_out} -b genome.fa -p 8 -L \
${ConditionOne},${ConditionTwo} ${merged_asm} \
${groupAcuffdiffvar} ${groupBcuffdiffvar}

# make directory for cuffnorm output
mkdir ${workingdirectory}/cuffnorm_out
mkdir ${workingdirectory}/cuffnorm_out/cuffnorm_out_${experiment}
# variable for experimental diff_out directory
cuffnorm_out=("${workingdirectory}/cuffnorm_out/cuffnorm_out_${experiment}")

###########################################
# cuffnorm
${cufflinks}/cuffnorm -o ${cuffnorm_out} -p 8 -L \
${ConditionOne},${ConditionTwo} ${merged_asm} \
${groupAcuffdiffvar} ${groupBcuffdiffvar}

# remove temp files and directory
rm -r ${temp}
