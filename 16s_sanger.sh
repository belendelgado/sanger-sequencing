#!/usr/bin/env bash
#SBATCH -J 16s_sanger.sh
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --constraint=cal
#SBATCH --error=16s_sanger.%J.err
#SBATCH --output=16s_sanger.%J.out

##############################################################
######################### PARAMETERS #########################
##############################################################

# IMPORTANT: merge2fastq.pl script needed

## Quality score ##
# less restrictive (20) or more restrictive (30)
export phred="20"

## Minimum overlap length for the contigs assembly ##
# This is the minimum overlap for reads to be merged when you have paired reads for the same 16S
export overlap="12"

## 16S database for blastn
# path/to/database/dbName
# blastn recognizes all the suffixes of the db, just provide the basename
# this one will work because we have reading permissions within the group, but feel free to use your own
export database="path/to/db/basename"

# Path to raw data (containing .seq, .qual, .raw.seq...)
export datadir="path/to/rawData"

##############################################################
##############################################################
######################### DON'T TOUCH ########################
##############################################################
##############################################################


##############################################################
###################### 1. seqQual2Fastq ######################
##############################################################

# Create directory for seq and qual files
mkdir -p 000_seqQual

# Copy only seq and qual from raw data
cp ${datadir}/*.seq 000_seqQual
cp ${datadir}/*.qual 000_seqQual
rm 000_seqQual/*raw.seq

SEQUAL_DIR="000_seqQual"

for PRIMER in `ls $SEQUAL_DIR/*GAGTT*`; do
    rename '_5_-GAGTTTGATGCTGGCTCAG_-3_' "" $PRIMER
done

for PRIMER in `ls $SEQUAL_DIR/*ACCTT*`; do
    rename '_5_-ACCTTGTTACGACTT-3_' "" $PRIMER
done

# Create a directory for fastq files
mkdir -p 001_merged_fastq

# Merge seq and qual to create a fastq file
# IMPORTANT: merge2fastq.pl script needed


############## ADDED 23-07 #################
cd 000_seqQual
for file in *.qual *.seq; do
    mv "$file" "${file#*+}"
done
cd ..
############################################


for SEQ in `ls 000_seqQual/*.seq`
do
    base_name=$(basename "$SEQ" .seq)
    perl merge2fastq.pl 000_seqQual/${base_name}.seq 000_seqQual/${base_name}.qual > "001_merged_fastq"/${base_name}.fastq
done

##############################################################
################ 2. Trimming, QC and merging #################
##############################################################

module load fastp/0.23.4
module load pear/0.9.6

# Create necessary directories
mkdir -p 003_contigs
mkdir -p 003_contigs/clean_pairs
mkdir -p 002_qc_reports
mkdir -p 004_fasta

# Directory containing reads
READS_DIR="001_merged_fastq"

# List reads
READS=$(ls $READS_DIR)

for READ in $READS; do
    # Extract sample names
     SAMPLE_NAME=$(echo $READ | sed -E 's/\+27F.fastq|\+1492R.fastq//g')
     echo "$SAMPLE_NAME" >> dup_sample_names.txt
done

sort dup_sample_names.txt | uniq > sample_names.txt
rm dup_sample_names.txt

while read SAMPLE_NAME; do
    FORWARD_READ="${READS_DIR}/${SAMPLE_NAME}+27F.fastq"
    REVERSE_READ="${READS_DIR}/${SAMPLE_NAME}+1492R.fastq"

    # Check if both forward and reverse reads exist
    if [[ -e $FORWARD_READ && -e $REVERSE_READ ]]; then
        echo "Forward and reverse reads found for sample ${SAMPLE_NAME}. Running trimming and assembly."
        fastp --in1 $FORWARD_READ --in2 $REVERSE_READ \
              --cut_front 50 --cut_tail 50 --cut_mean_quality $phred \
              --html 002_qc_reports/${SAMPLE_NAME}.html \
              --out1 003_contigs/clean_pairs/${SAMPLE_NAME}_1.fastq \
              --out2 003_contigs/clean_pairs/${SAMPLE_NAME}_2.fastq

        echo "Merging forward and reverse reads"
        pear -f 003_contigs/clean_pairs/${SAMPLE_NAME}_1.fastq \
             -r 003_contigs/clean_pairs/${SAMPLE_NAME}_2.fastq \
             -o 003_contigs/${SAMPLE_NAME} -v $overlap

        ASSEMBLED="003_contigs/${SAMPLE_NAME}.assembled.fastq"
        UNASSEMBLED_F="003_contigs/${SAMPLE_NAME}.unassembled.forward.fastq"
        UNASSEMBLED_R="003_contigs/${SAMPLE_NAME}.unassembled.reverse.fastq"

        if [[ -s $ASSEMBLED ]]; then
            echo "Using assembled reads for sample ${SAMPLE_NAME}."
            sed -n '1~4s/^@/>/p;2~4p' $ASSEMBLED > 004_fasta/${SAMPLE_NAME}.fasta
        else
            echo "Assembled reads are empty for sample ${SAMPLE_NAME}. Using unassembled reads."
            sed -n '1~4s/^@/>/p;2~4p' $UNASSEMBLED_F > 004_fasta/${SAMPLE_NAME}_F.fasta
            sed -n '1~4s/^@/>/p;2~4p' $UNASSEMBLED_R > 004_fasta/${SAMPLE_NAME}_R.fasta
        fi

    elif [[ -e $FORWARD_READ ]]; then
        echo "Forward read found for ${SAMPLE_NAME}. Running trimming."
        fastp --in1 $FORWARD_READ \
              --cut_front 50 --cut_tail 50 --cut_mean_quality $phred \
              --html 002_qc_reports/${SAMPLE_NAME}.html \
              --out1 003_contigs/${SAMPLE_NAME}_F.fastq
        sed -n '1~4s/^@/>/p;2~4p' 003_contigs/${SAMPLE_NAME}_F.fastq > 004_fasta/${SAMPLE_NAME}_F.fasta

    elif [[ -e $REVERSE_READ ]]; then
        echo "Reverse read found for ${SAMPLE_NAME}. Running trimming."
        fastp --in1 $REVERSE_READ \
              --cut_front 50 --cut_tail 50 --cut_mean_quality $phred \
              --html 002_qc_reports/${SAMPLE_NAME}.html \
              --out1 003_contigs/${SAMPLE_NAME}_R.fastq
        sed -n '1~4s/^@/>/p;2~4p' 003_contigs/${SAMPLE_NAME}_R.fastq > 004_fasta/${SAMPLE_NAME}_R.fasta

    else
        echo "No reads found for sample ${SAMPLE_NAME}."
    fi
done < sample_names.txt

#find 004_fasta -type f -empty -delete

mkdir -p 004_fasta/empty_low_quality
find 004_fasta -type f -empty -exec mv {} 004_fasta/empty_low_quality/ \;

# rename fasta headers
for file in 004_fasta/*.fasta; do
    filename=$(basename "$file" .fasta)
    awk -v name="$filename" '/^>/{print ">"name; next}{print}' "$file" > temp.fasta
    mv temp.fasta "$file"
done

############################################################
######################## 3. BLASTN #########################
############################################################

# MultiFASTA with all the sequences to be blasted

mkdir -p 005_blast
touch 005_blast/all_sequences.fasta
cat 004_fasta/*fasta >> "005_blast"/all_sequences.fasta

multifasta="005_blast/all_sequences.fasta"

# Run blast
time blastn -db $database -query $multifasta -num_threads 8 -out 005_blast/all_results.bls
time blastn -db $database -query $multifasta -num_threads 8 \
            -max_target_seqs 1 -outfmt "6 qseqid sseqid pident" \
            -out 005_blast/1match_table.bls
time blastn -db $database -query $multifasta -num_threads 8 \
            -max_target_seqs 5 -outfmt "6 qseqid sseqid pident" \
            -out 005_blast/5matches_table.bls

############################################################


