#!/usr/bin/env bash
#SBATCH -J makeblastdb.sh
#SBATCH --cpus-per-task=16
#SBATCH --mem=800gb
#SBATCH --time=2:00:00
#SBATCH --constraint=cal
#SBATCH --error=makeblastdb.%J.err
#SBATCH --output=makeblastdb.%J.out

input_fasta="path/to/input.fasta"
db_title="title"

time makeblastdb -in $input_fasta -blastdb_version 5 -title $db_title -dbtype nucl
