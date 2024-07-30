# sanger-sequencing
Pipeline for analyzing 16S Sanger sequencing data.
## Before running the pipeline

To run this pipeline, you will need to download some dependencies and databases first. You only need to do this once, unless you want to update the software or the databases. If you use it in a HPC cluster, some programs may be already installed. Dependencies are:
- [fastp](https://github.com/OpenGene/fastp)
- [PEAR](https://cme.h-its.org/exelixis/web/software/pear/index.html)
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). User installation of this one is recommended as it is updated frequently and compatibility problems with databases may arise.

### Install blast_plus

Download the last version of BLAST:

`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST`

Unzip the file:

`tar zxvpf ncbi-blast-2.14.1+-x64-linux.tar.gz`

Now you will get a directory called 'ncbi-blast-2.14.1+/'. The 'bin' directory is inside and contains all the programs.

Put BLAST plus in your PATH in the bashrc file. If you don't have a bashrc file, go to your home and create it:

`vi .bashrc`

And add BLAST plus (write the complete path to 'bin'; you can use 'pwd' command in 'bin' to get it):

`export PATH=$PATH:path/to/ncbi-blast-2.14.1+/bin`

Exit and login again. Now you can use all the programs included in BLAST plus.
Make sure that you don't have more than one version installed. Remove older versions of your path if you install a new one. Also do not load the program as a module if you are using your own installation.

### Download 16S blast DB - VERSION 5 (v5) - IMPORTANT! It doesn't work with v4

For its use with BLASTN, databases need to be formatted using blast commands. To create customized databases, use the script 'makeblastdb.sh'. Edit the script using a plain text editor and following the instructions inside. If you are using a HPC cluster, **NEVER** run it locally, use SLURM:

`sbatch makeblastdb.sh`

## Running the pipeline

The bash script '16s_sanger.sh' will run all the pipeline following the steps below.

Open the script and edit it with a plain text editor. Change parameters if you wish and paths to database and raw data. Remember to always provide the **full** path.

Be careful with naming convention. This pipeline is adapted to STABVIDA file naming. Sometimes, the company provides files with different names. Try to adapt it if possible.

If you are using a HPC cluster, do not run the script locally, use SLURM:

`sbatch 16s_sanger_sequencing.sh`

The steps included in the pipeline are:

### 1. seqQual2FASTQ

Perl script that merges seq file and qual file into FASTQ file, which is the standard file used as input for most of the programs used to analyze sequencing data. FASTQ files will be stored at '001_merged_fastq' folder.

### 2. Trimming, merging and QC

In this step, quality checks and trimming will be performed for all reads, single or paired. Paired-end reads will also be merged into a contig. QC reports will be stored at '002_qc_reports' folder and contigs will be processed at '003_contigs'. In this folder, you will find merged and unmerged reads. If the assembly of the reads was successful and a contig was obtained, the full contig will be used for further steps. If not, the forward and reverse reads will be used. If the quality of forward and/or reverse reads is too low, the reads will be discarded and the analysis of the sample will not be performed. Contigs and single reads will be stored at '004_fasta/' and low quality (empty) reads can be found at: '004_fasta/empty_low_quality/'.

### 3. FASTQ to multiFASTA

FASTQ files will be converted to FASTA files and then all sequences will be merged in a multiFASTA file that will be the input for BLAST.

### 5. BLAST

MultiFASTA will be blasted to the database you provided. You will find three outputs at '005_blastn': (1) 'all_results.bls', which contains all the results produced by BLAST, (2) '1match_table.bls', which contains the original ID, the first match and the percentage of identity; and (3) '5matches_table.bls', which contains the same as the previous but for five matches for each sequence. If you use SILVA database, you will also find a simplified table with the ID and the taxonomy.
