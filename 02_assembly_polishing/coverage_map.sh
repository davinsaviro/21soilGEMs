#!/bin/bash

######### SLURM OPTIONS
#SBATCH --job-name read_depth
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 6G
#SBATCH --time 00:10:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/read_map/002_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/read_map/002_%j.err

######### COMMANDS TO LOAD MODULES

user=dwijaya
organism=002
home_dir=/users/${user}/scratch/${user}
assembly_dir=${home_dir}/Oxford_Nanopore/read_map
illumina_dir=${home_dir}/Illumina_Reads/trimmed

source /dcsrsoft/spack/bin/setup_dcsrsoft

module load gcc
module load bwa
module load samtools

bwa index ${assembly_dir}/concat_${organism}.fasta
bwa mem ${assembly_dir}/concat_${organism}.fasta ${illumina_dir}/${organism}_*_R1_*_paired.fastq.gz ${illumina_dir}/${organism}_*_R2_*_paired.fastq.gz > ${assembly_dir}/${organism}/${organism}.sam
samtools view -b ${assembly_dir}/${organism}/${organism}.sam > ${assembly_dir}/${organism}/${organism}.bam
samtools sort ${assembly_dir}/${organism}/${organism}.bam > ${assembly_dir}/${organism}/${organism}_sorted.bam
samtools depth -a ${assembly_dir}/${organism}/${organism}_sorted.bam > ${assembly_dir}/${organism}/${organism}_depth.txt
