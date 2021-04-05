#!/bin/bash
#SBATCH --job-name trim
#SBATCH --partition wally
#SBATCH --time 20:00
#SBATCH --mem 4G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --output /users/dwijaya/scratch/dwijaya/Illumina_Reads/trimmed/002_S1.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Illumina_Reads/trimmed/002_S1.err

source /dcsrsoft/spack/bin/setup_dcsrsoft

module load gcc
module load trimmomatic 

sample_name=002_S1
username=dwijaya
working_directory=/users/${username}/scratch/${username}/Illumina_Reads/trimmed/
sample_directory=/users/${username}/scratch/${username}/Illumina_Reads/

trimmomatic PE -threads 8 ${sample_directory}${sample_name}_L001_R1_001.fastq.gz ${sample_directory}${sample_name}_L001_R2_001.fastq.gz \
${working_directory}${sample_name}_L001_R1_001_paired.fastq.gz ${working_directory}${sample_name}_L001_R1_001_unpaired.fastq.gz \
${working_directory}${sample_name}_L001_R2_001_paired.fastq.gz ${working_directory}${sample_name}_L001_R2_001_unpaired.fastq.gz \
ILLUMINACLIP:${sample_directory}AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 \
SLIDINGWINDOW:4:15 MINLEN:60 
