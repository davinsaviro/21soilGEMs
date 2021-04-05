#!/bin/bash
#SBATCH --job-name qc
#SBATCH --partition wally
#SBATCH --time 20:00
#SBATCH --mem 4G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --output /users/dwijaya/scratch/dwijaya/Illumina_Reads/002_S1/%x_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Illumina_Reads/002_S1/%x_%j.err

source /dcsrsoft/spack/bin/setup_dcsrsoft

sample_name=002_S1
username=dwijaya
working_directory=/users/${username}/scratch/${username}/Illumina_Reads/${sample_name}/
sample_directory=/users/${username}/scratch/${username}/Illumina_Reads/

module load gcc 
module load fastqc 

fastqc -o ${working_directory} ${sample_directory}${sample_name}_L001_R1_001.fastq.gz ${sample_directory}${sample_name}_L001_R2_001.fastq.gz
