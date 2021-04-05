#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name readLenth_rawReads
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 6G
#SBATCH --time 1:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/lengthraw/S002004/%x_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/lengthraw/S002004/%x_%j.err

####--------------------------------------
##preparation
##set you bash variables in order to quickly call them in the script
####--------------------------------------

username=dwijaya
personal_home_directory=/users/${username}/scratch/${username}
raw_data_directory=/scratch/wally/FAC/FBM/DMF/jvanderm/soil_metagenome/Soil_genomes/Nanopore_03_reads
sample_name=20181119_1648_002-004

####--------------------------------------
##Introduction to script
####--------------------------------------

start=$SECONDS
echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${sample_name}
echo "Step: Read length calculation"

####--------------------------------------
##Calculate raw read statistics
echo -e "1. First calculate the read length"
####--------------------------------------

mkdir -p ${personal_home_directory}/Oxford_Nanopore/lengthraw/S002004/log/

awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'  ${raw_data_directory}/${sample_name}.lnt_sel.fastq > ${personal_home_directory}/Oxford_Nanopore/lengthraw/S002004/log/${sample_name}_raw_meanlength.txt

less ${raw_data_directory}/${sample_name}.lnt_sel.fastq | awk '{if(NR%4==2) print length($1)}' > ${personal_home_directory}/Oxford_Nanopore/lengthraw/S002004/log/${sample_name}_raw_readLength.txt

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
