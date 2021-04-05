#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name readFiltering
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 2G
#SBATCH --time 2:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/filter/002004.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/filter/002004.err

username=dwijaya
personal_home_directory=/users/${username}/scratch/${username}/Oxford_Nanopore/filter/002004
raw_data_directory=/scratch/wally/FAC/FBM/DMF/jvanderm/soil_metagenome/Soil_genomes/Nanopore_03_reads
sample_name=20181119_1648_002-004

##!!!!!!!here you can adjust the parameters
##see how these bash variables are used in the script further bellow
MINIMUM_read_LENGTH=7000
min_mean_q=10
length_weight=10
target_bases=1230213300

####--------------------------------------
##modules
####--------------------------------------

source /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc
module load filtlong

####--------------------------------------
##Introduction to script
####--------------------------------------

start=$SECONDS
echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${sample_name}
echo "Step: Read filtering"

####--------------------------------------
##Filter reads
echo -e "1. First we filter the reads"
####--------------------------------------

mkdir -p ${personal_home_directory}/After_filtlong_trimming/
echo "=========================================================="
echo -e ${sample_name}

filtlong --min_length ${MINIMUM_read_LENGTH} \
        --min_mean_q ${min_mean_q} \
        --length_weight ${length_weight} \
        --target_bases ${target_bases}  \
        ${raw_data_directory}/${sample_name}.lnt_sel.fastq > \
        ${personal_home_directory}/After_filtlong_trimming/${sample_name}_filtered.fastq

####--------------------------------------
##Read length stats
echo -e "2. Secondly we calculate the read length distribution again"
####--------------------------------------

mkdir -p ${personal_home_directory}/log/

awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'  ${personal_home_directory}/After_filtlong_trimming/${sample_name}_filtered.fastq > ${personal_home_directory}/log/${sample_name}_filtered_meanlength.txt

less ${personal_home_directory}/After_filtlong_trimming/${sample_name}_filtered.fastq | awk '{if(NR%4==2) print length($1)}' > ${personal_home_directory}/log/${sample_name}_filtered_readLength.txt

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
