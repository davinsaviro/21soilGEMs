#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name racon_polishing
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 8G
#SBATCH --time 04:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/polishing/ONT/002004.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/polishing/ONT/002004.err

####--------------------------------------
##preparation, set variables
####--------------------------------------

user=dwijaya
sample=002004
sfile=20181119_1648_002-004
work_dir=/users/${user}/scratch/${user}/Oxford_Nanopore
sample_dir=${work_dir}/assembly/${sample}
out_dir=${work_dir}/polishing/ONT/${sample}

reference_fasta=${sample_dir}/${sample}assembly.fasta
ont_reads=${work_dir}/filter/${sample}/After_filtlong_trimming/${sfile}_filtered.fastq

####--------------------------------------
##modules
####--------------------------------------

module load HPC/Software 
module load SequenceAnalysis/SequenceAlignment/graphmap/0.5.2
module load UHTS/Assembler/racon/1.0.1

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: Racon polishing"

echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${sample}

###===========================
##Priming the variables used and starting the for-loop
echo -e "-------0. Priming the variables used and starting the for-loop"
###===========================   

overall_polishing_counter=1

for ONT_polishing_counter in $(echo "first second")
do 
echo -e "----"${ONT_polishing_counter}" round of Racon polishing"

###===========================
##graphmap mapping
echo -e "-------1. First we map with graphmap"
###===========================      

mkdir -p ${out_dir}/0${overall_polishing_counter}_${ONT_polishing_counter}Racon_graphmap/{Graphmap_ONT,assembly}

graphmap align --rebuild-index --circular  \
  -r ${reference_fasta} \
  -d ${ont_reads} \
  -o ${out_dir}/0${overall_polishing_counter}_${ONT_polishing_counter}Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}Racon_ONT2FinalAssembly_sorted.sam 

###===========================
##racon polishing
echo -e "-------2. Second polish with Racon"
###===========================   

racon ${ont_reads} \
  ${out_dir}/0${overall_polishing_counter}_${ONT_polishing_counter}Racon_graphmap/Graphmap_ONT/${ONT_polishing_counter}Racon_ONT2FinalAssembly_sorted.sam \
  ${reference_fasta} > \
  ${out_dir}/0${overall_polishing_counter}_${ONT_polishing_counter}Racon_graphmap/assembly/${sample}_${ONT_polishing_counter}_Racon_Assembly.fasta

###===========================
##reset the parameters
echo -e "-------3. Resetting the parameters"
###===========================   

reference_fasta=$(echo ${out_dir}/0${overall_polishing_counter}_${ONT_polishing_counter}Racon_graphmap/assembly/${sample}_${ONT_polishing_counter}_Racon_Assembly.fasta)
overall_polishing_counter=$((overall_polishing_counter+1)) 

done #polishing round


####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
