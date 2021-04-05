#!/bin/bash
#SBATCH --job-name FlyeAssembly
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 12G
#SBATCH --time 03:00:00
#SBATCH --export NONE
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/assembly/002004.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/assembly/002004.err

sample=002004
sfile=20181119_1648_002-004
username=dwijaya
work_dir=/users/dwijaya/scratch/dwijaya/Oxford_Nanopore/assembly/${sample}
sample_dir=/users/dwijaya/scratch/dwijaya/Oxford_Nanopore/filter/${sample}

GENOME_SIZE=12.2m

start=$SECONDS

module load HPC/Software 
module load UHTS/Assembler/flye/2.7.1

echo "=========================================================="
date +"START : %a %b %e %H:%M:%S %Z %Y"
echo -e "Sample Name: "${sample}
echo "Step: Genome assembly"

flye --threads 8  --iterations 5 --genome-size ${GENOME_SIZE} --nano-raw \
      ${sample_dir}/After_filtlong_trimming/${sfile}_filtered.fastq \
      --out-dir ${work_dir}/${sfile}

date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
