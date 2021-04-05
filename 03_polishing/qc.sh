#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name assembly_QC
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 8G
#SBATCH --time 01:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/polishing/QC/002_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/polishing/QC/002_%j.err

user=dwijaya
sample=002004
organism=002
home_dir=/users/${user}/scratch/${user}
work_dir=${home_dir}/Oxford_Nanopore
out_dir=${work_dir}/polishing/QC/${organism}

module load HPC/Software 
module load UHTS/Analysis/samtools/1.8
module load UHTS/Quality_control/quast/4.6.0
module load UHTS/Analysis/mummer/3.9.4alpha

start=$SECONDS
echo "Step: running QC"

echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${organism}

###===========================
##Quast
echo -e "-------1. First run Quast"
###===========================

mkdir -p ${work_dir}/polishing/QC/${organism}/quast

python3.6 /software/UHTS/Quality_control/quast/4.6.0/bin/quast.py \
      -l "FlyeAssembly,first_ONT_racon,second_ONT_racon,first_Illumina_pilon,second_Illumina_pilon" \
      -R ${work_dir}/polishing/Illumina/${organism}/05_third_pillon_bowtie2/assembly/${sample}_third_Pilon_Assembly.fasta \
      ${work_dir}/assembly/${sample}/${sample}assembly.fasta \
      ${work_dir}/polishing/ONT/${sample}/01_firstRacon_graphmap/assembly/${sample}_first_Racon_Assembly.fasta \
      ${work_dir}/polishing/ONT/${sample}/02_secondRacon_graphmap/assembly/${sample}_second_Racon_Assembly.fasta \
      ${work_dir}/polishing/Illumina/${organism}/03_first_pillon_bowtie2/assembly/${sample}_first_Pilon_Assembly.fasta \
      ${work_dir}/polishing/Illumina/${organism}/04_second_pillon_bowtie2/assembly/${sample}_second_Pilon_Assembly.fasta \
      -o ${out_dir}/quast/

###===========================
##check mapping
echo -e "-------2. Second check mapping"
###===========================   

mkdir -p ${out_dir}/mapping

samtools flagstat ${work_dir}/polishing/Illumina/${organism}/04_second_pillon_bowtie2/bowtie2/second_pillon_ONT2FinalAssembly_sorted.bam >  ${out_dir}/mapping/${organism}_mapping.txt

grep "mapped (" ${out_dir}/mapping/${organism}_mapping.txt >> ${out_dir}/mapping/Final_mapping.txt

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
