#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name pilon_polishing
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 20G
#SBATCH --time 01:30:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/polishing/Illumina/002_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/polishing/Illumina/002_%j.err

user=dwijaya
sample=002004
organism=002
threads=8
home_dir=/users/${user}/scratch/${user}
work_dir=/users/${user}/scratch/${user}/Oxford_Nanopore
out_dir=${work_dir}/polishing/Illumina/${organism}

reference_fasta=${work_dir}/polishing/ONT/${sample}/02_secondRacon_graphmap/assembly/${sample}_second_Racon_Assembly.fasta
raw_data_dir_illumina=${home_dir}/Illumina_Reads/trimmed

####--------------------------------------
##modules
####--------------------------------------

module load HPC/Software 
module load UHTS/Analysis/pilon/1.22
module load UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Analysis/samtools/1.8

####--------------------------------------
##start of script
####--------------------------------------

start=${SECONDS}
echo "Step: Pilon polishing"

echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${organism}

###===========================
##Priming the variables used and starting the for-loop
echo -e "-------0. Priming the variables used and starting the for-loop"
###===========================   

overall_polishing_counter=3

for Illumina_polishing_counter in $(echo "first second third")
do 

echo -e "----"${Illumina_polishing_counter}" round of pilon polishing"

###===========================
##bowtie2 mapping
echo -e "-------1. First we map with bowtie2"
###===========================      
mkdir -p ${out_dir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/{bowtie2,assembly}

bowtie2-build --quiet $reference_fasta $reference_fasta

bowtie2 -1 ${raw_data_dir_illumina}/${organism}_*_R1*_paired.fastq.gz \
        -2 ${raw_data_dir_illumina}/${organism}_*_R2*_paired.fastq.gz \
        -x ${reference_fasta} \
        --threads ${threads} \
        --local --very-sensitive-local | samtools sort -O BAM -o ${out_dir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/bowtie2/${Illumina_polishing_counter}_pillon_ONT2FinalAssembly_sorted.bam - 
     
###--------------------
##index bam file
###--------------------     
samtools index ${out_dir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/bowtie2/${Illumina_polishing_counter}_pillon_ONT2FinalAssembly_sorted.bam

###===========================
##pilon polishing
echo -e "-------2. Second polish with pilon"
###===========================   

pilon --threads ${threads} --genome $reference_fasta \
  --frags ${out_dir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/bowtie2/${Illumina_polishing_counter}_pillon_ONT2FinalAssembly_sorted.bam \
  --output ${sample}_${Illumina_polishing_counter}_Pilon_Assembly \
  --outdir ${out_dir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/assembly/ \
  --changes

###===========================
##reset the parameters
echo -e "-------3. Resetting the parameters"
###===========================   

reference_fasta=$(echo ${out_dir}/0${overall_polishing_counter}_${Illumina_polishing_counter}_pillon_bowtie2/assembly/${sample}_${Illumina_polishing_counter}_Pilon_Assembly.fasta)
echo $reference_fasta
overall_polishing_counter=$((overall_polishing_counter+1)) 

done #polishing round

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
