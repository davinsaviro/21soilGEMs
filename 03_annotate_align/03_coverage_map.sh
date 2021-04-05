#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name read_mapping_complete_genome
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 4G
#SBATCH --time 03:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/prokka/mapping/002_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/prokka/mapping/002_%j.err

user=dwijaya
sample=002
home_dir=/users/${user}/scratch/${user}
ref_fasta=${home_dir}/Oxford_Nanopore/prokka/start_align/${sample}/Final_${sample}_startaligned.fasta
out_dir=${home_dir}/Oxford_Nanopore/prokka/mapping/${sample}
illumina_dir=${home_dir}/Illumina_Reads/trimmed

module load HPC/Software
module load Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.17
module load UHTS/Analysis/samtools/1.8
module load Development/gcc/9.2.1

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: read mapping against complete genome"

#done
echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${sample}


###===========================
##calculating depth of the coverage at every location
echo -e "-------1. calculating depth of the coverage at every location"
###===========================

mkdir -p ${home_dir}/Oxford_Nanopore/prokka/mapping/${sample}
echo ">concat_${sample}" > ${out_dir}/concat_${sample}.fasta
cat ${ref_fasta} | grep -v '>' >> ${out_dir}/concat_${sample}.fasta
ref_fasta=${out_dir}/concat_${sample}.fasta

bwa index ${ref_fasta}
bwa mem ${ref_fasta} ${illumina_dir}/${sample}_*_R1_*_paired.fastq.gz ${illumina_dir}/${sample}_*_R2_*_paired.fastq.gz > ${out_dir}/${sample}_complete.sam

samtools view -b ${out_dir}/${sample}_complete.sam > ${out_dir}/${sample}_complete.bam
samtools sort ${out_dir}/${sample}_complete.bam > ${out_dir}/${sample}_complete_sorted.bam
samtools depth -a ${out_dir}/${sample}_complete_sorted.bam > ${out_dir}/${sample}_complete_depth.txt

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
