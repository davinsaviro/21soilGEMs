#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name startAlign
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 12G
#SBATCH --time 03:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/prokka/start_align/002_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/prokka/start_align/002_%j.err

####--------------------------------------
##preparation, setting bash variables
####--------------------------------------

user=dwijaya
sample=002
threads=4
home_dir=/users/${user}/scratch/${user}
out_dir=${home_dir}/Oxford_Nanopore/prokka/start_align/${sample}

ref_fasta=${home_dir}/Oxford_Nanopore/prokka/start_align/${sample}assembly.fasta
raw_illumina=${home_dir}/Illumina_Reads/trimmed

bac_contig_nam=contig_2_C
Startalign_location=3396129
####--------------------------------------
##modules
####--------------------------------------

module load HPC/Software 
module load UHTS/Analysis/pilon/1.22
module load UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Analysis/samtools/1.8
module load Emboss/EMBOSS/6.6.0
module load UHTS/Analysis/BEDTools/2.29.2

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: startalign Reverse strand"

echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${sample}

###===========================
##seperating bacterial contig
echo -e "-------1. Extracting the correct contig"
###===========================   
mkdir -p ${out_dir}/tmp

samtools faidx ${ref_fasta} ${bac_contig_nam} > ${out_dir}/only_bacteria_contig.fasta

###--------------------
##extract remaining contigs (e.g. plasmid) into another file (this is done automatically)
#these contigs are added later again to the assembly
###-------------------- 
for plasmid_contigs in $(grep ">" ${ref_fasta}| sed 's/>//g' |grep -v "${bac_contig_nam}")
do
echo -e ${plasmid_contigs}"is a small remaining contig which is not start aligned"
samtools faidx ${ref_fasta} ${plasmid_contigs} >> ${out_dir}/without_bacteria_contig.fasta

done

###===========================
##start align
echo -e "-------2. Second step: start algin"
###===========================      

###--------------------
##genome length
###--------------------   
complete_Genome_length=$(cat ${out_dir}/only_bacteria_contig.fasta | awk '$0 ~ ">" {if (NR > 1) c=0 } $0 !~ ">" {c+=length($0);} END { print c }')
echo -e "calculated genome length of :" ${complete_Genome_length} "bp"

###--------------------
##extract fragement one
###-------------------- 
echo -e "Fragment one is located from :" ${bac_contig_nam}"\t"${Startalign_location}"\t"${complete_Genome_length}

echo -e ${bac_contig_nam}"\t"${Startalign_location}"\t"${complete_Genome_length} > ${out_dir}/tmp/fragment_1.bed
bedtools getfasta -fi ${out_dir}/only_bacteria_contig.fasta -bed ${out_dir}/tmp/fragment_1.bed > ${out_dir}/tmp/fragment_1.fasta

###--------------------
##extract fragement two
###--------------------  
echo -e "Fragment two is located from :" ${bac_contig_nam}"\t1\t"${Startalign_location}

echo -e ${bac_contig_nam}"\t1\t"${Startalign_location} > ${out_dir}/tmp/fragment_2.bed
bedtools getfasta -fi ${out_dir}/only_bacteria_contig.fasta -bed ${out_dir}/tmp/fragment_2.bed > ${out_dir}/tmp/fragment_2.fasta

###--------------------
##merge and thereby startalign fragments
###--------------------  
echo -e "merge fragments"

echo ">"${bac_contig_nam} > ${out_dir}/only_bacteria_contig_startaligned.fasta

grep -v ">" ${out_dir}/tmp/fragment_1.fasta >> ${out_dir}/only_bacteria_contig_startaligned.fasta
grep -v ">" ${out_dir}/tmp/fragment_2.fasta >> ${out_dir}/only_bacteria_contig_startaligned.fasta

###===========================
##reverse complement 
echo -e "-------3. third step: reverse complement"
###===========================   

revseq -sequence ${out_dir}/only_bacteria_contig_startaligned.fasta -outseq ${out_dir}/${sample}_contig_startaligned_oriented.fasta -notag 

###--------------------
##adding back the plasmid contigs
###--------------------
cat ${out_dir}/without_bacteria_contig.fasta >> ${out_dir}/${sample}_contig_startaligned_oriented.fasta

###===========================
##polish
echo -e "-------4. fourth step: polish genome"
###===========================   


###--------------------
##bowtie2 mapping bam file
###--------------------
echo -e "bowtie2 mapping"

mkdir -p ${out_dir}/startAligned_pillon_bowtie2_polishing/{bowtie2,assembly}

bowtie2-build --quiet ${out_dir}/${sample}_contig_startaligned_oriented.fasta ${out_dir}/${sample}_contig_startaligned_oriented.fasta

bowtie2 -1 ${raw_illumina}/${sample}_*_R1*_paired.fastq.gz \
        -2 ${raw_illumina}/${sample}_*_R2*_paired.fastq.gz \
        -x  ${out_dir}/${sample}_contig_startaligned_oriented.fasta \
        --threads ${threads} \
        --local --very-sensitive-local | samtools sort -O BAM -o ${out_dir}/startAligned_pillon_bowtie2_polishing/bowtie2/startAligned_pillon_Illumina2FinalAssembly_sorted.bam - 
     
###--------------------
##index bam file
###--------------------     
samtools index ${out_dir}/startAligned_pillon_bowtie2_polishing/bowtie2/startAligned_pillon_Illumina2FinalAssembly_sorted.bam

###--------------------
##pilon polishing
###--------------------
echo -e "pilon polishing"

pilon --threads ${threads} --genome ${out_dir}/${sample}_contig_startaligned_oriented.fasta \
  --frags $out_dir/startAligned_pillon_bowtie2_polishing/bowtie2/startAligned_pillon_Illumina2FinalAssembly_sorted.bam \
  --output Final_${sample}_startaligned \
  --outdir ${out_dir}/ \
  --changes

####--------------------------------------
##End of script
####--------------------------------------

date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
