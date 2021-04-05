#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name Prokka_startAligned
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 6G
#SBATCH --time 01:00:00
#SBATCH --output /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/prokka/annotation/002_%j.out
#SBATCH --error /users/dwijaya/scratch/dwijaya/Oxford_Nanopore/prokka/annotation/002_%j.err

####--------------------------------------
##preparation
##set you bash variables in order to quickly call them in the script
####--------------------------------------

user=dwijaya
sample=002
threads=8
home_dir=/users/${user}/scratch/${user}
out_dir=${home_dir}/Oxford_Nanopore/prokka/annotation/${sample}
ref_fasta=${home_dir}/Oxford_Nanopore/prokka/start_align/${sample}/Final_${sample}_startaligned.fasta

####---------------------
##modules
####--------------------------------------

module load HPC/Software 
module load UHTS/Analysis/prokka/1.13
module load UHTS/Analysis/barrnap/0.8

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: PROKKA annotation"

#done
echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${sample}


###===========================
##PROKKA
echo -e "-------1. PROKKA annotation"
###===========================   

prokka --compliant  --outdir  ${out_dir}/ --locustag ${sample} --proteins proteins.faa --evalue 0.001 --addgenes ${ref_fasta}

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"
