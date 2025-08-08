#!/bin/bash
#SBATCH --job-name=trimming
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=20:00:00  # Job time limit
#SBATCH -o slurm-trimming.out  # %j = job ID
#SBATCH -e slurm-trimming.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/Cleaned_reads

#load modules 
module load uri/main
module load fastqc/0.12.1
module load MultiQC/1.12-foss-2021b
module load fastp/0.23.2-GCC-11.2.0

# Make an array of sequences to trim in raw data directory 

cd /project/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/raw_data_genohub

array1=($(ls *R1_001.fastq.gz))

echo "Read trimming of adapters started." $(date)

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/Cleaned_reads/trim.${i} \
        --out2 /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/Cleaned_reads/trim.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 

done

echo "Read trimming of adapters completed." $(date)
