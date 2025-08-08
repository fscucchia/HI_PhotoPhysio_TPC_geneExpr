#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=36:00:00  # Job time limit
#SBATCH -o slurm-fastqc_raw.out  # %j = job ID
#SBATCH -e slurm-fastqc_raw.err  # %j = job ID
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /project/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/raw_data_genohub

#load modules 
module load uri/main
module load fastqc/0.12.1
module load MultiQC/1.12-foss-2021b

#run fastqc on raw data
fastqc /project/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/raw_data_genohub/*.fastq.gz -o /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/QC_rawData

#generate multiqc report
multiqc /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/QC_rawData --filename multiqc_report_raw.html 

echo "Initial QC of raw seq data complete." $(date)
