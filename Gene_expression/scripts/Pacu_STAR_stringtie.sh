#!/bin/bash
#SBATCH --job-name=Pacu-stringtie
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Pacu-stringite.out  # %j = job ID
#SBATCH -e slurm-Pacu-stringtie.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pacu


# load modules needed
module load uri/main HISAT2/2.2.1-gompi-2021b #Alignment to reference genome: HISAT2
module load samtools/1.19.2 #Preparation of alignment for assembly: SAMtools
#load packages
module load uri/main StringTie/2.2.1-GCC-11.2.0

#Specify working directory
W="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pacu"
mkdir -p Stringtie

  #StringTie reference-guided assembly
  #These BAM files contain both forward and reverse reads
  array1=($(ls $W/*.bam))

  for i in ${array1[@]}; do
        stringtie -p 8 --fr -G /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.genes.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
        mv ${i}.gtf /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Pacu/Stringtie
        echo "StringTie-assembly-to-ref ${i}" $(date)
  done