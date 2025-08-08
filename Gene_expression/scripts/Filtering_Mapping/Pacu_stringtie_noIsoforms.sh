#!/bin/bash
#SBATCH --job-name=Pacu-stringtie
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Pacu-stringtie_noIso.out  # %j = job ID
#SBATCH -e slurm-Pacu-stringtie_noIso.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pacu


# load modules needed
module load uri/main HISAT2/2.2.1-gompi-2021b #Alignment to reference genome: HISAT2
module load samtools/1.19.2 #Preparation of alignment for assembly: SAMtools
#load packages
module load uri/main StringTie/2.2.1-GCC-11.2.0

#Specify working directory
W="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pacu"
mkdir -p Stringtie_noIsoforms

  #StringTie reference-guided assembly
  #These BAM files contain both forward and reverse reads
  array1=($(ls $W/*.bam))

  for i in ${array1[@]}; do
        stringtie -p 8 --rf -e -G /work/pi_hputnam_uri_edu/snRNA_analysis/references/Pocillopora_acuta_HIv2_modified.gtf -A ${i}.gene_abund_noIsoforms.tab -o ${i}_noIsoforms.gtf ${i}
        mv ${i}_noIsoforms.gtf /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pacu/Stringtie_noIsoforms
        mv ${i}.gene_abund_noIsoforms.tab /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pacu/Stringtie_noIsoforms
        echo "StringTie-assembly-to-ref ${i}" $(date)
  done

