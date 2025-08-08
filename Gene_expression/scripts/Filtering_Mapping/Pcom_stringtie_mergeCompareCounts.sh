#!/bin/bash
#SBATCH --job-name=Pcom-mergeCompare
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Pcom-mergeCompare.out  # %j = job ID
#SBATCH -e slurm-Pcom-mergeCompare.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom


# load modules needed
#module load uri/main HISAT2/2.2.1-gompi-2021b #Alignment to reference genome: HISAT2
#module load samtools/1.19.2 #Preparation of alignment for assembly: SAMtools
#load packages
module load uri/main StringTie/2.2.1-GCC-11.2.0
module load uri/main GffCompare/0.12.6-GCC-11.2.0

## Merge the GTF files generated from the assembly to assess how well the predicted transcripts track to the
## reference annotation gff file. Stringtie --merge to create a comprehensive annotation (including novel isoforms).

#stringtie --merge -p 8 -G /work/pi_hputnam_uri_edu/snRNA_analysis/references/Porites_compressa_HIv_modified.gtf -o ./stringtie_merged.gtf list_to_merge.txt
#stringtie --merge -p 8 -G /work/pi_hputnam_uri_edu/snRNA_analysis/references/Porites_compressa_HIv_modified.gtf list_to_merge_noIso.txt -o ./stringtie_merged_noIso.gtf 

## Use the program gffcompare to compare the merged GTF files to the reference genome

#gffcompare -r /work/pi_hputnam_uri_edu/snRNA_analysis/references/Porites_compressa_HIv_modified.gtf stringtie_merged.gtf -o ./compared 
gffcompare -r /work/pi_hputnam_uri_edu/snRNA_analysis/references/Porites_compressa_HIv_modified.gtf stringtie_merged_noIso.gtf -o ./compared_noIso 
