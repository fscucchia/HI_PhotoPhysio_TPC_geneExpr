#!/bin/bash
#SBATCH --job-name=Pacu-align
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-align.out  # %j = job ID
#SBATCH -e slurm-align.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Hisat_Pacu


# load modules needed
module load uri/main HISAT2/2.2.1-gompi-2021b #Alignment to reference genome: HISAT2
module load samtools/1.19.2 #Preparation of alignment for assembly: SAMtools
module load uri/main MultiQC/1.12-foss-2021b

# index the reference genome for Pocillopora acuta output index to working directory
#hisat2-build -f /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.assembly.fasta ./PacuHisat_ref # called the reference genome (scaffolds)
#echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
# The R1 in array1 is changed to R2 in the for loop. SAM files are of both forward and reverse reads
# Bam files are created and sorted, since Stringtie takes sorted file as input
# The sam file is removed at the end since it is not needed anymore
# The command --summary-file ${i}.txt reates a summary file per sample, which can be used by multiqc

# F="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output/Cleaned_reads"
# array1=($(ls $F/trim.Pacu*_R1_001.fastq.gz)) # Aligning paired end reads

# for i in ${array1[@]}; do
#     # Extract the base filename without the path and extension
#     base_name=$(basename ${i} _R1_001.fastq.gz)

#     # Construct output file paths
#     sam_output="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Hisat_Pacu/${base_name}.sam"
#     bam_output="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Hisat_Pacu/${base_name}.bam"
#     summary_output="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Hisat_Pacu/${base_name}.txt"

#     # Run HISAT2
#     hisat2 -p 8 --new-summary --rf --dta -q -x PacuHisat_ref \
#         -1 ${i} -2 $(echo ${i} | sed s/_R1/_R2/) \
#         -S ${sam_output} --summary-file ${summary_output}

#     # Sort BAM file
#     samtools sort -@ 8 -o ${bam_output} ${sam_output}

#     # Remove SAM file
#     rm ${sam_output}

#     # Log progress
#     echo "HISAT2 PE ${base_name}" $(date)
# done

H="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Hisat_Pacu"
mkdir -p $H/multiqc_hisat; multiqc $H