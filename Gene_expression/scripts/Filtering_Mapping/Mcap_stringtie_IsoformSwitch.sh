#!/bin/bash
#SBATCH --job-name=Mcap-stringtie_IsoformSwitch
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Mcap-stringtie_IsoformSwitch.out  # %j = job ID
#SBATCH -e slurm-Mcap-stringtie_IsoformSwitch.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap


# load modules needed
module load uri/main HISAT2/2.2.1-gompi-2021b #Alignment to reference genome: HISAT2
module load samtools/1.19.2 #Preparation of alignment for assembly: SAMtools
#load packages
module load uri/main StringTie/2.2.1-GCC-11.2.0


# #Specify working directory
# W="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Mcap"
# mkdir -p Stringtie

#   #StringTie reference-guided assembly
#   #These BAM files contain both forward and reverse reads
#   # Run StringTie without -e to assemble and discover new isoforms
#   array1=($(ls $W/*.bam))

#   for i in ${array1[@]}; do
#         stringtie -p 8 --rf -G /work/pi_hputnam_uri_edu/snRNA_analysis/references/Montipora_capitata_HIv3.gtf -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
#         mv ${i}.gtf /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_Mcap/Stringtie
#         echo "StringTie-assembly-to-ref ${i}" $(date)
#   done


#### 2 pass mode STAR

# #Specify working directory
# W="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap"
# mkdir -p Stringtie_IsoformSwitch

#   #StringTie reference-guided assembly
#   #These BAM files contain both forward and reverse reads
#   array1=($(ls $W/*.bam))

#   for i in ${array1[@]}; do
#         stringtie -p 8 --rf -G /work/pi_hputnam_uri_edu/snRNA_analysis/references/Montipora_capitata_HIv3.gtf -B -A ${i}.IsoformSwitch_gene_abund.tab -o ${i}_IsoformSwitch.gtf ${i}
#         mv ${i}_IsoformSwitch.gtf /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap/Stringtie_IsoformSwitch
#         echo "StringTie-assembly-to-ref IsoformSwitch ${i}" $(date)
#   done

## Re-run StringTie with -e and -G merged.gtf for each sample, using the merged annotation as reference. This will quantify both known and novel isoforms across all samples.
## After this re-run with the -e option, run directly prepDE.py3 on these new GTFs (prepDe needs files generated with -e option). 
## Do not run again --merge and gffcompare.

# #Specify working directory
W="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap"

  #StringTie reference-guided assembly
  #These BAM files contain both forward and reverse reads
  array1=($(ls $W/*.bam))

  for i in ${array1[@]}; do
        stringtie -p 8 --rf -e -G stringtie_merged_IsoformSwitch.gtf -B -A ${i}.gene_abund_IsoformSwitch_final.tab -o ${i}_IsoformSwitch_final.gtf ${i}
        mv ${i}_IsoformSwitch_final.gtf /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap/Stringtie_IsoformSwitch_final
        echo "StringTie-assembly-to-ref _IsoformSwitch_final ${i}" $(date)
  done