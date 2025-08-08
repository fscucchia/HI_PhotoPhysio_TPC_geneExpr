#!/bin/bash
#SBATCH --job-name=Mcap-prepDe
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Mcap-prepDe_IsoformSwitch.out  # %j = job ID
#SBATCH -e slurm-Mcap-prepDe_IsoformSwitch.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Mcap


# If you want StringTie to discover new isoforms (novel transcripts), do not use the -e option. The -e option restricts StringTie to only quantifying transcripts present in the reference annotation, so it will not report novel isoforms.
# However, the prepDE.py3 script (used for generating count matrices for DE analysis) expects GTF files generated with -e. If you want to include novel isoforms in your downstream analysis, you need to:

# Run StringTie without -e to assemble and discover new isoforms.
# Merge all sample GTFs using stringtie --merge to create a comprehensive annotation (including novel isoforms).
# Re-run StringTie with -e and -G merged.gtf for each sample, using the merged annotation as reference. This will quantify both known and novel isoforms across all samples.
# Run prepDE.py3 on these new GTFs.

# load modules needed
#module load python/3.11.7
module load python/3.12.3 
#load packages
module load uri/main StringTie/2.2.1-GCC-11.2.0

# # run prepDe script
# python /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/scripts/prepDE.py3 -g ./gene_count_matrix_noIso_Mcap.csv -i ./sample_list.txt
# #IsoformSwitch
python /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/scripts/prepDE.py3 -g ./gene_count_matrix_IsoformSwitch_Mcap.csv -i ./sample_list_IsoformSwitch.txt
