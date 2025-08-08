#!/bin/bash
#SBATCH --job-name=Pfam_Pcom
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Pfam_Pcom.out  # %j = job ID
#SBATCH -e slurm-Pfam_Pcom.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files


# load modules needed
module load conda/latest
conda activate hisat

pfam_scan.pl \
  -fasta /scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/isoformSwitchAnalyzeR_isoform_AA.fasta \
  -dir /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/PfamScan \
  -outfile Pcom_pfam_results_complete.txt

