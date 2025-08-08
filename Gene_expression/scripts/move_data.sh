#!/bin/bash
#SBATCH --job-name=move_cleaned_seq_data
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 10:00:00
#SBATCH -o slurm-move_cleaned_seq_data.out  # %j = job ID
#SBATCH -e slurm-move_cleaned_seq_data.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/output

mv Cleaned_reads /scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output

echo "cleaned seq data moving complete." $(date)
