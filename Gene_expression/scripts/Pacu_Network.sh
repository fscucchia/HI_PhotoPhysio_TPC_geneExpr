#!/bin/bash
#SBATCH --job-name=Pacu-networkBlue
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Pacu-networkBlue.out  # %j = job ID
#SBATCH -e slurm-Pacu-networkBlue.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu


# load modules needed
module load r/4.4.0

#Run R script
Rscript Network_Pacu.r > network_density_results_blue35_Pacu.txt