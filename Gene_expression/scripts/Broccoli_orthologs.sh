#!/bin/bash
#SBATCH --job-name=Mcap-broccoli
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=500G  # Requested Memory
#SBATCH -t 48:00:00
#SBATCH -o slurm-Mcap-broccoli.out  # %j = job ID
#SBATCH -e slurm-Mcap-broccoli.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Broccoli_orthologs/Broccoli-master

conda activate env-broccoli
module load diamond/2.1.10
python3 broccoli.py -dir /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Broccoli_orthologs