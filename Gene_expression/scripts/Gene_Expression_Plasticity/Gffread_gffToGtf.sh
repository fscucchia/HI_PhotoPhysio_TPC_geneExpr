#!/bin/bash
#SBATCH --job-name=CladoDur-gffread
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-CladoDur-gffread.out  # %j = job ID
#SBATCH -e slurm-CladoDur-gffread.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Symbiont_Genomes_concat/Cladocopium_Durusdinium

#load modules
echo "Loading programs" $(date)
module load gffread/0.12.7

## use gffread to convert GFF3 to GTF format

gffread CladocDurusd_concat.gff -T -o CladocDurusd_concat.gtf
