#!/bin/bash
#SBATCH --job-name=Pacu-InterproScan2
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=500G  # Requested Memory
#SBATCH -t 48:00:00
#SBATCH -o slurm-Pacu-InterproScan2.out  # %j = job ID
#SBATCH -e slurm-Pacu-InterproScan2.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pacu


# load modules needed
#module load diamond/2.1.10
#module load blast-plus/2.14.1
module load uri/main all/InterProScan/5.60-92.0-foss-2021b

#diamond blastp -d /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pcom/nr.dmnd -q /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.genes.pep.faa -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pacu/Pacu.diamondBlastpNCBInr_moreSensitive -f 6 -b 2 --more-sensitive --threads 8 -e 0.00001 -k1       
#Use -f 6 to get the diamond results in tab format

# Run InterProScan
interproscan.sh -f TSV -i /work/pi_hputnam_uri_edu/HI_Genomes/PacutaV2/Pocillopora_acuta_HIv2.genes.pep.faa -b pacu.interpro -iprlookup -goterms 
