#!/bin/bash
#SBATCH --job-name=Mcap-blastP
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=500G  # Requested Memory
#SBATCH -t 48:00:00
#SBATCH -o slurm-Mcap-blastP.out  # %j = job ID
#SBATCH -e slurm-Mcap-blastP.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Mcap


# load modules needed
#module load diamond/2.1.10
module load blast-plus/2.14.1
#module load uri/main all/InterProScan/5.60-92.0-foss-2021b

#diamond blastp -d /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pcom/nr.dmnd -q /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes.pep.faa -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Mcap/Mcap.diamondBlastpNCBInr -f 6 -b 4 --more-sensitive -e 0.00001 -k1       
#Use -f 6 to get the diamond results in tab format

# Run blastp


# Run InterProScan
#interproscan.sh -f TSV -i /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes.pep.faa -b mcap.interpro -iprlookup -goterms -pa 
