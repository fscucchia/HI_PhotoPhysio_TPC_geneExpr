#!/bin/bash
#SBATCH --job-name=Pcom-blastP
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -t 5-24:00:00
#SBATCH -q long
#SBATCH -o slurm-Pcom-blastP.out  # %j = job ID
#SBATCH -e slurm-Pcom-blastP.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pcom


# load modules needed
module load diamond/2.1.10
#module load blast-plus/2.14.1
#module load uri/main all/InterProScan/5.60-92.0-foss-2021b

# #Download the nr database from NCBI
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz #download nr database in fasta format, I downloaded version of 2024-02-07
# diamond makedb --in nr.gz -d nr  

#diamond blastp -d /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pcom/nr.dmnd -q /work/pi_hputnam_uri_edu/HI_Genomes/Pcompressa/Porites_compressa_HIv1.genes.pep.faa -o /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pcom/Pcom.diamondBlastpNCBInr_moreSensitive --more-sensitive -f 6 -b 2 -e 0.00001 -k1       
# #Use -f 6 to get the diamond results in tab format

## Run BLASTp
#Prepare the nr database for Blastp
gunzip nr.gz
#makeblastdb -in nr.fasta -dbtype prot

#blastp -query /work/pi_hputnam_uri_edu/HI_Genomes/Pcompressa/Porites_compressa_HIv1.genes.pep.faa -db nr.fasta -num_threads 15 -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out Past_annot_blastp



## Run InterProScan
#interproscan.sh -f TSV -i /work/pi_hputnam_uri_edu/HI_Genomes/Pcompressa/Porites_compressa_HIv1.genes.pep.faa -b pcom.interpro -iprlookup -goterms