#!/bin/bash
#SBATCH --job-name=GeneLength
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-GeneLength.out  # %j = job ID
#SBATCH -e slurm-GeneLength.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom

#load modules
echo "Loading programs" $(date)
module load bedtools2/2.31.1

#This removes trailing tabs and spaces from each line before further processing and excludes comments and keeps only lines where the feature is "exon".
#Extract gene ID and exon coordinates
#Extracts gene_id from the attributes column
#Converts GTF coordinates (1-based) to BED format (0-based start, 1-based end)
#Sort for bedtools
#Sorts by chromosome and start position — required for bedtools merge.
#Merge overlapping exons per gene
#Merges overlapping or adjacent exon intervals per gene
#Ensures that overlapping bases are not double-counted
#Output: gene_lengths.tsv, gene_id and length of merged exon span (end - start)

grep -v "^#" stringtie_merged_noIso.gtf | \
sed 's/[ \t]*$//' | \
awk '$3=="exon"' | \
awk '{
    match($0, /gene_id "([^"]+)"/, arr);
    if (arr[1] != "") {
        print $1"\t"($4-1)"\t"$5"\t"arr[1]
    }
}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4 -o distinct -d 0 -i - | \
awk '{print $4"\t"$3 - $2}' > gene_lengths_Pcom.tsv
