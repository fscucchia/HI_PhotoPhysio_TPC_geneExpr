#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=36:00:00  # Job time limit
#SBATCH -o slurm-fastqc_raw.out  # %j = job ID
#SBATCH -e slurm-fastqc_raw.err  # %j = job ID
#SBATCH -D /project/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/raw_data_genohub


# Generate new MD5 checksums for all .fastq.gz files
#for file in /project/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/raw_data_genohub/*.fastq.gz; do
#    md5sum "$file" > "${file}.new.md5"
#done

# Compare the new MD5 checksums with the original ones
for file in /project/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/raw_data_genohub/*.fastq.gz; do
    original_md5="${file}.md5"
    new_md5="${file}.new.md5"
    if [ -f "$original_md5" ]; then
        diff "$original_md5" "$new_md5" > /dev/null
        if [ $? -eq 0 ]; then
            echo "MD5 checksum matches for $file" >> md5_verification_results.txt
        else
            echo "MD5 checksum mismatch for $file" >> md5_verification_results.txt
        fi
    else
        echo "Original MD5 file not found for $file" >> md5_verification_results.txt
    fi
done

echo "MD5 checksum verification complete." $(date) >> md5_verification_results.txt