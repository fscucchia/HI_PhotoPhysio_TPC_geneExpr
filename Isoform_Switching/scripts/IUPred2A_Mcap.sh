#!/bin/bash
#SBATCH --job-name=IUPred2A_Mcap
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-IUPred2A_Mcap.out  # %j = job ID
#SBATCH -e slurm-IUPred2A_Mcap.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files


# load modules needed
module load conda/latest
conda activate hisat

# Set up working directory and paths
IUPRED_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/IUPred2A/iupred2a"
SPLIT_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files"

cd ${IUPRED_DIR}

# Process each AA split file with disorder prediction
for aa_file in ${SPLIT_DIR}/isoformSwitchAnalyzeR_isoform_AA_*.fasta; do
    if [ -f "$aa_file" ]; then
        file_num=$(basename "$aa_file" | sed 's/isoformSwitchAnalyzeR_isoform_AA_//' | sed 's/.fasta//')
        
        echo "Processing IUPred2A for file: $aa_file (part $file_num)"
        
        # Output file for this part
        output_file="${SPLIT_DIR}/iupred_long_results_part_${file_num}.txt"
        > "$output_file"  # Clear output file
        
        # Create temporary directory for this batch
        temp_dir="${SPLIT_DIR}/temp_iupred_${file_num}"
        mkdir -p "$temp_dir"
        
        # Split FASTA file into individual sequences
        awk '/^>/ {if(filename) close(filename); filename=sprintf("'$temp_dir'/seq_%04d.fasta", ++count)} {print > filename}' "$aa_file"
        
        # Process each individual sequence file
        for seq_file in ${temp_dir}/seq_*.fasta; do
            if [ -f "$seq_file" ]; then
                # Extract sequence ID from FASTA header
                seq_id=$(grep "^>" "$seq_file" | head -1 | sed 's/>//')
                
                # Run IUPred2A on this single sequence file
                python3 iupred2a.py "$seq_file" long | \
                awk -v id="$seq_id" 'NF>=3 && !/^#/ {print id "\t" $0}' >> "$output_file"
            fi
        done
        
        # Clean up temporary files
        rm -rf "$temp_dir"
        
        echo "Completed processing $aa_file ($(wc -l < "$output_file") lines)"
    fi
done