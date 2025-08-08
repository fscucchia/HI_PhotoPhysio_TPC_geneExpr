#!/bin/bash
#SBATCH --job-name=DeepTMHMM_Mcap
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-DeepTMHMM_Mcap.out  # %j = job ID
#SBATCH -e slurm-DeepTMHMM_Mcap.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files


# load modules needed
module load conda/latest
conda activate hisat

# Set up paths
SPLIT_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files"

# Process each AA file with DeepTMHMM
for aa_file in ${SPLIT_DIR}/isoformSwitchAnalyzeR_isoform_AA_*.fasta; do
    if [ -f "$aa_file" ]; then
        file_num=$(basename "$aa_file" | sed 's/isoformSwitchAnalyzeR_isoform_AA_//' | sed 's/.fasta//')
        
        echo "Processing DeepTMHMM for file: $aa_file (part $file_num)"
        
        # Create output directory for this part
        output_dir="${SPLIT_DIR}/deeptmhmm_results_part_${file_num}"
        mkdir -p "$output_dir"
        
        # Change to output directory and run DeepTMHMM with correct syntax
        cd "$output_dir"
        
        # Correct biolib command syntax with --local flag
        biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta "$aa_file"
        
        # The results will be in the current directory (output_dir)
        echo "DeepTMHMM results for part $file_num saved in $output_dir"
    fi
done

# Combine DeepTMHMM results
echo "Combining DeepTMHMM results..."

# Find all result files (DeepTMHMM typically creates predicted_topologies.3line files)
find ${SPLIT_DIR}/deeptmhmm_results_part_* -name "*.3line" -exec cat {} \; > ${SPLIT_DIR}/deeptmhmm_combined_results.3line

# Also combine the detailed output if available
find ${SPLIT_DIR}/deeptmhmm_results_part_* -name "*.gff3" -exec cat {} \; > ${SPLIT_DIR}/deeptmhmm_combined_results.gff3

echo "DeepTMHMM processing complete"
