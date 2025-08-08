#!/bin/bash
#SBATCH --job-name=Deeploc2_Mcap
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-Deeploc2_Mcap.out  # %j = job ID
#SBATCH -e slurm-Deeploc2_Mcap.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files


# load modules needed
module load conda/latest
conda activate hisat

# Add local bin to PATH (where pip --user installs executables)
export PATH="/home/federica_scucchia_uri_edu/.local/bin:$PATH"

# Set up paths
SPLIT_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files"
OUTPUT_DIR="${SPLIT_DIR}/deeploc2_results"

# Process each AA split file
for aa_file in ${SPLIT_DIR}/isoformSwitchAnalyzeR_isoform_AA_*.fasta; do
    if [ -f "$aa_file" ]; then
        file_num=$(basename "$aa_file" | sed 's/isoformSwitchAnalyzeR_isoform_AA_//' | sed 's/.fasta//')
        
        echo "Processing DeepLoc2 for file: $aa_file (part $file_num)"
        
        # Run DeepLoc2 with Fast model
        deeploc2 \
          -f "$aa_file" \
          -o "${OUTPUT_DIR}/deeploc2_results_part_${file_num}" \
          -m Fast
        
        echo "Completed processing part $file_num"
    fi
done

# Combine all results
echo "Combining DeepLoc2 results..."
cat ${OUTPUT_DIR}/deeploc2_results_part_*/Results.tsv > ${OUTPUT_DIR}/coral_deeploc2_combined_results.tsv
  
# Count results
echo "Number of result files created:"
ls -1 ${OUTPUT_DIR}/deeploc2_results_part_* | wc -l

echo "Lines in combined results file:"
wc -l ${OUTPUT_DIR}/coral_deeploc2_combined_results.tsv