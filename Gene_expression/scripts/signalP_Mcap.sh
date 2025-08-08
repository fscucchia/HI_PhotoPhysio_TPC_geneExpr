#!/bin/bash
#SBATCH --job-name=SignalP_Mcap
#SBATCH --nodes=1 --cpus-per-task=15
#SBATCH --mem=200G  # Requested Memory
#SBATCH -t 24:00:00
#SBATCH -o slurm-SignalP_Mcap.out  # %j = job ID
#SBATCH -e slurm-SignalP_Mcap.err  # %j = job ID
#SBATCH --mail-type=END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=federica.scucchia@uri.edu #your email to send notifications
#SBATCH -D /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files


# load modules needed
module load conda/latest
conda activate hisat

# Set up paths
SIGNALP_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/signalp-5.0b"
SPLIT_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files"

# Fix SignalP directory structure if needed
if [ ! -f "${SIGNALP_DIR}/bin/bin/signalp" ]; then
    echo "Creating expected SignalP directory structure..."
    mkdir -p ${SIGNALP_DIR}/bin/bin
    cp ${SIGNALP_DIR}/bin/signalp ${SIGNALP_DIR}/bin/bin/signalp
    chmod +x ${SIGNALP_DIR}/bin/bin/signalp
fi

# Test SignalP installation
echo "Testing SignalP installation..."
cd ${SIGNALP_DIR}
./bin/signalp -version
echo "SignalP version check complete, proceeding with analysis..."

# Process each AA split file with SignalP
for aa_file in ${SPLIT_DIR}/isoformSwitchAnalyzeR_isoform_AA_*.fasta; do
    if [ -f "$aa_file" ]; then
        file_num=$(basename "$aa_file" | sed 's/isoformSwitchAnalyzeR_isoform_AA_//' | sed 's/.fasta//')
        
        echo "Processing SignalP for file: $aa_file (part $file_num)"
        
        # Run SignalP from its own directory
        cd ${SIGNALP_DIR}
        
        # SignalP 5.0 with -stdout option to get summary output for IsoformSwitchAnalyzeR
        ./bin/signalp \
            -fasta "$aa_file" \
            -org euk \
            -format short \
            -prefix "SignalP_subset_${file_num}" \
            -batch 1000 \
            -stdout > "${SPLIT_DIR}/SignalP_subset_${file_num}.txt"
        
        # Check if SignalP succeeded and created output
        if [ $? -eq 0 ] && [ -f "${SPLIT_DIR}/SignalP_subset_${file_num}.txt" ]; then
            echo "SignalP completed for part $file_num - output saved to stdout"
        else
            echo "Warning: SignalP failed for part $file_num"
            
            # Check for alternative output files that SignalP 5.0 creates
            if [ -f "SignalP_subset_${file_num}_summary.signalp5" ]; then
                mv "SignalP_subset_${file_num}_summary.signalp5" "${SPLIT_DIR}/SignalP_subset_${file_num}.txt"
                echo "Found and moved summary file for part $file_num"
            else
                echo "No output files found for part $file_num"
                ls -la SignalP_subset_${file_num}* 2>/dev/null || echo "No files with that prefix"
            fi
        fi
        
        # Clean up temporary files
        rm -f SignalP_subset_${file_num}*.gff3 2>/dev/null
        rm -f SignalP_subset_${file_num}*.signalp5 2>/dev/null
        
        # Return to working directory
        cd ${SPLIT_DIR}
    fi
done

# Combine all SignalP summary results
echo "Combining SignalP results..."

# Find first valid output file for header
first_file=""
for file in ${SPLIT_DIR}/SignalP_subset_*.txt; do
    if [ -f "$file" ] && [ -s "$file" ]; then  # Check file exists and is not empty
        first_file="$file"
        break
    fi
done

if [ -n "$first_file" ]; then
    # Create combined file with header from first file
    head -1 "$first_file" > "${SPLIT_DIR}/signalp_summary_combined_results.txt"
    
    # Combine all data (skip headers)
    for file in ${SPLIT_DIR}/SignalP_subset_*.txt; do
        if [ -f "$file" ] && [ -s "$file" ]; then
            tail -n +2 "$file" >> "${SPLIT_DIR}/signalp_summary_combined_results.txt"
        fi
    done
    
    echo "SignalP summary analysis complete"
    echo "Combined results in: signalp_summary_combined_results.txt"
    
    # Show summary
    total_lines=$(wc -l < "${SPLIT_DIR}/signalp_summary_combined_results.txt")
    echo "Total sequences processed: $((total_lines - 1))"
    
    # Show sample of results
    echo "Sample results:"
    head -3 "${SPLIT_DIR}/signalp_summary_combined_results.txt"
    
else
    echo "Error: No valid SignalP output files found!"
    echo "Available files:"
    ls -la ${SPLIT_DIR}/SignalP_subset_* 2>/dev/null || echo "No SignalP output files"
fi