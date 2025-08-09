
############## Analysis of changes in genome-wide patterns of alternative splicing and its functional consequences
########## using IsoformSwitchAnalyzeR

# Extract transcript lengths for normalization
library(GenomicFeatures)
library(txdbmaker)
library(rtracklayer)

setwd("/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom")

# Read GTF and fix strand information
gtf <- import("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom/stringtie_merged_IsoformSwitch.gtf")

# Fix strand information - convert "." to "+"
strand(gtf)[strand(gtf) == "."] <- "+"
strand(gtf)[strand(gtf) == "*"] <- "+"

# Export cleaned GTF to temporary file
temp_gtf <- tempfile(fileext = ".gtf")
export(gtf, temp_gtf, format = "gtf")

# Build TxDb from cleaned GTF
txdb <- makeTxDbFromGFF(temp_gtf, format="gtf")

# Get exons grouped by transcript
exons_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)

# Calculate actual transcript length as sum of exon lengths
transcript_lengths_proper <- sum(width(exons_by_tx))

# Convert to data frame
transcript_lengths <- data.frame(
    transcript_id = names(transcript_lengths_proper),
    length = as.integer(transcript_lengths_proper)
)

# Load count matrix and normalize
library(IsoformSwitchAnalyzeR)
library(edgeR)

transcriptCountMatrix <- read.table("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom/transcript_count_matrix.csv", 
                                   header = TRUE, row.names = 1, sep = ",")

# Match transcript lengths to the count matrix
matched_lengths <- transcript_lengths$length[match(rownames(transcriptCountMatrix), transcript_lengths$transcript_id)]

# Create DGEList and normalize for library size
dge <- DGEList(counts = transcriptCountMatrix)
dge <- calcNormFactors(dge)

# Calculate FPKM (accounts for both library size and transcript length)
transcriptFPKM <- rpkm(dge, gene.length = matched_lengths, normalized.lib.sizes = TRUE, log = FALSE)

## create a transcript FASTA file from StringTie GTF using gffread
#run in terminal
module load gffread/0.12.7
gffread -w transcripts.fa -g /work/pi_hputnam_uri_edu/HI_Genomes/Pcompressa/Porites_compressa_HIv1.assembly.fasta /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom/stringtie_merged_IsoformSwitch.gtf

# Create sampleTable for IsoformSwitchAnalyzeR
sampleTable <- data.frame(
    sampleID = c("B10", "F10", "G10", "H10",    # Control samples
                 "F7", "G7", "H7",        # 35°C samples
                 "B9", "F9", "G9", "H9"),       # 30°C samples
    condition = c("control", "control", "control", "control",
                  "t35", "t35", "t35",
                  "t30", "t30", "t30", "t30"),
    stringsAsFactors = FALSE
)

# Filter count matrix to keep only the samples needed
samples_to_keep <- sampleTable$sampleID
cat("Samples to keep:", paste(samples_to_keep, collapse = ", "), "\n")

# Remove Pcom_/Pcomp_ prefix from column names, keeping only what's after the underscore
colnames(transcriptCountMatrix) <- gsub("^Pcom[p]?_", "", colnames(transcriptCountMatrix))
colnames(transcriptFPKM) <- gsub("^Pcom[p]?_", "", colnames(transcriptFPKM))

# Check the result
cat("Updated column names:\n")
print(colnames(transcriptCountMatrix))

# Filter the count matrix and FPKM matrix to keep only needed samples
transcriptCountMatrix_filtered <- transcriptCountMatrix[, samples_to_keep]
transcriptFPKM_filtered <- transcriptFPKM[, samples_to_keep]

# Check that all sampleIDs match the filtered matrix columns
all_match <- all(sampleTable$sampleID %in% colnames(transcriptCountMatrix_filtered))
cat("Do all sampleTable IDs match filtered matrix?", all_match, "\n")

# Create switchList with normalized data
switchList <- importRdata(
    isoformCountMatrix = transcriptCountMatrix_filtered,  # Use filtered matrix
    isoformRepExpression = transcriptFPKM_filtered,       # Use filtered matrix
    designMatrix = sampleTable,
    isoformNtFasta = "transcripts.fa",
    isoformExonAnnoation = "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Mapping_twopassMode_Pcom/stringtie_merged_IsoformSwitch.gtf",
    showProgress = TRUE
)
# Step 10 of 10: Guestimating differential usage...
#     The GUESSTIMATED number of genes with differential isoform usage are:
#       comparison estimated_genes_with_dtu
# 1 control vs t30                538 - 897
# 2 control vs t35                345 - 576
# 3     t30 vs t35                447 - 745


# Filter the switchList to remove irrelevant genes/isoforms
switchListFiltered <- preFilter(
    switchAnalyzeRlist = switchList,
    geneExpressionCutoff = 0.1,        # Remove genes with low expression
    isoformExpressionCutoff = 0.1,   # Remove lowly expressed isoforms, Removes isoforms with expression < 0.1 FPKM
    removeSingleIsoformGenes = TRUE, # Remove genes with only one isoform (default)
    IFcutoff = 0.01                   # Remove isoforms contributing <1% to gene expression 
)
#The filtering removed 65042 ( 46.81% of ) transcripts. There is now 73921 isoforms left


# Testing for Isoform Switches via DEXSeq
# DEXSeq is the default and recommended test in IsoformSwitchAnalyzeR
# It controls FDR well and provides effect size corrected for confounding effects using limma

switchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = switchListFiltered,
    reduceToSwitchingGenes = TRUE,  # Reduce to genes with at least one differential isoform
    alpha = 0.05,  # FDR cutoff for significance
    dIFcutoff = 0.01,  # Minimum change in isoform fraction (1%)
    showProgress = TRUE
)
#Isoform switch analysis was performed for 53346 gene comparisons (100%).

# Summarize the isoform switching results
switchSummary <- extractSwitchSummary(switchListAnalyzed)
print("Isoform Switch Summary:")
print(switchSummary)
#       Comparison nrIsoforms nrSwitches nrGenes
# 1 control vs t30        627        821     513
# 2 control vs t35         85        155      76
# 3     t30 vs t35        466        659     372
# 4       Combined       1099       1578     871

# Get detailed switch analysis results
switchResults <- extractTopSwitches(
    switchAnalyzeRlist = switchListAnalyzed,
    filterForConsequences = FALSE,  # Don't filter for functional consequences yet
    n = Inf,  # Return all significant switches
    sortByQvals = TRUE  # Sort by q-values (FDR)
)

print(paste("Total number of significant isoform switches:", nrow(switchResults)))
#[1] "Total number of significant isoform switches: 961"

# Save the analyzed switchList for further analysis
#save(switchListAnalyzed, file = "switchListAnalyzed_Pcom.RData")

# Next Step: Extract Nucleotide and Amino Acid Sequences
# Since we imported from GTF/GFF file, ORF information should already be available
# We need to extract sequences for functional consequence analysis

# Check if ORF information is available
print(paste("Number of isoforms with ORF info:", sum(!is.na(switchListAnalyzed$orfAnalysis$orfTranscriptStart))))
#[1] "Number of isoforms with ORF info: 0" 
#ORF info not avaiable, need to extract sequences from GTF

### Analyzing Known and Novel Isoforms: Extract sequences for isoforms with ORF information

# First: Add CDS information from reference annotation (for known isoforms)
# You'll need the original reference GTF file that was used with StringTie
switchListAnalyzed <- addORFfromGTF(
    switchAnalyzeRlist = switchListAnalyzed,
    pathToGTF = "/work/pi_hputnam_uri_edu/snRNA_analysis/references/Porites_compressa_HIv_modified.gtf"
)
# Step 1 of 2: importing GTF (this may take a while)...
# Step 2 of 2: Adding ORF...
#     Added ORF info (incl info about isoforms annotated as not having an ORF) to 757 isoforms.
#         This correspond to 11.87% of isoforms in the switchAnalyzeRlist.

#Only 11.87% of isoforms got ORF info - This is normal for StringTie because most of the isoforms are novel transcripts discovered by StringTie
# 757 isoforms got ORF annotations from known genes
# The remaining isoforms are novel discoveries that need ORF prediction

## Second: predict ORFs for novel isoforms discovered by StringTie
# The analyzeORF() function will automatically analyze all isoforms that don't already have ORF information. Since the CDS information for 
# known isoforms is already added using addORFfromGTF(), the function will only work on the remaining novel isoforms that need ORF prediction.
switchListAnalyzed <- analyzeORF(
    switchAnalyzeRlist = switchListAnalyzed,
    showProgress = TRUE
    #method = "longest.AnnotatedWhenPossible"  # Default and recommended method
)
#6284 putative ORFs were identified, analyzed and added.

# Check ORF availability again
print("After ORF analysis:")
print(paste("Number of isoforms with ORF info:", sum(!is.na(switchListAnalyzed$orfAnalysis$orfTransciptStart))))
#[1] "Number of isoforms with ORF info: 6284"

# proceed with sequence extraction
switchListAnalyzed <- extractSequence(
    switchAnalyzeRlist = switchListAnalyzed,
    pathToOutput = "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom",
    writeToFile = TRUE
)

# Re-extract sequences with splitting for webserver limitations
switchListAnalyzed <- extractSequence(
    switchAnalyzeRlist = switchListAnalyzed,
    pathToOutput = "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files",
    writeToFile = TRUE,
    removeLongAAseq=TRUE,  # Remove long AA sequences 
    alsoSplitFastaFile = TRUE     # Split into multiple files for webserver upload
)

#Use R to split the NT file
library(Biostrings)

# Read the NT FASTA file
nt_seqs <- readDNAStringSet("/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files/isoformSwitchAnalyzeR_isoform_nt.fasta")

# Split into chunks of 13000 sequences
chunk_size <- 13000
n_chunks <- ceiling(length(nt_seqs) / chunk_size)

for(i in 1:n_chunks) {
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, length(nt_seqs))
  
  chunk_seqs <- nt_seqs[start_idx:end_idx]
  
  writeXStringSet(chunk_seqs, 
                  filepath = paste0("/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files/switchList_nt_part_", i, ".fasta"))
}

# Get the total number of genes in Pcom genome
# Option 1: From reference GTF file
ref_gtf <- import("/work/pi_hputnam_uri_edu/snRNA_analysis/references/Porites_compressa_HIv_modified.gtf")
total_genes_genome <- length(unique(ref_gtf$gene_id))

genes_with_switching <- 871  # From Combined nrGenes result

# Calculate global percentages for cross-species comparison
percentage_switching <- (genes_with_switching / total_genes_genome) * 100

print(paste("Genes with temperature-responsive isoform switching:", genes_with_switching))
#[1]  "Genes with temperature-responsive isoform switching: 871"
print(paste("Total genes in genome:", total_genes_genome))
#[1]  "Total genes in genome: 44130"
print(paste("Percentage of genome showing isoform switching:", round(percentage_switching, 2), "%"))
#[1] "Percentage of genome showing isoform switching: 1.97 %"

total_isoforms_analyzed <- nrow(switchListAnalyzed$isoformFeatures)
switching_isoforms <- 1099  # From the total switches

print(paste("Switching isoforms out of total analyzed:", 
           round((switching_isoforms/total_isoforms_analyzed)*100, 2), "%"))
#[1] "Switching isoforms out of total analyzed: 5.74%"



####### Running External Sequence Analysis Tools and Downloading Results

# Workflow:
# Start with CPC2 → Filter to coding isoforms
# Run Pfam → Identify functional domains
# Run SignalP → Find secreted/membrane proteins
# Run IUPred2A → Identify disordered regions

# 1. CPC2 - Coding Potential Prediction
# Determines which of the 19,168 ORFs are likely protein-coding vs non-coding
# Use: Nucleotide FASTA file (_nt.fasta)
# Webserver: http://cpc2.gao-lab.org/
# Critical for: Filtering out non-coding transcripts before functional analysis

# 2. Pfam - Protein Domain Prediction
# Identifies functional protein domains in the switching isoforms
# Use: Amino acid FASTA file (_AA.fasta)
# ran as batch job, see Pfam_Mcap.sh
# Critical for: Understanding functional consequences of isoform switches

conda install bioconda::pfam_scan

# Download pfam_scan.pl script
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz
tar -xzf PfamScan.tar.gz
cd PfamScan

# Install required Perl modules
conda install -c bioconda perl-moose perl-bioperl

# Download Pfam database (this is large - ~3GB compressed, ~9GB uncompressed)
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz

# Uncompress database files
gunzip Pfam-A.hmm.gz
gunzip Pfam-A.hmm.dat.gz  
gunzip active_site.dat.gz

# Index the HMM database (required step)
hmmpress Pfam-A.hmm

# Basic command structure - run as batch job, see Pfam_Mcap.sh
pfam_scan.pl \
  -fasta /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/isoformSwitchAnalyzeR_isoform_AA.fasta \
  -dir /scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/PfamScan \
  -outfile coral_pfam_results_complete.txt \
  -cpu 8


# 3. SignalP - Signal Peptide Prediction
# Prediction of Signal Peptides — a short N-terminal sequence of a peptide indicating where a protein should be membrane bound or secreted.
# Use: Amino acid FASTA file (_AA.fasta)
# run as batch job, see signalP_Mcap.sh

# 4. IUPred2A - Intrinsically Disordered Regions
# The parts of a protein which does not have a fixed three-dimensional structure (as opposite protein domains)
# Use: Amino acid FASTA file (_AA.fasta)
# Webserver: https://iupred2a.elte.hu/
#  - run as batch job, see IUPred2A_Mcap.sh


# 5. DeepLoc2 - Subcellular Localization
# Understanding where switching proteins function in the cell
#  - run as batch job, see Deeploc2_Mcap.sh

# Navigate to where you downloaded the DeepLoc package
cd /path/to/downloaded/deeploc2

# Install DeepLoc 2.1
pip install DeepLoc-2.1.0.tar.gz

# OR if you're in the deeploc2_package directory:
pip install . --user  ## within conda environment

# Add the local bin to PATH
export PATH="/home/federica_scucchia_uri_edu/.local/bin:$PATH"

# Test installation
deeploc2 -f test.fasta


# # # or 5. DeepTHMMM - Subcellular Localization
# # # Understanding where switching proteins function in the cell
# # # run as batch job, see DeepTHMMM_Mcap.sh

# pip3 install --upgrade pybiolib

# # Set up paths
# SPLIT_DIR="/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files"

# # Process each AA file with DeepTMHMM
# for aa_file in ${SPLIT_DIR}/isoformSwitchAnalyzeR_isoform_AA_*.fasta; do
#     if [ -f "$aa_file" ]; then
#         file_num=$(basename "$aa_file" | sed 's/isoformSwitchAnalyzeR_isoform_AA_//' | sed 's/.fasta//')
        
#         echo "Processing DeepTMHMM for file: $aa_file (part $file_num)"
        
#         # Create output directory for this part
#         output_dir="${SPLIT_DIR}/deeptmhmm_results_part_${file_num}"
#         mkdir -p "$output_dir"
        
#         # Change to output directory and run DeepTMHMM with correct syntax
#         cd "$output_dir"
        
#         # Correct biolib command syntax with --local flag
#         biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta "$aa_file"
        
#         # The results will be in the current directory (output_dir)
#         echo "DeepTMHMM results for part $file_num saved in $output_dir"
#     fi
# done

# # Combine DeepTMHMM results
# echo "Combining DeepTMHMM results..."

# # Find all result files (DeepTMHMM typically creates predicted_topologies.3line files)
# find ${SPLIT_DIR}/deeptmhmm_results_part_* -name "*.3line" -exec cat {} \; > ${SPLIT_DIR}/deeptmhmm_combined_results.3line

# # Also combine the detailed output if available
# find ${SPLIT_DIR}/deeptmhmm_results_part_* -name "*.gff3" -exec cat {} \; > ${SPLIT_DIR}/deeptmhmm_combined_results.gff3

# echo "DeepTMHMM processing complete"





##### Combine split results

# Read CPC2 results
cpc2_combined <- read.table("/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files/result_cpc2.txt", header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)

# Write the file
write.table(cpc2_combined, 
           "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files/cpc2_results.txt",
           sep = "\t", row.names = FALSE, quote = FALSE)




# # Combine SignalP results (all 32 files)
# signalp_files <- list.files(
#   "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/",
#   pattern = "SignalP_subset_.*\\.txt$",
#   full.names = TRUE
# )

# cat("Found", length(signalp_files), "SignalP files\n")
# #Found 32 SignalP files

# # Read and combine SignalP results
# read_signalp_file <- function(file_path) {
#   tryCatch({
#     # Read the file and handle the header properly
#     # Skip the first comment line, use the second line as header
#     df <- read.table(file_path, header = FALSE, sep = "\t", skip = 1, 
#                      stringsAsFactors = FALSE, comment.char = "#")
    
#     # Set proper column names based on SignalP format
#     # From the header line: "# ID\tPrediction\tSP(Sec/SPI)\tOTHER\tCS Position"
#     colnames(df) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
    
#     return(df)
#   }, error = function(e) {
#     cat("Error reading", basename(file_path), ":", e$message, "\n")
#     return(NULL)
#   })
# }

# # Read all SignalP files
# signalp_list <- lapply(signalp_files, read_signalp_file)

# # Remove NULL entries (failed reads)
# signalp_list <- signalp_list[!sapply(signalp_list, is.null)]

# if(length(signalp_list) > 0) {
#   # Combine all SignalP files
#   signalp_combined <- do.call(rbind, signalp_list)
  
#   # Write combined file
#   write.table(signalp_combined,
#              "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/signalp_combined_results.txt",
#              sep = "\t", row.names = FALSE, quote = FALSE)
  
#   cat("SignalP files combined successfully. Rows:", nrow(signalp_combined), "\n")
#   cat("Column names:", paste(colnames(signalp_combined), collapse = ", "), "\n")
  
#   # Show first few rows
#   cat("First few rows:\n")
#   print(head(signalp_combined, 3))
# } else {
#   cat("No SignalP files could be read successfully\n")
# }
# # SignalP files combined successfully. Rows: 15641 
# # Column names: ID, Prediction, SP_Sec_SPI, OTHER, CS_Position 
# # First few rows:
# #                                         ID Prediction SP_Sec_SPI    OTHER CS_Position
# # 1 Montipora_capitata_HIv3___RNAseq.11286_t      OTHER   0.000739 0.999261            
# # 2 Montipora_capitata_HIv3___RNAseq.12942_t      OTHER   0.000561 0.999439            
# # 3 Montipora_capitata_HIv3___RNAseq.13784_t      OTHER   0.001418 0.998582   

# # Convert "SP(Sec/SPI)" to "SP" to standardize signal peptide predictions
# signalp_combined$Prediction[signalp_combined$Prediction == "SP(Sec/SPI)"] <- "SP"

# # Total signal peptides should now be 436 + 983 = 1,419
# total_sp <- sum(signalp_combined$Prediction == "SP")
# cat("Total signal peptides after standardization:", total_sp, "\n")

# # Write the corrected SignalP file
# write.table(signalp_combined,
#            paste0(base_path, "signalp_corrected_results.txt"),
#            sep = "\t", row.names = FALSE, quote = FALSE)

# # # Load SignalP data
# base_path <- "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/"

# cat("Recreating SignalP combined file with correct format...\n")

# # Find all SignalP output files
# signalp_files <- list.files(base_path, 
#                            pattern = "SignalP_subset_.*\\.txt$", 
#                            full.names = TRUE)

# cat("Found", length(signalp_files), "SignalP files\n")

# # Create new combined file with EXACT format from working file
# combined_file <- paste0(base_path, "signalp_exact_format_combined.txt")

# # Write the EXACT header format from the working file
# cat("# SignalP-5.0\tOrganism: euk\tTimestamp: 20250725150847\n", file = combined_file)
# cat("# ID\tPrediction\tSP(Sec/SPI)\tOTHER\tCS Position\n", file = combined_file, append = TRUE)

# # Combine data from all files (skip headers in each file)
# total_sequences <- 0
# for(file in signalp_files) {
#     if(file.exists(file)) {
#         file_lines <- readLines(file)
        
#         # Find where data starts (after the header lines)
#         data_start <- which(!grepl("^#", file_lines))[1]
        
#         if(!is.na(data_start) && data_start <= length(file_lines)) {
#             data_lines <- file_lines[data_start:length(file_lines)]
            
#             # Remove empty lines
#             data_lines <- data_lines[data_lines != "" & !is.na(data_lines)]
            
#             # Append to combined file
#             if(length(data_lines) > 0) {
#                 cat(data_lines, file = combined_file, sep = "\n", append = TRUE)
#                 total_sequences <- total_sequences + length(data_lines)
#             }
#         }
#     }
# }

# cat("Created combined file with", total_sequences, "sequences\n")

# # Check the format of the new combined file
# cat("New combined file format:\n")
# print(readLines(combined_file, n = 10))

# # Test SignalP analysis with the corrected format
# cat("Testing SignalP with corrected format (Organism: euk)...\n")
# tryCatch({
#     switchListAnalyzed <- analyzeSignalP(
#         switchAnalyzeRlist = switchListAnalyzed,
#         pathToSignalPresultFile = combined_file
#     )
#     cat("SUCCESS: SignalP analysis worked with corrected format!\n")
    
#     # Check results
#     if("signalPeptideAnalysis" %in% names(switchListAnalyzed)) {
#         sp_count <- length(unique(switchListAnalyzed$signalPeptideAnalysis$isoform_id))
#         cat("Signal peptides added to", sp_count, "transcripts\n")
        
#         # Show sample results
#         cat("Sample SignalP results:\n")
#         print(head(switchListAnalyzed$signalPeptideAnalysis, 3))
#     }
    
# }, error = function(e) {
#     cat("SignalP analysis still failed:", e$message, "\n")
    
#     # Compare formats side by side
#     cat("\nFormat comparison:\n")
#     cat("Working single file header:\n")
#     working_header <- readLines(paste0(base_path, "SignalP_subset_subset_8_of_32.txt"), n = 2)
#     print(working_header)
    
#     cat("New combined file header:\n")
#     combined_header <- readLines(combined_file, n = 2)
#     print(combined_header)
    
#     cat("Are headers identical?", identical(working_header, combined_header), "\n")
# })


 # Fix IUPred2A analysis with ID mapping ----  DOESN'T WORK
# base_path <- "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/"

# # Read the IUPred2A combined results
# iupred_data <- read.table(paste0(base_path, "iupred_long_results_part_complete.txt"), 
#                          header = FALSE, sep = "\t", stringsAsFactors = FALSE,
#                          comment.char = "#", fill = TRUE)

# cat("Current command-line format:\n")
# cat("Columns:", ncol(iupred_data), "\n")
# print("Sample data:")
# print(head(iupred_data, 5))

# # Convert to web server format (remove sequence IDs, add proper headers)
# # Web server format should be: position, amino_acid, iupred_score, anchor_score

# # If the data has 4 columns (seq_id, pos, aa, score), convert to web format
# if(ncol(iupred_data) == 4) {
#     # Remove sequence ID column and convert to web server format
#     web_format_data <- iupred_data[, 2:4]  # Keep position, amino_acid, score
    
#     # Add a dummy ANCHOR score column (if not available)
#     web_format_data$anchor_score <- 0.5000  # Default anchor score
    
#     # Reorder columns to match web server: POS, AMINO ACID, IUPRED SCORE, ANCHOR SCORE
#     web_format_data <- web_format_data[, c(1, 2, 3, 4)]
    
# } else if(ncol(iupred_data) == 5) {
#     # If you have anchor scores, use them
#     web_format_data <- iupred_data[, 2:5]  # Remove seq_id, keep the rest
    
# } else {
#     cat("Unexpected number of columns:", ncol(iupred_data), "\n")
#     stop("Please check the IUPred2A data format")
# }

# # Create web server format file with proper headers
# web_server_file <- paste0(base_path, "iupred_webserver_format.txt")

# # Write the web server format file with headers (like the online tool)
# cat("# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding\n", 
#     file = web_server_file)
# cat("# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi\n", 
#     file = web_server_file, append = TRUE)
# cat("# Nucleic Acids Research 2018, Submitted\n", 
#     file = web_server_file, append = TRUE)
# cat("# IUPred2 type: long\n", 
#     file = web_server_file, append = TRUE)
# cat("# POS\tAMINO ACID\tIUPRED SCORE\tANCHOR SCORE\n", 
#     file = web_server_file, append = TRUE)

# # Append the data without headers
# write.table(web_format_data, 
#            file = web_server_file, 
#            sep = "\t", 
#            row.names = FALSE, 
#            col.names = FALSE, 
#            quote = FALSE, 
#            append = TRUE)

# cat("Created web server format file\n")
# cat("First few lines of web server format:\n")
# print(readLines(web_server_file, n = 10))

# # Test IUPred2A analysis with web server format
# cat("Testing IUPred2A with web server format...\n")
# tryCatch({
#     switchListAnalyzed <- analyzeIUPred2A(
#         switchAnalyzeRlist = switchListAnalyzed,
#         pathToIUPred2AresultFile = web_server_file,
#         showProgress = TRUE
#     )
#     cat("SUCCESS: IUPred2A analysis worked with web server format!\n")
    
#     # Check results
#     if("idrAnalysis" %in% names(switchListAnalyzed)) {
#         idr_count <- length(unique(switchListAnalyzed$idrAnalysis$isoform_id))
#         cat("IDR regions added to", idr_count, "transcripts\n")
        
#         # Show sample results
#         cat("Sample IDR analysis:\n")
#         print(head(switchListAnalyzed$idrAnalysis, 3))
#     }
    
# }, error = function(e) {
#     cat("Web server format still failed:", e$message, "\n")
    
#     # Try alternative: split by sequence and create separate files
#     cat("Trying to split into individual sequence files...\n")
    
#     # Group by sequence ID and create separate files for each
#     unique_seqs <- unique(iupred_data[,1])
#     cat("Found", length(unique_seqs), "unique sequences\n")
    
#     # Create a single combined web server file that processes all sequences
#     combined_web_file <- paste0(base_path, "iupred_all_sequences_web_format.txt")
    
#     # Write headers once
#     cat("# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding\n", 
#         file = combined_web_file)
#     cat("# Multiple sequences processed\n", 
#         file = combined_web_file, append = TRUE)
#     cat("# POS\tAMINO ACID\tIUPRED SCORE\tANCHOR SCORE\n", 
#         file = combined_web_file, append = TRUE)
    
#     # Add all data regardless of sequence ID
#     write.table(web_format_data, 
#                file = combined_web_file, 
#                sep = "\t", 
#                row.names = FALSE, 
#                col.names = FALSE, 
#                quote = FALSE, 
#                append = TRUE)
    
#     # Test again
#     tryCatch({
#         switchListAnalyzed <- analyzeIUPred2A(
#             switchAnalyzeRlist = switchListAnalyzed,
#             pathToIUPred2AresultFile = combined_web_file,
#             showProgress = TRUE
#         )
#         cat("SUCCESS: Combined web format worked!\n")
        
#     }, error = function(e2) {
#         cat("All IUPred2A attempts failed:", e2$message, "\n")
#         cat("Proceeding without IUPred2A analysis.\n")
#     })
# })




# # Combine DeepLoc2 results from subfolders ---  DOESN'T WORK
# deeploc_dirs <- list.dirs(
#   "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/deeploc2_results/",
#   recursive = FALSE
# )

# # Find all Results.tsv files in subdirectories
# deeploc_files <- list.files(
#   "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/deeploc2_results/",
#   pattern = "results_.*\\.csv$",  # Look for results_TIMESTAMP.csv files
#   recursive = TRUE,
#   full.names = TRUE
# )

# cat("Found", length(deeploc_files), "DeepLoc2 files\n")
# #Found 33 DeepLoc2 files

# if(length(deeploc_files) > 0) {
#   cat("DeepLoc2 files found:\n")
#   print(basename(deeploc_files)[1:min(5, length(deeploc_files))])  # Show first 5 files
  
#   # Check first file structure
#   cat("\nFirst DeepLoc2 file structure:\n")
#   first_file <- read.table(deeploc_files[1], header = TRUE, sep = ",", nrows = 3)  # CSV files use comma separator
#   print(head(first_file))
  
#   # Read and combine DeepLoc2 results (CSV format)
#   deeploc_combined <- do.call(rbind, lapply(deeploc_files, function(x) {
#     tryCatch({
#       read.table(x, header = TRUE, sep = ",", stringsAsFactors = FALSE)  # Use comma separator for CSV
#     }, error = function(e) {
#       cat("Error reading", basename(x), ":", e$message, "\n")
#       return(NULL)
#     })
#   }))
  

#   # Remove NULL entries if any failed
#   deeploc_combined <- deeploc_combined[!sapply(deeploc_combined, is.null)]
  
#   if(nrow(deeploc_combined) > 0) {
#     write.table(deeploc_combined,
#                "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/deeploc2_combined_results.tsv",
#                sep = "\t", row.names = FALSE, quote = FALSE)
    
#     cat("DeepLoc2 files combined successfully. Rows:", nrow(deeploc_combined), "\n")
#     cat("Column names:", paste(colnames(deeploc_combined), collapse = ", "), "\n")
    
#     # Show first few rows
#     cat("First few rows:\n")
#     print(head(deeploc_combined, 3))
#   } else {
#     cat("No DeepLoc2 data could be combined\n")
#   }
# } else {
#   cat("No DeepLoc2 files found with pattern 'results_*.csv'\n")
  
#   # Alternative search - list all files in deeploc2_results directory
#   cat("All files in deeploc2_results directory:\n")
#   all_files <- list.files(
#     "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/deeploc2_results/",
#     recursive = TRUE,
#     full.names = TRUE
#   )
#   print(basename(all_files))
# }


# # The command-line DeepLoc2 has 5 extra columns that the online tool doesn't have:
# # Membrane.types
# # Peripheral
# # Transmembrane
# # Lipid.anchor
# # Soluble
# # IsoformSwitchAnalyzeR probably expects the online tool format (13 columns) rather than the command-line format (18 columns).

# # Fix DeepLoc2 format by removing extra columns to match online tool format
# base_path <- "/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Alternative_splicing/Mcap/split_files/"

# # Define the columns that the online tool has (and IsoformSwitchAnalyzeR expects)
# expected_columns <- c(
#     "Protein_ID", "Localizations", "Signals", "Cytoplasm", "Nucleus", 
#     "Extracellular", "Cell.membrane", "Mitochondrion", "Plastid", 
#     "Endoplasmic.reticulum", "Lysosome.Vacuole", "Golgi.apparatus", "Peroxisome"
# )

# # Read the combined DeepLoc2 file
# deeploc_combined <- read.table(paste0(base_path, "deeploc2_combined_results.tsv"), 
#                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# cat("Original DeepLoc2 file:\n")
# cat("Columns:", ncol(deeploc_combined), "\n")
# cat("Rows:", nrow(deeploc_combined), "\n")

# # Remove the extra columns to match online format
# deeploc_fixed <- deeploc_combined[, expected_columns]

# cat("Fixed DeepLoc2 file (online format):\n")
# cat("Columns:", ncol(deeploc_fixed), "\n")
# cat("Column names:", paste(colnames(deeploc_fixed), collapse = ", "), "\n")
# #Column names: Protein_ID, Localizations, Signals, Cytoplasm, Nucleus, Extracellular, Cell.membrane, Mitochondrion, Plastid, Endoplasmic.reticulum, Lysosome.Vacuole, Golgi.apparatus, Peroxisome 

# # Apply ID mapping using the same approach as SignalP
# deeploc_mapped <- merge(deeploc_fixed, id_mapping, 
#                        by.x = "Protein_ID", by.y = "original_id", 
#                        all.x = FALSE)

# # Replace Protein_ID with internal IDs
# deeploc_mapped$Protein_ID <- deeploc_mapped$internal_id
# deeploc_mapped$internal_id <- NULL

# cat("DeepLoc2 entries after ID mapping:", nrow(deeploc_mapped), "\n")

# # Write the fixed and mapped DeepLoc2 file as CSV (online tool format)
# write.csv(deeploc_mapped,
#           paste0(base_path, "deeploc2_online_format.csv"),
#           row.names = FALSE)

# # Test DeepLoc2 analysis with online format
# cat("Testing DeepLoc2 with online tool format...\n")
# tryCatch({
#     switchListAnalyzed <- analyzeDeepLoc2(
#         switchAnalyzeRlist = switchListAnalyzed,
#         pathToDeepLoc2resultFile = paste0(base_path, "deeploc2_online_format.csv")
#     )
#     cat("SUCCESS: DeepLoc2 analysis works with online format!\n")
    
#     # Check results
#     if("subcellularOrganelleAnalysis" %in% names(switchListAnalyzed)) {
#         subcell_count <- length(unique(switchListAnalyzed$subcellularOrganelleAnalysis$isoform_id))
#         cat("Subcellular localization added to", subcell_count, "transcripts\n")
#     }
    
# }, error = function(e) {
#     cat("DeepLoc2 still failed:", e$message, "\n")
    
#     # Try TSV format as backup
#     write.table(deeploc_mapped,
#                paste0(base_path, "deeploc2_online_format.tsv"),
#                sep = "\t", row.names = FALSE, quote = FALSE)
    
#     tryCatch({
#         switchListAnalyzed <- analyzeDeepLoc2(
#             switchAnalyzeRlist = switchListAnalyzed,
#             pathToDeepLoc2resultFile = paste0(base_path, "deeploc2_online_format.tsv")
#         )
#         cat("SUCCESS: DeepLoc2 works with TSV format!\n")
        
#     }, error = function(e2) {
#         cat("Both CSV and TSV failed:", e2$message, "\n")
#     })
# })






### import all results into IsoformSwitchAnalyzeR

# Set base path for convenience
base_path <- "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/split_files/"

# Add CPC2 analysis (coding potential prediction)
switchListAnalyzed <- analyzeCPC2(
    switchAnalyzeRlist = switchListAnalyzed,
    pathToCPC2resultFile = paste0(base_path, "cpc2_results.txt"),
    removeNoncodinORFs = TRUE  # Remove ORFs from non-coding transcripts
    )
#Added coding potential to 5170 (81.06%) transcripts

# Add Pfam analysis (protein domains)
switchListAnalyzed <- analyzePFAM(
    switchAnalyzeRlist = switchListAnalyzed,
    pathToPFAMresultFile = paste0(base_path, "Pcom_pfam_results_complete.txt"),
    showProgress = TRUE
)
#Added domain information to 3752 (58.83%) transcripts

# # # Add SignalP analysis (signal peptides) ---  DOESN'T WORK AT ALL
# switchListAnalyzed <- analyzeSignalP(
#     switchAnalyzeRlist = switchListAnalyzed,
#     pathToSignalPresultFile = paste0(base_path, "SignalP_subset_subset_8_of_32.txt")
# )

# # Add IUPred2A analysis (intrinsically disordered regions) ---  DOESN'T WORK AT ALL
# switchListAnalyzed <- analyzeIUPred2A(
#     switchAnalyzeRlist = switchListAnalyzed,
#     pathToIUPred2AresultFile = paste0(base_path, "coral_iupred_combined_results.txt")
# )

# Add DeepLoc2 analysis (subcellular localization)  ---  DOESN'T WORK AT ALL
# switchListAnalyzed <- analyzeDeepLoc2(
#     switchAnalyzeRlist = switchListAnalyzed,
#     pathToDeepLoc2resultFile = paste0(base_path, "deeploc2_combined_results.tsv")
#     )



### Alternative Splicing Analysis

# Analyze alternative splicing patterns in the switching isoforms
switchListAnalyzed <- analyzeAlternativeSplicing(
    switchAnalyzeRlist = switchListAnalyzed,
    quiet = FALSE  # Set to FALSE to see progress
)

# Check what types of alternative splicing events were identified
if("AlternativeSplicingAnalysis" %in% names(switchListAnalyzed)) {
    # Overview of intron retention events
    cat("Intron Retention (IR) events:\n")
    ir_table <- table(switchListAnalyzed$AlternativeSplicingAnalysis$IR)
    print(ir_table)
    
    # Overview of other alternative splicing events
    cat("\nAlternative splicing event summary:\n")
    
    # Check for different types of events
    as_events <- c("IR", "A5", "A3", "ES", "MES", "ATSS", "ATTS")
    as_summary <- list()
    
    for(event in as_events) {
        if(event %in% colnames(switchListAnalyzed$AlternativeSplicingAnalysis)) {
            event_count <- sum(switchListAnalyzed$AlternativeSplicingAnalysis[[event]] > 0, na.rm = TRUE)
            as_summary[[event]] <- event_count
        }
    }
    
    cat("Alternative splicing events detected:\n")
    for(event in names(as_summary)) {
        cat("-", event, ":", as_summary[[event]], "isoforms\n")
    }
    
} else {
    cat("Alternative splicing analysis not found in switchList\n")
}
# Intron Retention (IR) events:

#    0    1    2    3    4    5    6    9   10 
# 3665 1104  297   78   16    6    2    1    1 

# Alternative splicing event summary:
# Alternative splicing events detected:
# - IR : 1505 isoforms
# - A5 : 2069 isoforms
# - A3 : 2445 isoforms
# - ES : 2420 isoforms
# - MES : 1510 isoforms
# - ATSS : 3286 isoforms
# - ATTS : 2988 isoforms


#### Analyze Switch Consequences

# Re-extract sequences to ensure amino acid sequences are available
switchListAnalyzed <- extractSequence(
    switchAnalyzeRlist = switchListAnalyzed,
    pathToOutput = "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom",
    writeToFile = TRUE,
    removeLongAAseq = FALSE,  # Keep all AA sequences for analysis
    alsoSplitFastaFile = FALSE  # Don't split for now
)

# Define consequences of interest
# Focus on the most relevant functional consequences
valid_consequences <- c(
    'intron_retention',        # Alternative splicing - intron retention
    'coding_potential',        # Coding vs non-coding switches  
    'NMD_status',             # Nonsense-mediated decay
    'domains_identified',      # Protein domain changes (from Pfam)
    'ORF_seq_similarity',     # Protein sequence similarity
    'ORF_length',             # ORF length changes
    'last_exon',              # Last exon changes
    'isoform_seq_similarity', # Isoform sequence similarity
    'exon_number',            # Number of exons
    'isoform_length'          # Isoform length
)

# Check which consequences can be analyzed with the current data
available_consequences <- c()

# Alternative splicing consequences
if("AlternativeSplicingAnalysis" %in% names(switchListAnalyzed)) {
    available_consequences <- c(available_consequences, "intron_retention", "last_exon")
}

# Coding potential consequences
if("codingPotential" %in% colnames(switchListAnalyzed$isoformFeatures)) {
    available_consequences <- c(available_consequences, "coding_potential", "NMD_status")
}

# Pfam domain consequences
if("pfamAnalysis" %in% names(switchListAnalyzed)) {
    available_consequences <- c(available_consequences, "domains_identified")
}

# ORF-based consequences
if("orfAnalysis" %in% names(switchListAnalyzed)) {
    available_consequences <- c(available_consequences, "ORF_seq_similarity", "ORF_length")
}

# Basic isoform consequences (always available)
available_consequences <- c(available_consequences, "isoform_seq_similarity", "exon_number", "isoform_length")

# Use only valid and available consequences
final_consequences <- intersect(valid_consequences, available_consequences)
cat("Valid consequences to analyze:", paste(final_consequences, collapse = ", "), "\n")

#### Perform switch consequence analysis with valid parameters
# isoforms are divided into the isoforms that increase their contribution to gene expression (positive dIF values larger than dIFcutoff) and the isoforms that decrease their contribution (negative dIF values smaller than -dIFcutoff). 
# The isoforms with increased contribution are then (in a pairwise manner) compared to the isoform with decreasing contribution. In each of these comparisons the isoforms compared are analyzed for differences in their annotation 
switchListAnalyzed <- analyzeSwitchConsequences(
    switchAnalyzeRlist = switchListAnalyzed,
    consequencesToAnalyze = final_consequences,
    dIFcutoff = 0.1,      # Standard cutoff (10% change in isoform fraction)
    alpha = 0.05,         # FDR cutoff
    ntCutoff = 50,        # Nucleotide length difference cutoff
    AaCutoff = 10,        # Amino acid length difference cutoff
    showProgress = TRUE
)

##### Extract and summarize switch consequences
# Summary without filtering for consequences
summary_all <- extractSwitchSummary(
    switchListAnalyzed, 
    dIFcutoff = 0.1, 
    filterForConsequences = FALSE
)
cat("All isoform switches (dIF > 0.1):\n")
print(summary_all)
#       Comparison nrIsoforms nrSwitches nrGenes
# 1 control vs t30        627        821     513
# 2 control vs t35         85        155      76
# 3     t30 vs t35        466        659     372
# 4       Combined       1099       1578     871


# Summary filtering for functional consequences
summary_consequences <- extractSwitchSummary(
    switchListAnalyzed, 
    dIFcutoff = 0.1, 
    filterForConsequences = TRUE
)
cat("\nIsoform switches with functional consequences:\n")
print(summary_consequences)
#       Comparison nrIsoforms nrSwitches nrGenes
# 1 control vs t30        597        806     485
# 2 control vs t35         84        154      75
# 3     t30 vs t35        456        657     363
# 4       Combined       1062       1561     837

# Calculate percentage of switches with consequences
if(nrow(summary_all) > 0 && nrow(summary_consequences) > 0) {
    # Calculate percentages for each comparison
    for(i in 1:nrow(summary_all)) {
        comparison <- summary_all$Comparison[i]
        total_switches <- summary_all$nrSwitches[i]
        total_genes <- summary_all$nrGenes[i]
        
        consequence_row <- which(summary_consequences$Comparison == comparison)
        if(length(consequence_row) > 0) {
            consequence_switches <- summary_consequences$nrSwitches[consequence_row]
            consequence_genes <- summary_consequences$nrGenes[consequence_row]
            
            switch_pct <- round((consequence_switches / total_switches) * 100, 1)
            gene_pct <- round((consequence_genes / total_genes) * 100, 1)
            
            cat(sprintf("%s: %d/%d switches (%s%%) in %d/%d genes (%s%%) have functional consequences\n",
                       comparison, consequence_switches, total_switches, switch_pct,
                       consequence_genes, total_genes, gene_pct))
        }
    }
}
# control vs t30: 806/821 switches (98.2%) in 485/513 genes (94.5%) have functional consequences
# control vs t35: 154/155 switches (99.4%) in 75/76 genes (98.7%) have functional consequences
# t30 vs t35: 657/659 switches (99.7%) in 363/372 genes (97.6%) have functional consequences
# Combined: 1561/1578 switches (98.9%) in 837/871 genes (96.1%) have functional consequences

# Get detailed consequence results
if("switchConsequence" %in% names(switchListAnalyzed)) {
    consequence_details <- switchListAnalyzed$switchConsequence
    
    cat("\nDetailed consequence breakdown:\n")
    cat("Total consequence entries:", nrow(consequence_details), "\n")
    
    # Count different types of consequences
    if("featureCompared" %in% colnames(consequence_details)) {
        consequence_types <- table(consequence_details$featureCompared)
        cat("Consequence types:\n")
        print(consequence_types)
    }
    
    # Show sample consequences
    cat("\nSample functional consequences:\n")
    print(head(consequence_details, 3))
}


############ Analysis of Individual Isoform Switching

# Extract top switching genes by q-value (most statistically significant)
top_genes_qval <- extractTopSwitches(
    switchListAnalyzed, 
    filterForConsequences = TRUE,  # Only genes with functional consequences
    n = 10,                       # Top 10 genes
    sortByQvals = TRUE,           # Sort by statistical significance
    extractGenes = TRUE           # Extract gene-level results
)
print(top_genes_qval[, c("gene_id", "gene_name", "condition_1", "condition_2", "gene_switch_q_value", "switchConsequencesGene")])

# Extract top switching genes by effect size (largest dIF changes)
top_genes_dif <- extractTopSwitches(
    switchListAnalyzed, 
    filterForConsequences = TRUE,  # Only genes with functional consequences
    n = 10,                       # Top 10 genes
    sortByQvals = FALSE,          # Sort by effect size
    extractGenes = TRUE           # Extract gene-level results
)
print(top_genes_dif[, c("gene_id", "gene_name", "condition_1", "condition_2", "combinedDIF", "switchConsequencesGene")])



### Get top 10 switching genes by effect size for each temperature comparison separately

comparisons <- c("control vs t30", "control vs t35", "t30 vs t35")

# Store results for each comparison
top_genes_by_comparison <- list()

# First, get ALL switching genes with consequences
all_switching_genes <- extractTopSwitches(
    switchListAnalyzed,
    filterForConsequences = TRUE,  # Only genes with functional consequences
    n = Inf,                      # Get all significant switches
    sortByQvals = FALSE,          # Sort by effect size
    extractGenes = TRUE           # Extract gene-level results
)

cat("Total switching genes found:", nrow(all_switching_genes), "\n")


# Now filter by each comparison
for(comparison in comparisons) {
    cat(sprintf("\n=== %s ===\n", comparison))
    
    # Filter for this specific comparison
    comp_genes <- all_switching_genes[
        paste(all_switching_genes$condition_1, "vs", all_switching_genes$condition_2) == comparison, 
    ]
    
    if(nrow(comp_genes) > 0) {
        # Sort by effect size and take top 10
        comp_genes <- comp_genes[order(comp_genes$combinedDIF, decreasing = TRUE), ]
        top_comp_dif <- comp_genes[1:min(10, nrow(comp_genes)), ]
        
        cat(sprintf("Top 10 genes with largest effect sizes for %s:\n", comparison))
        
        # Display the results with key columns
        display_cols <- c("gene_id", "gene_name", "combinedDIF", "gene_switch_q_value", "switchConsequencesGene")
        print(top_comp_dif[, display_cols])
        
        # Store results for later use
        top_genes_by_comparison[[comparison]] <- top_comp_dif
        
        # Summary statistics for this comparison
        cat(sprintf("\nSummary for %s:\n", comparison))
        cat(sprintf("- Number of genes found: %d\n", nrow(comp_genes)))
        cat(sprintf("- Top 10 genes shown: %d\n", nrow(top_comp_dif)))
        cat(sprintf("- Largest effect size: %.3f\n", max(top_comp_dif$combinedDIF, na.rm = TRUE)))
        cat(sprintf("- Smallest effect size in top 10: %.3f\n", min(top_comp_dif$combinedDIF, na.rm = TRUE)))
        cat(sprintf("- Mean effect size (top 10): %.3f\n", mean(top_comp_dif$combinedDIF, na.rm = TRUE)))
        
        # Show genes with names vs unknown
        named_genes <- sum(!is.na(top_comp_dif$gene_name) & top_comp_dif$gene_name != "")
        cat(sprintf("- Genes with annotation: %d/%d\n", named_genes, nrow(top_comp_dif)))
        
    } else {
        cat(sprintf("No significant switches found for %s\n", comparison))
    }
    
    cat("\n", paste(rep("-", 80), collapse = ""), "\n")
}

# Create a comparative summary table
if(length(top_genes_by_comparison) > 0) {
    # Create summary table
    comparison_summary <- data.frame()
    
    for(comp_name in names(top_genes_by_comparison)) {
        comp_data <- top_genes_by_comparison[[comp_name]]
        
        if(nrow(comp_data) > 0) {
            summary_row <- data.frame(
                Comparison = comp_name,
                Max_Effect_Size = round(max(comp_data$combinedDIF, na.rm = TRUE), 3),
                Min_Effect_Size = round(min(comp_data$combinedDIF, na.rm = TRUE), 3),
                Mean_Effect_Size = round(mean(comp_data$combinedDIF, na.rm = TRUE), 3),
                Genes_with_Names = sum(!is.na(comp_data$gene_name) & comp_data$gene_name != ""),
                Most_Significant_Gene = comp_data$gene_id[which.min(comp_data$gene_switch_q_value)],
                Largest_Effect_Gene = comp_data$gene_id[which.max(comp_data$combinedDIF)]
            )
            comparison_summary <- rbind(comparison_summary, summary_row)
        }
    }
    
    print(comparison_summary)
}

# Identify genes that appear in multiple comparisons

all_top_genes <- c()
for(comp_name in names(top_genes_by_comparison)) {
    comp_genes <- top_genes_by_comparison[[comp_name]]$gene_id
    all_top_genes <- c(all_top_genes, comp_genes)
}

# Find genes that appear multiple times
gene_counts <- table(all_top_genes)
multi_comparison_genes <- gene_counts[gene_counts > 1]

if(length(multi_comparison_genes) > 0) {
    cat("Genes appearing in multiple top-10 lists:\n")
    for(gene in names(multi_comparison_genes)) {
        count <- multi_comparison_genes[gene]
        
        # Get gene name
        gene_name <- "Unknown"
        for(comp_data in top_genes_by_comparison) {
            if(gene %in% comp_data$gene_id) {
                gene_name <- comp_data$gene_name[comp_data$gene_id == gene][1]
                if(!is.na(gene_name) && gene_name != "") break
            }
        }
        
        cat(sprintf("- %s (%s): appears in %d comparisons\n", gene, gene_name, count))
        
        # Show which comparisons
        appearances <- c()
        for(comp_name in names(top_genes_by_comparison)) {
            if(gene %in% top_genes_by_comparison[[comp_name]]$gene_id) {
                comp_data <- top_genes_by_comparison[[comp_name]]
                gene_row <- comp_data[comp_data$gene_id == gene, ]
                effect_size <- round(gene_row$combinedDIF, 3)
                appearances <- c(appearances, paste0(comp_name, " (dIF=", effect_size, ")"))
            }
        }
        cat(sprintf("  Comparisons: %s\n", paste(appearances, collapse = ", ")))
    }
} else {
    cat("No genes appear in multiple top-10 lists\n")
}

# Save results to files
for(comp_name in names(top_genes_by_comparison)) {
    file_name <- gsub(" vs ", "_vs_", comp_name)
    file_name <- gsub(" ", "_", file_name)
    
    write.csv(top_genes_by_comparison[[comp_name]], 
              paste0("/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/top10_genes_effect_size_", file_name, ".csv"),
              row.names = FALSE)
}





# Get top 10 switching ISOFORMS by effect size for each temperature comparison separately

comparisons <- c("control vs t30", "control vs t35", "t30 vs t35")

# Store results for each comparison
top_isoforms_by_comparison <- list()

# First, get ALL switching isoforms with consequences
all_switching_isoforms_effect <- extractTopSwitches(
    switchListAnalyzed,
    filterForConsequences = TRUE,  # Only isoforms with functional consequences
    n = Inf,                      # Get all significant switches
    sortByQvals = FALSE,          # Sort by effect size (dIF)
    extractGenes = FALSE          # Extract ISOFORM-level results (not genes)
)

cat("Total switching isoforms found:", nrow(all_switching_isoforms_effect), "\n")

# Now filter by each comparison
for(comparison in comparisons) {
    cat(sprintf("\n=== %s ===\n", comparison))
    
    # Filter for this specific comparison
    comp_isoforms <- all_switching_isoforms_effect[
        paste(all_switching_isoforms_effect$condition_1, "vs", all_switching_isoforms_effect$condition_2) == comparison, 
    ]
    
    if(nrow(comp_isoforms) > 0) {
        # Sort by effect size (dIF) and take top 10
        comp_isoforms <- comp_isoforms[order(abs(comp_isoforms$dIF), decreasing = TRUE), ]  # Use absolute dIF for ranking
        top_comp_isoforms <- comp_isoforms[1:min(10, nrow(comp_isoforms)), ]
        
        cat(sprintf("Top 10 isoforms with largest effect sizes for %s:\n", comparison))
        
        # Display the results with key columns
        display_cols <- c("isoform_id", "gene_id", "gene_name", "dIF", "IF1", "IF2", "isoform_switch_q_value")
        print(top_comp_isoforms[, display_cols])
        
        # Store results for later use
        top_isoforms_by_comparison[[comparison]] <- top_comp_isoforms
        
        # Summary statistics for this comparison
        cat(sprintf("\nSummary for %s:\n", comparison))
        cat(sprintf("- Number of isoforms found: %d\n", nrow(comp_isoforms)))
        cat(sprintf("- Top 10 isoforms shown: %d\n", nrow(top_comp_isoforms)))
        cat(sprintf("- Largest effect size (|dIF|): %.3f\n", max(abs(top_comp_isoforms$dIF), na.rm = TRUE)))
        cat(sprintf("- Smallest effect size in top 10 (|dIF|): %.3f\n", min(abs(top_comp_isoforms$dIF), na.rm = TRUE)))
        cat(sprintf("- Mean effect size (|dIF|): %.3f\n", mean(abs(top_comp_isoforms$dIF), na.rm = TRUE)))
        
        # Show isoforms with gene names vs unknown
        named_isoforms <- sum(!is.na(top_comp_isoforms$gene_name) & top_comp_isoforms$gene_name != "")
        cat(sprintf("- Isoforms with gene annotation: %d/%d\n", named_isoforms, nrow(top_comp_isoforms)))
        
        # Show direction of changes
        positive_dif <- sum(top_comp_isoforms$dIF > 0)
        negative_dif <- sum(top_comp_isoforms$dIF < 0)
        cat(sprintf("- Increased usage (positive dIF): %d isoforms\n", positive_dif))
        cat(sprintf("- Decreased usage (negative dIF): %d isoforms\n", negative_dif))
        
    } else {
        cat(sprintf("No significant switching isoforms found for %s\n", comparison))
    }
    
    cat("\n", paste(rep("-", 80), collapse = ""), "\n")
}

# Create a comparative summary table for isoforms

if(length(top_isoforms_by_comparison) > 0) {
    # Create summary table
    isoform_comparison_summary <- data.frame()
    
    for(comp_name in names(top_isoforms_by_comparison)) {
        comp_data <- top_isoforms_by_comparison[[comp_name]]
        
        if(nrow(comp_data) > 0) {
            summary_row <- data.frame(
                Comparison = comp_name,
                Max_Effect_Size = round(max(abs(comp_data$dIF), na.rm = TRUE), 3),
                Min_Effect_Size = round(min(abs(comp_data$dIF), na.rm = TRUE), 3),
                Mean_Effect_Size = round(mean(abs(comp_data$dIF), na.rm = TRUE), 3),
                Isoforms_with_Names = sum(!is.na(comp_data$gene_name) & comp_data$gene_name != ""),
                Most_Significant_Isoform = comp_data$isoform_id[which.min(comp_data$isoform_switch_q_value)],
                Largest_Effect_Isoform = comp_data$isoform_id[which.max(abs(comp_data$dIF))],
                Positive_Changes = sum(comp_data$dIF > 0),
                Negative_Changes = sum(comp_data$dIF < 0)
            )
            isoform_comparison_summary <- rbind(isoform_comparison_summary, summary_row)
        }
    }
    
    print(isoform_comparison_summary)
}

# Identify isoforms that appear in multiple comparisons

all_top_isoforms <- c()
for(comp_name in names(top_isoforms_by_comparison)) {
    comp_isoforms <- top_isoforms_by_comparison[[comp_name]]$isoform_id
    all_top_isoforms <- c(all_top_isoforms, comp_isoforms)
}

# Find isoforms that appear multiple times
isoform_counts <- table(all_top_isoforms)
multi_comparison_isoforms <- isoform_counts[isoform_counts > 1]

if(length(multi_comparison_isoforms) > 0) {
    cat("Isoforms appearing in multiple top-10 lists:\n")
    for(isoform in names(multi_comparison_isoforms)) {
        count <- multi_comparison_isoforms[isoform]
        
        # Get gene name
        gene_name <- "Unknown"
        gene_id <- "Unknown"
        for(comp_data in top_isoforms_by_comparison) {
            if(isoform %in% comp_data$isoform_id) {
                isoform_row <- comp_data[comp_data$isoform_id == isoform, ][1, ]
                gene_name <- ifelse(!is.na(isoform_row$gene_name) && isoform_row$gene_name != "", 
                                   isoform_row$gene_name, "Unknown")
                gene_id <- isoform_row$gene_id
                break
            }
        }
        
        cat(sprintf("- %s (Gene: %s, %s): appears in %d comparisons\n", 
                   isoform, gene_id, gene_name, count))
        
        # Show which comparisons and effect sizes
        appearances <- c()
        for(comp_name in names(top_isoforms_by_comparison)) {
            if(isoform %in% top_isoforms_by_comparison[[comp_name]]$isoform_id) {
                comp_data <- top_isoforms_by_comparison[[comp_name]]
                isoform_row <- comp_data[comp_data$isoform_id == isoform, ]
                effect_size <- round(isoform_row$dIF, 3)
                appearances <- c(appearances, paste0(comp_name, " (dIF=", effect_size, ")"))
            }
        }
        cat(sprintf("  Comparisons: %s\n", paste(appearances, collapse = ", ")))
    }
} else {
    cat("No isoforms appear in multiple top-10 lists\n")
}

# Identify genes with multiple top isoforms in the same comparison

for(comp_name in names(top_isoforms_by_comparison)) {
    comp_data <- top_isoforms_by_comparison[[comp_name]]
    
    # Count isoforms per gene in this comparison
    gene_isoform_counts <- table(comp_data$gene_id)
    multi_isoform_genes_in_comp <- gene_isoform_counts[gene_isoform_counts > 1]
    
    if(length(multi_isoform_genes_in_comp) > 0) {
        cat(sprintf("\n%s - Genes with multiple top-10 isoforms:\n", comp_name))
        for(gene in names(multi_isoform_genes_in_comp)) {
            gene_isoforms <- comp_data[comp_data$gene_id == gene, ]
            gene_name <- unique(gene_isoforms$gene_name)[1]
            gene_name <- ifelse(is.na(gene_name) || gene_name == "", "Unknown", gene_name)
            
            cat(sprintf("- %s (%s): %d isoforms\n", gene, gene_name, multi_isoform_genes_in_comp[gene]))
            
            # Show the isoforms and their effects
            for(i in 1:nrow(gene_isoforms)) {
                iso_row <- gene_isoforms[i, ]
                cat(sprintf("  %s: dIF=%.3f, IF1=%.3f, IF2=%.3f\n", 
                           iso_row$isoform_id, iso_row$dIF, iso_row$IF1, iso_row$IF2))
            }
        }
    } else {
        cat(sprintf("\n%s - No genes with multiple top-10 isoforms\n", comp_name))
    }
}

# Save results to files
for(comp_name in names(top_isoforms_by_comparison)) {
    file_name <- gsub(" vs ", "_vs_", comp_name)
    file_name <- gsub(" ", "_", file_name)
    
    write.csv(top_isoforms_by_comparison[[comp_name]], 
              paste0("/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/top10_isoforms_effect_size_", file_name, ".csv"),
              row.names = FALSE)
}



# Generate switch plots for top switching genes per temperature comparison
# Create output directory for switch plots
plot_output_dir <- "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/isoforms_plots/"
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

# Generate plots for each temperature comparison separately (top 1 gene each)
comparisons <- c("control vs t30", "control vs t35", "t30 vs t35")

for(comparison in comparisons) {
    cat(sprintf("\nGenerating plot for top 1 gene in %s...\n", comparison))
    
    # Create comparison-specific directory
    comp_dir <- gsub(" vs ", "_vs_", comparison)
    comp_dir <- gsub(" ", "_", comp_dir)
    comparison_output_dir <- paste0(plot_output_dir, comp_dir, "/")
    dir.create(comparison_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Filter isoformFeatures for this comparison
    comparison_isoforms <- switchListAnalyzed$isoformFeatures[
        paste(switchListAnalyzed$isoformFeatures$condition_1, "vs", switchListAnalyzed$isoformFeatures$condition_2) == comparison, 
    ]
    
    if(nrow(comparison_isoforms) > 0) {
        # Create a subset switchList for this comparison
        switchList_subset <- switchListAnalyzed
        switchList_subset$isoformFeatures <- comparison_isoforms
        
        # Also subset other relevant components if they exist
        if("switchTestResults" %in% names(switchListAnalyzed)) {
            subset_switch_results <- switchListAnalyzed$switchTestResults[
                paste(switchListAnalyzed$switchTestResults$condition_1, "vs", switchListAnalyzed$switchTestResults$condition_2) == comparison, 
            ]
            switchList_subset$switchTestResults <- subset_switch_results
        }
        
        # Generate plot for TOP 1 gene in this comparison
        switchPlotTopSwitches(
            switchAnalyzeRlist = switchList_subset,
            n = 1,                        # TOP 1 GENE ONLY
            filterForConsequences = TRUE,  # Only genes with functional consequences
            fileType = "pdf",             # PDF format
            pathToOutput = comparison_output_dir,
            sortByQvals = FALSE           # Sort by effect size (largest dIF)
        )
        
        cat(sprintf("Top 1 gene plot saved to: %s\n", comparison_output_dir))
        
        # Also create PNG version for easy viewing
        switchPlotTopSwitches(
            switchAnalyzeRlist = switchList_subset,
            n = 1,                        # TOP 1 GENE ONLY
            filterForConsequences = TRUE,
            fileType = "png",             # PNG format for easy viewing
            pathToOutput = comparison_output_dir,
            sortByQvals = FALSE
        )
        
    } else {
        cat(sprintf("No switching isoforms found for %s\n", comparison))
    }
}


# Create summary of which genes were plotted with PFAM DOMAINS
for(comp_name in names(top_genes_by_comparison)) {
    if(nrow(top_genes_by_comparison[[comp_name]]) > 0) {
        top_gene <- top_genes_by_comparison[[comp_name]][1, ]  # Get the #1 ranked gene
        
        gene_id <- top_gene$gene_id
        gene_name <- ifelse(is.na(top_gene$gene_name) || top_gene$gene_name == "", "Unknown", top_gene$gene_name)
        effect_size <- round(top_gene$combinedDIF, 3)
        q_value <- top_gene$gene_switch_q_value
        
        cat(sprintf("%s:\n", comp_name))
        cat(sprintf("  Gene ID: %s\n", gene_id))
        cat(sprintf("  Gene Name: %s\n", gene_name))
        cat(sprintf("  Effect Size (combinedDIF): %s\n", effect_size))
        cat(sprintf("  Q-value: %.2e\n", q_value))
        
        # Check if this gene has Pfam domains
        if("pfamAnalysis" %in% names(switchListAnalyzed)) {
            gene_pfam <- switchListAnalyzed$pfamAnalysis[
                switchListAnalyzed$pfamAnalysis$isoform_id %in% 
                switchListAnalyzed$isoformFeatures$isoform_id[
                    switchListAnalyzed$isoformFeatures$gene_id == gene_id], ]
            
            if(nrow(gene_pfam) > 0) {
                domains <- unique(gene_pfam$pfamName)
                cat(sprintf("  Pfam domains: %s\n", paste(domains, collapse = ", ")))
                cat(sprintf("  Domain tracks will show: %d different domains\n", length(domains)))
            } else {
                cat("  No Pfam domains found for this gene\n")
            }
        }
        
        cat(sprintf("  Plot files: %s_1.pdf and %s_1.png\n", 
                   gsub(" vs ", "_vs_", comp_name), gsub(" vs ", "_vs_", comp_name)))
        cat("\n")
    }
}





######### Genome-Wide Analysis of Alternative Splicing

# 1. EXTRACT SPLICING SUMMARY
# Basic splicing summary (absolute numbers)
splicing_summary_absolute <- extractSplicingSummary(
    switchListAnalyzed,
    asFractionTotal = FALSE,  # Show absolute numbers
    plotGenes = FALSE         # Count events, not genes
)

cat("Alternative splicing events (absolute numbers):\n")
print(splicing_summary_absolute)

# Splicing summary as fractions of total
splicing_summary_fraction <- extractSplicingSummary(
    switchListAnalyzed,
    asFractionTotal = TRUE,   # Show as fractions
    plotGenes = FALSE         # Count events, not genes
)

cat("\nAlternative splicing events (as fractions):\n")
print(splicing_summary_fraction)

# # Gene-level splicing summary
# splicing_summary_genes <- extractSplicingSummary(
#     switchListAnalyzed,
#     asFractionTotal = FALSE,  # Show absolute numbers
#     plotGenes = TRUE          # Count genes with events
# )

# cat("\nGenes with alternative splicing events:\n")
# print(splicing_summary_genes)

# 2. SPLICING ENRICHMENT ANALYSIS
# Test for enrichment of splicing events in switching isoforms vs all isoforms
splicing_enrichment <- extractSplicingEnrichment(
    switchListAnalyzed,
    splicingToAnalyze = 'all',     # Analyze all types of splicing
    alpha = 0.05,                  # Significance threshold of enrichment test
    dIFcutoff = 0.1,              # dIF cutoff for switching, difference in isoform fraction (dIF) - effect size
    onlySigIsoforms = FALSE,     # Isee below for explanation
    returnResult = TRUE            # Return results as data frame
)

cat("Splicing enrichment in switching isoforms:\n")
print(splicing_enrichment)
#save as 4x16 landscape 	splicing_enrichment_allisoforms_Pcom.pdf

# onlySigIsoforms: A logic indicating whether to only consider significant isoforms, meaning only analyzing genes where at least 
# two isoforms which both have significant usage changes in opposite direction (quite strict) Naturally this only works if the isoform switch test used have isoform resolution 
# (which the build in isoformSwitchTestDEXSeq has). If FALSE all isoforms with an absolute dIF value larger than dIFcutoff in a gene with significant switches (defined by alpha and dIFcutoff) 
# are included in the pairwise comparison. Default is FALSE (non significant isoforms are also considered based on the logic that if one isoform changes it contribution - there must be an equivalent 
# opposite change in usage in the other isoforms from that gene).


############# Global Analysis of Alternative Splicing Patterns (All Isoforms)

# get basic statistics about complete dataset
total_isoforms <- nrow(switchListAnalyzed$isoformFeatures)
#19130
total_genes <- length(unique(switchListAnalyzed$isoformFeatures$gene_id))
#969

##  Global Alternative Splicing Event Counts
if("AlternativeSplicingAnalysis" %in% names(switchListAnalyzed)) {
    as_data <- switchListAnalyzed$AlternativeSplicingAnalysis
    
    # Basic dataset info
    total_isoforms <- nrow(switchListAnalyzed$isoformFeatures)
    total_genes <- length(unique(switchListAnalyzed$isoformFeatures$gene_id))
    total_as_isoforms <- nrow(as_data)
    
    cat("Dataset Overview:\n")
    cat("- Total isoforms analyzed:", total_isoforms, "\n")
    cat("- Total genes analyzed:", total_genes, "\n")
    cat("- Isoforms with AS analysis:", total_as_isoforms, "\n")
    cat("- Percentage with AS data:", round((total_as_isoforms / total_isoforms) * 100, 1), "%\n")
    
} else {
    cat("Alternative splicing analysis not available\n")
    stop("Cannot proceed without alternative splicing data")
}

# GLOBAL ALTERNATIVE SPLICING EVENT COUNTS
splicing_events <- c("IR", "A5", "A3", "ES", "MES", "ATSS", "ATTS")

global_as_summary <- data.frame()

for(event in splicing_events) {
    if(event %in% colnames(as_data)) {
        # Count isoforms with this event
        isoforms_with_event <- sum(as_data[[event]] > 0, na.rm = TRUE)
        # Total events
        total_events <- sum(as_data[[event]], na.rm = TRUE)
        # Percentage affected
        pct_isoforms <- round((isoforms_with_event / total_as_isoforms) * 100, 1)
        # Average events per affected isoform
        avg_events <- ifelse(isoforms_with_event > 0, 
                            round(total_events / isoforms_with_event, 2), 0)
        
        global_as_summary <- rbind(global_as_summary, data.frame(
            Event_Type = event,
            Isoforms_Affected = isoforms_with_event,
            Total_Events = total_events,
            Percentage_Isoforms = pct_isoforms,
            Events_Per_Isoform = avg_events
        ))
    }
}

print(global_as_summary)
#   Event_Type Isoforms_Affected Total_Events Percentage_Isoforms Events_Per_Isoform
# 1         IR              1505         2057                29.1               1.37
# 2         A5              2069         2688                40.0               1.30
# 3         A3              2445         3687                47.3               1.51
# 4         ES              2420         3513                46.8               1.45
# 5        MES              1510         1799                29.2               1.19
# 6       ATSS              3286         3286                63.6               1.00
# 7       ATTS              2988         2988                57.8               1.00

write.csv(global_as_summary, "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/global_splicing_events.csv", row.names = FALSE)

# #### Switching vs Non-Switching Comparison  
# # Get switching isoform IDs
# switching_isoform_ids <- unique(switchListAnalyzed$isoformFeatures$isoform_id[
#     switchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05 & 
#     abs(switchListAnalyzed$isoformFeatures$dIF) >= 0.1
# ])

# # Split data
# switching_as <- as_data[as_data$isoform_id %in% switching_isoform_ids, ]
# non_switching_as <- as_data[!as_data$isoform_id %in% switching_isoform_ids, ]

# cat("Data split:\n")
# cat("- Switching isoforms with AS data:", nrow(switching_as), "\n")
# cat("- Non-switching isoforms with AS data:", nrow(non_switching_as), "\n")

# # Create comparison table
# comparison_summary <- data.frame()

# for(event in splicing_events) {
#     if(event %in% colnames(as_data)) {
#         # Switching isoforms
#         switch_affected <- sum(switching_as[[event]] > 0, na.rm = TRUE)
#         switch_pct <- round((switch_affected / nrow(switching_as)) * 100, 1)
        
#         # Non-switching isoforms
#         non_switch_affected <- sum(non_switching_as[[event]] > 0, na.rm = TRUE)
#         non_switch_pct <- round((non_switch_affected / nrow(non_switching_as)) * 100, 1)
        
#         # Enrichment ratio
#         enrichment <- ifelse(non_switch_pct > 0, round(switch_pct / non_switch_pct, 2), NA)
        
#         comparison_summary <- rbind(comparison_summary, data.frame(
#             Event_Type = event,
#             Switching_Affected = switch_affected,
#             Switching_Percent = switch_pct,
#             NonSwitching_Affected = non_switch_affected,
#             NonSwitching_Percent = non_switch_pct,
#             Enrichment_Ratio = enrichment
#         ))
#     }
# }

# print(comparison_summary)




# INTEGRATION: WGCNA HUB GENES vs ISOFORM SWITCHING ANALYSIS

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(VennDiagram)

# 1. EXTRACT SWITCHING GENES FROM ISOFORM ANALYSIS - FILTERED BY WGCNA GENES

# Get all genes used in WGCNA analysis
wgcna_genes <- colnames(datExpr)  # All genes in WGCNA expression matrix
#28378 genes

# Get switching genes for each comparison - FILTERED by WGCNA genes
switching_genes_control_vs_t30_all <- unique(switchListAnalyzed$isoformFeatures$gene_id[
    paste(switchListAnalyzed$isoformFeatures$condition_1, "vs", 
          switchListAnalyzed$isoformFeatures$condition_2) == "control vs t30" &
    switchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05 &
    abs(switchListAnalyzed$isoformFeatures$dIF) >= 0.1
])

switching_genes_control_vs_t35_all <- unique(switchListAnalyzed$isoformFeatures$gene_id[
    paste(switchListAnalyzed$isoformFeatures$condition_1, "vs", 
          switchListAnalyzed$isoformFeatures$condition_2) == "control vs t35" &
    switchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05 &
    abs(switchListAnalyzed$isoformFeatures$dIF) >= 0.1
])

switching_genes_t30_vs_t35_all <- unique(switchListAnalyzed$isoformFeatures$gene_id[
    paste(switchListAnalyzed$isoformFeatures$condition_1, "vs", 
          switchListAnalyzed$isoformFeatures$condition_2) == "t30 vs t35" &
    switchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05 &
    abs(switchListAnalyzed$isoformFeatures$dIF) >= 0.1
])

# FILTER switching genes to only include those in WGCNA analysis
switching_genes_control_vs_t30 <- intersect(switching_genes_control_vs_t30_all, wgcna_genes)
switching_genes_control_vs_t35 <- intersect(switching_genes_control_vs_t35_all, wgcna_genes)
switching_genes_t30_vs_t35 <- intersect(switching_genes_t30_vs_t35_all, wgcna_genes)

# All switching genes combined 
all_switching_genes_wgcna <- unique(c(switching_genes_control_vs_t30, 
                                     switching_genes_control_vs_t35, 
                                     switching_genes_t30_vs_t35))

# Print comparison of switching genes before and after filtering
cat("Switching genes summary (before vs after WGCNA filtering):\n")
cat("Control vs 30°C: ", length(switching_genes_control_vs_t30_all), " -> ", length(switching_genes_control_vs_t30), "\n")
#Control vs 30°C:  513  ->  457 
cat("Control vs 35°C: ", length(switching_genes_control_vs_t35_all), " -> ", length(switching_genes_control_vs_t35), "\n")
#Control vs 35°C:  77  ->  69
cat("30°C vs 35°C: ", length(switching_genes_t30_vs_t35_all), " -> ", length(switching_genes_t30_vs_t35), "\n")
#30°C vs 35°C:  373  ->  338 
cat("Total unique switching genes: ", length(unique(c(switching_genes_control_vs_t30_all, switching_genes_control_vs_t35_all, switching_genes_t30_vs_t35_all))), " -> ", length(all_switching_genes_wgcna), "\n")
#Total unique switching genes:  872  ->  785 

# Calculate what percentage of switching genes are in WGCNA
all_switching_original <- unique(c(switching_genes_control_vs_t30_all, switching_genes_control_vs_t35_all, switching_genes_t30_vs_t35_all))
overlap_pct <- round((length(all_switching_genes_wgcna) / length(all_switching_original)) * 100, 1)
cat("Percentage of switching genes in WGCNA: ", overlap_pct, "%\n")
#Percentage of switching genes in WGCNA:  90 %


# 2. LOAD HUB GENES FROM WGCNA ANALYSIS

#cluster1
cluster1_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster1.csv")
cluster1_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster1.csv")
cluster1_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster1.csv")

#cluster2
cluster2_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster2.csv")
cluster2_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster2.csv")
cluster2_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster2.csv")

#cluster3
cluster3_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster3.csv")
cluster3_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster3.csv")
cluster3_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster3.csv")

#cluster9
cluster9_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster9.csv")
cluster9_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster9.csv")
cluster9_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster9.csv")

# Extract gene vectors from the CSV files you loaded
cluster1_hub_genes_control <- cluster1_hub_gene_counts_control$Gene  
cluster1_hub_genes_30 <- cluster1_hub_gene_counts_30$Gene
cluster1_hub_genes_35 <- cluster1_hub_gene_counts_35$Gene

cluster2_hub_genes_control <- cluster2_hub_gene_counts_control$Gene
cluster2_hub_genes_30 <- cluster2_hub_gene_counts_30$Gene
cluster2_hub_genes_35 <- cluster2_hub_gene_counts_35$Gene

cluster3_hub_genes_control <- cluster3_hub_gene_counts_control$Gene
cluster3_hub_genes_30 <- cluster3_hub_gene_counts_30$Gene
cluster3_hub_genes_35 <- cluster3_hub_gene_counts_35$Gene

cluster9_hub_genes_control <- cluster9_hub_gene_counts_control$Gene
cluster9_hub_genes_30 <- cluster9_hub_gene_counts_30$Gene
cluster9_hub_genes_35 <- cluster9_hub_gene_counts_35$Gene

# Get unique hub genes for each temp (non-shared genes)
cluster1_hub_genes_30_unique <- setdiff(
  cluster1_hub_genes_30,
  union(cluster1_hub_genes_control, cluster1_hub_genes_35)
)

cluster1_hub_genes_35_unique <- setdiff(
  cluster1_hub_genes_35,
  union(cluster1_hub_genes_control, cluster1_hub_genes_30)
)

cluster1_hub_genes_control_unique <- setdiff(
  cluster1_hub_genes_control,
  union(cluster1_hub_genes_30, cluster1_hub_genes_35)
)

cluster2_hub_genes_30_unique <- setdiff(
  cluster2_hub_genes_30,
  union(cluster2_hub_genes_control, cluster2_hub_genes_35)
)

cluster2_hub_genes_35_unique <- setdiff(
  cluster2_hub_genes_35,
  union(cluster2_hub_genes_control, cluster2_hub_genes_30)
)

cluster2_hub_genes_control_unique <- setdiff(
  cluster2_hub_genes_control,
  union(cluster2_hub_genes_30, cluster2_hub_genes_35)
)

cluster3_hub_genes_30_unique <- setdiff(
  cluster3_hub_genes_30,
  union(cluster3_hub_genes_control, cluster3_hub_genes_35)
)
cluster3_hub_genes_35_unique <- setdiff(
  cluster3_hub_genes_35,
  union(cluster3_hub_genes_control, cluster3_hub_genes_30)
)
cluster3_hub_genes_control_unique <- setdiff(
  cluster3_hub_genes_control,
  union(cluster3_hub_genes_30, cluster3_hub_genes_35)
)

cluster9_hub_genes_30_unique <- setdiff(
  cluster9_hub_genes_30,
  union(cluster9_hub_genes_control, cluster9_hub_genes_35)
)

cluster9_hub_genes_35_unique <- setdiff(
  cluster9_hub_genes_35,
  union(cluster9_hub_genes_control, cluster9_hub_genes_30)
)

cluster9_hub_genes_control_unique <- setdiff(
  cluster9_hub_genes_control,
  union(cluster9_hub_genes_30, cluster9_hub_genes_35)
)

all_hub_genes_control <- unique(c(cluster1_hub_genes_control_unique, cluster2_hub_genes_control_unique, cluster3_hub_genes_control_unique, cluster9_hub_genes_control_unique))
all_hub_genes_30 <- unique(c(cluster1_hub_genes_30_unique, cluster2_hub_genes_30_unique, cluster3_hub_genes_30_unique, cluster9_hub_genes_30_unique))
all_hub_genes_35 <- unique(c(cluster1_hub_genes_35_unique, cluster2_hub_genes_35_unique, cluster3_hub_genes_35_unique, cluster9_hub_genes_35_unique))


# 3. OVERLAP ANALYSIS: UNIQUE HUB GENES vs SWITCHING GENES (WGCNA-filtered)
overlap_results_unique_wgcna <- data.frame()

# Control-specific hub genes
if(length(all_hub_genes_control) > 0) {
    overlap_control_all <- intersect(all_hub_genes_control, all_switching_genes_wgcna)
    overlap_control_vs_t30 <- intersect(all_hub_genes_control, switching_genes_control_vs_t30)
    overlap_control_vs_t35 <- intersect(all_hub_genes_control, switching_genes_control_vs_t35)
    overlap_t30_vs_t35_ctrl <- intersect(all_hub_genes_control, switching_genes_t30_vs_t35)
    
    overlap_results_unique_wgcna <- rbind(overlap_results_unique_wgcna, data.frame(
        temperature_specific = "Control",
        total_unique_hub_genes = length(all_hub_genes_control),
        overlap_any_switching = length(overlap_control_all),
        overlap_control_vs_t30 = length(overlap_control_vs_t30),
        overlap_control_vs_t35 = length(overlap_control_vs_t35),
        overlap_t30_vs_t35 = length(overlap_t30_vs_t35_ctrl),
        pct_overlap_any = round((length(overlap_control_all) / length(all_hub_genes_control)) * 100, 1),
        pct_overlap_control_t30 = round((length(overlap_control_vs_t30) / length(all_hub_genes_control)) * 100, 1),
        pct_overlap_control_t35 = round((length(overlap_control_vs_t35) / length(all_hub_genes_control)) * 100, 1),
        pct_overlap_t30_t35 = round((length(overlap_t30_vs_t35_ctrl) / length(all_hub_genes_control)) * 100, 1),
        stringsAsFactors = FALSE
    ))
}

# 30°C-specific hub genes
if(length(all_hub_genes_30) > 0) {
    overlap_30_all <- intersect(all_hub_genes_30, all_switching_genes_wgcna)
    overlap_30_vs_control_t30 <- intersect(all_hub_genes_30, switching_genes_control_vs_t30)
    overlap_30_vs_control_t35 <- intersect(all_hub_genes_30, switching_genes_control_vs_t35)
    overlap_30_vs_t30_t35 <- intersect(all_hub_genes_30, switching_genes_t30_vs_t35)
    
    overlap_results_unique_wgcna <- rbind(overlap_results_unique_wgcna, data.frame(
        temperature_specific = "30C",
        total_unique_hub_genes = length(all_hub_genes_30),
        overlap_any_switching = length(overlap_30_all),
        overlap_control_vs_t30 = length(overlap_30_vs_control_t30),
        overlap_control_vs_t35 = length(overlap_30_vs_control_t35),
        overlap_t30_vs_t35 = length(overlap_30_vs_t30_t35),
        pct_overlap_any = round((length(overlap_30_all) / length(all_hub_genes_30)) * 100, 1),
        pct_overlap_control_t30 = round((length(overlap_30_vs_control_t30) / length(all_hub_genes_30)) * 100, 1),
        pct_overlap_control_t35 = round((length(overlap_30_vs_control_t35) / length(all_hub_genes_30)) * 100, 1),
        pct_overlap_t30_t35 = round((length(overlap_30_vs_t30_t35) / length(all_hub_genes_30)) * 100, 1),
        stringsAsFactors = FALSE
    ))
}

# 35°C-specific hub genes
if(length(all_hub_genes_35) > 0) {
    overlap_35_all <- intersect(all_hub_genes_35, all_switching_genes_wgcna)
    overlap_35_vs_control_t30 <- intersect(all_hub_genes_35, switching_genes_control_vs_t30)
    overlap_35_vs_control_t35 <- intersect(all_hub_genes_35, switching_genes_control_vs_t35)
    overlap_35_vs_t30_t35 <- intersect(all_hub_genes_35, switching_genes_t30_vs_t35)
    
    overlap_results_unique_wgcna <- rbind(overlap_results_unique_wgcna, data.frame(
        temperature_specific = "35C",
        total_unique_hub_genes = length(all_hub_genes_35),
        overlap_any_switching = length(overlap_35_all),
        overlap_control_vs_t30 = length(overlap_35_vs_control_t30),
        overlap_control_vs_t35 = length(overlap_35_vs_control_t35),
        overlap_t30_vs_t35 = length(overlap_35_vs_t30_t35),
        pct_overlap_any = round((length(overlap_35_all) / length(all_hub_genes_35)) * 100, 1),
        pct_overlap_control_t30 = round((length(overlap_35_vs_control_t30) / length(all_hub_genes_35)) * 100, 1),
        pct_overlap_control_t35 = round((length(overlap_35_vs_control_t35) / length(all_hub_genes_35)) * 100, 1),
        pct_overlap_t30_t35 = round((length(overlap_35_vs_t30_t35) / length(all_hub_genes_35)) * 100, 1),
        stringsAsFactors = FALSE
    ))
}

cat("\nOverlap results (WGCNA-filtered switching genes):\n")
print(overlap_results_unique_wgcna)


# 4. DETAILED ANALYSIS OF UNIQUE HUB-SWITCHING GENES
# Store the actual overlapping genes for each temperature
unique_hub_switching_genes <- list()

if(length(all_hub_genes_control) > 0) {
    control_overlap <- intersect(all_hub_genes_control, all_switching_genes_wgcna)  
    unique_hub_switching_genes[["control_specific"]] <- control_overlap
    cat("Control-specific hub genes with switching:", length(control_overlap), "\n")
    if(length(control_overlap) > 0) {
        cat("Examples:", paste(head(control_overlap, 3), collapse = ", "), "\n")
    }
}

if(length(all_hub_genes_30) > 0) {
    t30_overlap <- intersect(all_hub_genes_30, all_switching_genes_wgcna)  
    unique_hub_switching_genes[["30C_specific"]] <- t30_overlap
    cat("30°C-specific hub genes with switching:", length(t30_overlap), "\n")
    if(length(t30_overlap) > 0) {
        cat("Examples:", paste(head(t30_overlap, 3), collapse = ", "), "\n")
    }
}

if(length(all_hub_genes_35) > 0) {
    t35_overlap <- intersect(all_hub_genes_35, all_switching_genes_wgcna)  
    unique_hub_switching_genes[["35C_specific"]] <- t35_overlap
    cat("35°C-specific hub genes with switching:", length(t35_overlap), "\n")
    if(length(t35_overlap) > 0) {
        cat("Examples:", paste(head(t35_overlap, 3), collapse = ", "), "\n")
    }
}

# 5. VISUALIZATIONS

if(nrow(overlap_results_unique_wgcna) > 0) {  
    library(ggplot2)
    library(reshape2)
    
    # Plot 1: Percentage of unique hub genes with switching 
    p1_wgcna <- ggplot(overlap_results_unique_wgcna, aes(x = temperature_specific, y = pct_overlap_any)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_text(aes(label = paste0(pct_overlap_any, "%\n(", overlap_any_switching, "/", total_unique_hub_genes, ")")), 
                  vjust = -0.5) +
        labs(title = "Percentage of Unique Hub Genes with Isoform Switching",
             subtitle = "Temperature-specific hub genes that undergo alternative splicing (WGCNA-filtered)",
             x = "Temperature-Specific Hub Genes",
             y = "% with Isoform Switching") +
        theme_bw() +
        ylim(0, max(overlap_results_unique_wgcna$pct_overlap_any) * 1.2)
    
    print(p1_wgcna)
    
    # Plot 2: Detailed comparison breakdown
    overlap_matrix_wgcna <- overlap_results_unique_wgcna %>%
        select(temperature_specific, pct_overlap_control_t30, pct_overlap_control_t35, pct_overlap_t30_t35) %>%
        reshape2::melt(id.vars = "temperature_specific") %>%
        mutate(comparison = case_when(
            variable == "pct_overlap_control_t30" ~ "Control vs 30°C",
            variable == "pct_overlap_control_t35" ~ "Control vs 35°C",
            variable == "pct_overlap_t30_t35" ~ "30°C vs 35°C"
        ))
    
    p2_wgcna <- ggplot(overlap_matrix_wgcna, aes(x = temperature_specific, y = comparison, fill = value)) +
        geom_tile(color = "white") +
        geom_text(aes(label = paste0(value, "%")), color = "black") +
        scale_fill_gradient(low = "white", high = "darkred", name = "% Overlap") +
        labs(title = "Unique Hub Gene - Switching Gene Overlap by Comparison",
             subtitle = "Which switching comparisons involve temperature-specific hub genes (WGCNA-filtered)",
             x = "Temperature-Specific Hub Genes",
             y = "Switching Comparison") +
        theme_bw()
    
    print(p2_wgcna)
    
    # # Plot 3: Absolute numbers (WGCNA-filtered)
    # p3_wgcna <- ggplot(overlap_results_unique_wgcna, aes(x = temperature_specific, y = overlap_any_switching)) +
    #     geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
    #     geom_text(aes(label = overlap_any_switching), vjust = -0.5) +
    #     labs(title = "Number of Unique Hub Genes with Isoform Switching",
    #          subtitle = "Absolute count of temperature-specific hub genes undergoing splicing (WGCNA-filtered)",
    #          x = "Temperature-Specific Hub Genes",
    #          y = "Number of Genes") +
    #     theme_bw() +
    #     ylim(0, max(overlap_results_unique_wgcna$overlap_any_switching) * 1.2)
    
    # print(p3_wgcna)
}

# 6. SAVE RESULTS

# Save overlap summary
write.csv(overlap_results_unique, "unique_hub_switching_overlap_summary_Pcom.csv", row.names = FALSE)


# 7. SLOPE GRAPH - Shows how overlaps change across comparisons
p7_temp_ordered <- overlap_matrix_wgcna %>%  # Changed from overlap_matrix to overlap_matrix_wgcna
    mutate(temp_order = case_when(
        temperature_specific == "Control" ~ 1,
        temperature_specific == "30C" ~ 2,
        temperature_specific == "35C" ~ 3
    )) %>%
    ggplot(aes(x = temp_order, y = value, color = comparison)) +
    geom_line(aes(group = comparison), size = 1.5, alpha = 0.8) +
    geom_point(size = 4) +
    geom_text(aes(label = paste0(value, "%")), 
              nudge_y = 2, size = 3, fontface = "bold") +
    scale_color_manual(values = c("Control vs 30°C" = "#d9e31aff", 
                                 "Control vs 35°C" = "#ff0d00ff", 
                                 "30°C vs 35°C" = "#ef750aff"),
                      name = "Switching\nComparison") +
    scale_x_continuous(breaks = 1:3, 
                      labels = c("Control\nHub Genes", "30°C\nHub Genes", "35°C\nHub Genes"),
                      expand = c(0.1, 0.1)) +
    labs(title = "Temperature-Specific Hub Gene Involvement in Switching",
         subtitle = "Percentage of each temperature's unique hub genes that undergo switching (WGCNA-filtered)",  # Updated subtitle
         x = "Temperature-Specific Hub Genes",
         y = "% Overlap with Switching Genes",
         color = "Switching\nComparison") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11, hjust = 0.5),
          legend.position = "bottom",
          panel.grid.minor.x = element_blank())

print(p7_temp_ordered) #save as landscape 6x4 	switch_isof_hub_genes_Pcom.pdf




# STATISTICAL SIGNIFICANCE TESTING FOR HUB-SWITCHING OVERLAPS
## test only the annotated genes (do so ones present in the genome found with WGCNA)
library(stats)

# Example Fisher table - switching in Control vs 30°C and also hub genes at control
# A: genes switching in control vs 30 and are control hub genes
# B: genes not switching in control vs 30 (total number of genes with isoforms minus those switching at control) and are control hub genes (all genes in DetExpr from WGCNA minus control hubs)
# C: genes switching in control vs 30 and are NOT control hub genes
# D: genes not switching in control vs 30 and are NOT control hub genes

# Extract all switching genes from the analysis
all_switching_genes_raw <- extractTopSwitches(
    switchListAnalyzed,
    filterForConsequences = FALSE,  # Get all switching genes, not just those with consequences
    n = Inf,                       # Get all significant switches
    sortByQvals = FALSE,           # Don't sort
    extractGenes = TRUE            # Extract gene-level results
)

# Extract switching genes for each specific comparison - BEFORE WGCNA filtering
switching_genes_control_vs_30_raw <- unique(all_switching_genes_raw[
    paste(all_switching_genes_raw$condition_1, "vs", all_switching_genes_raw$condition_2) == "control vs t30", 
    "gene_id"
])

switching_genes_control_vs_35_raw <- unique(all_switching_genes_raw[
    paste(all_switching_genes_raw$condition_1, "vs", all_switching_genes_raw$condition_2) == "control vs t35", 
    "gene_id"
])

switching_genes_30_vs_35_raw <- unique(all_switching_genes_raw[
    paste(all_switching_genes_raw$condition_1, "vs", all_switching_genes_raw$condition_2) == "t30 vs t35", 
    "gene_id"
])


# Get all genes used in WGCNA analysis
all_WGCNA_genes <- colnames(datExpr) # all genes in WGCNA expression matrix

# FILTER switching genes to only include those in WGCNA analysis
switching_genes_control_vs_30 <- intersect(switching_genes_control_vs_30_raw, all_WGCNA_genes)
switching_genes_control_vs_35 <- intersect(switching_genes_control_vs_35_raw, all_WGCNA_genes)
switching_genes_30_vs_35 <- intersect(switching_genes_30_vs_35_raw, all_WGCNA_genes)

# Print filtering results
cat("Control vs 30°C: ", length(switching_genes_control_vs_30_raw), " -> ", length(switching_genes_control_vs_30), "\n")
cat("Control vs 35°C: ", length(switching_genes_control_vs_35_raw), " -> ", length(switching_genes_control_vs_35), "\n")
cat("30°C vs 35°C: ", length(switching_genes_30_vs_35_raw), " -> ", length(switching_genes_30_vs_35), "\n")


# Get all genes with isoforms to determine non-switching genes - USE WGCNA GENES
# Use only genes that are both in the isoform analysis AND in WGCNA
all_genes_in_isoform_analysis <- unique(switchListAnalyzed$isoformFeatures$gene_id)
all_genes_in_analysis <- intersect(all_genes_in_isoform_analysis, all_WGCNA_genes)

cat("Total genes with isoforms: ", length(all_genes_in_analysis), "\n")
#Total genes with isoforms:  861 

# Extract non-switching genes for each comparison
not_switching_genes_control_vs_30 <- setdiff(all_genes_in_analysis, switching_genes_control_vs_30)
not_switching_genes_control_vs_35 <- setdiff(all_genes_in_analysis, switching_genes_control_vs_35)
not_switching_genes_30_vs_35 <- setdiff(all_genes_in_analysis, switching_genes_30_vs_35)

# Set hub genes variables (you'll need to define these based on the hub gene analysis)
hub_genes_control <- all_hub_genes_control
hub_genes_30 <- all_hub_genes_30
hub_genes_35 <- all_hub_genes_35

# Total number of genes
ncol(datExpr) #WGCNA expression matrix
#> [1] 28378

all_WGCNA_genes <- colnames(datExpr) # all genes in WGCNA expression matrix

# Create non-hub gene sets (from WGCNA genes)
not_hub_genes_control <- setdiff(all_WGCNA_genes, hub_genes_control)
not_hub_genes_30 <- setdiff(all_WGCNA_genes, hub_genes_30)
not_hub_genes_35 <- setdiff(all_WGCNA_genes, hub_genes_35)


#### Fisher's exact tests for all three comparisons
#testing whether switching genes (with consequences) are enriched among hub genes, the background is the total WGCNA genes

# For each comparison, create non-switching genes from the FULL WGCNA background
not_switching_genes_control_vs_30_corrected <- setdiff(all_WGCNA_genes, switching_genes_control_vs_30)
not_switching_genes_control_vs_35_corrected <- setdiff(all_WGCNA_genes, switching_genes_control_vs_35)
not_switching_genes_30_vs_35_corrected <- setdiff(all_WGCNA_genes, switching_genes_30_vs_35)

switching_genes_control_vs_30 <- switching_genes_control_vs_t30  # Use existing variable
switching_genes_control_vs_35 <- switching_genes_control_vs_t35  # Use existing variable
switching_genes_30_vs_35 <- switching_genes_t30_vs_t35          # Use existing variable

# 1. Control vs 30°C comparison with control hub genes
A1_corrected <- length(intersect(switching_genes_control_vs_t30, hub_genes_control))  
B1_corrected <- length(intersect(not_switching_genes_control_vs_30_corrected, hub_genes_control))
C1_corrected <- length(intersect(switching_genes_control_vs_t30, not_hub_genes_control))  
D1_corrected <- length(intersect(not_switching_genes_control_vs_30_corrected, not_hub_genes_control))

#                     | Hub Gene | Not Hub Gene |
# Switching genes     |    A1    |      C1      |
# Not switching genes |    B1    |      D1      |

# Verify the totals add up correctly
total_check_1 <- A1_corrected + B1_corrected + C1_corrected + D1_corrected
cat("Control vs 30°C - Total genes in Fisher test:", total_check_1, "Should equal:", length(all_WGCNA_genes), "\n")

fisher_matrix_1_corrected <- matrix(c(A1_corrected, B1_corrected, C1_corrected, D1_corrected), nrow = 2, byrow = TRUE,
                        dimnames = list(
                          Switching = c("Switching", "Not Switching"),
                          Hub = c("Control Hub", "Not Control Hub")
                        ))

print(fisher_matrix_1_corrected)
fisher_result_1_corrected <- fisher.test(fisher_matrix_1_corrected)
print(fisher_result_1_corrected)
# Fisher's Exact Test for Count Data
# data:  fisher_matrix_1_corrected
# p-value = 0.4642
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.7516671 1.6914220
# sample estimates:
# odds ratio 
#   1.148815 



# 2. Control vs 30°C comparison with 30°C hub genes
A2_corrected <- length(intersect(switching_genes_control_vs_t30, hub_genes_30))  # ← Use t30
B2_corrected <- length(intersect(not_switching_genes_control_vs_30_corrected, hub_genes_30))
C2_corrected <- length(intersect(switching_genes_control_vs_t30, not_hub_genes_30))  # ← Use t30
D2_corrected <- length(intersect(not_switching_genes_control_vs_30_corrected, not_hub_genes_30))

fisher_matrix_2_corrected <- matrix(c(A2_corrected, B2_corrected, C2_corrected, D2_corrected), nrow = 2, byrow = TRUE,
                        dimnames = list(
                          Switching = c("Switching", "Not Switching"),
                          Hub = c("30°C Hub", "Not 30°C Hub")
                        ))

cat("\nCORRECTED Fisher's test for Control vs 30°C with 30°C hub genes:\n")
print(fisher_matrix_2_corrected)
fisher_result_2_corrected <- fisher.test(fisher_matrix_2_corrected)
print(fisher_result_2_corrected)
# 	Fisher's Exact Test for Count Data
# data:  fisher_matrix_2_corrected
# p-value = 0.6046
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.5367687 1.3364337
# sample estimates:
# odds ratio 
#  0.8681611 



# 3. Control vs 35°C comparison with 35°C hub genes  
A3_corrected <- length(intersect(switching_genes_control_vs_t35, hub_genes_35))  # ← Use t35
B3_corrected <- length(intersect(not_switching_genes_control_vs_35_corrected, hub_genes_35))
C3_corrected <- length(intersect(switching_genes_control_vs_t35, not_hub_genes_35))  # ← Use t35
D3_corrected <- length(intersect(not_switching_genes_control_vs_35_corrected, not_hub_genes_35))

fisher_matrix_3_corrected <- matrix(c(A3_corrected, B3_corrected, C3_corrected, D3_corrected), nrow = 2, byrow = TRUE,
                        dimnames = list(
                          Switching = c("Switching", "Not Switching"),
                          Hub = c("35°C Hub", "Not 35°C Hub")
                        ))

cat("\nCORRECTED Fisher's test for Control vs 35°C with 35°C hub genes:\n")
print(fisher_matrix_3_corrected)
fisher_result_3_corrected <- fisher.test(fisher_matrix_3_corrected)
print(fisher_result_3_corrected)
# data:  fisher_matrix_3_corrected
# p-value = 0.2782
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.581320 3.795876
# sample estimates:
# odds ratio 
#   1.646378


# 4. Control vs 35°C comparison with control hub genes 
A4_corrected <- length(intersect(switching_genes_control_vs_t35, hub_genes_control))  # ← Use t35
B4_corrected <- length(intersect(not_switching_genes_control_vs_35_corrected, hub_genes_control))
C4_corrected <- length(intersect(switching_genes_control_vs_t35, not_hub_genes_control))  # ← Use t35
D4_corrected <- length(intersect(not_switching_genes_control_vs_35_corrected, not_hub_genes_control))

fisher_matrix_4_corrected <- matrix(c(A4_corrected, B4_corrected, C4_corrected, D4_corrected), nrow = 2, byrow = TRUE,
                        dimnames = list(
                          Switching = c("Switching", "Not Switching"),
                          Hub = c("Control Hub", "Not Control Hub")
                        ))

cat("\nCORRECTED Fisher's test for Control vs 35°C with Control hub genes:\n")
print(fisher_matrix_4_corrected)
fisher_result_4_corrected <- fisher.test(fisher_matrix_4_corrected)
print(fisher_result_4_corrected)
# data:  fisher_matrix_4_corrected
# p-value = 0.2746
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.5914035 3.8621018
# sample estimates:
# odds ratio 


# 5. 30°C vs 35°C comparison with 30°C hub genes
A5_corrected <- length(intersect(switching_genes_t30_vs_t35, hub_genes_30))  # ← Use t30_vs_t35
B5_corrected <- length(intersect(not_switching_genes_30_vs_35_corrected, hub_genes_30))
C5_corrected <- length(intersect(switching_genes_t30_vs_t35, not_hub_genes_30))  # ← Use t30_vs_t35
D5_corrected <- length(intersect(not_switching_genes_30_vs_35_corrected, not_hub_genes_30))

fisher_matrix_5_corrected <- matrix(c(A5_corrected, B5_corrected, C5_corrected, D5_corrected), nrow = 2, byrow = TRUE,
                        dimnames = list(
                          Switching = c("Switching", "Not Switching"),
                          Hub = c("30°C Hub", "Not 30°C Hub")
                        ))

cat("\nCORRECTED Fisher's test for 30°C vs 35°C with 30°C hub genes:\n")
print(fisher_matrix_5_corrected)
fisher_result_5_corrected <- fisher.test(fisher_matrix_5_corrected)
print(fisher_result_5_corrected)
# data:  fisher_matrix_5_corrected
# p-value = 0.02988
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.018788 2.314729
# sample estimates:
# odds ratio 
#   1.563871 


# 6. 30°C vs 35°C comparison with 35°C hub genes
A6_corrected <- length(intersect(switching_genes_t30_vs_t35, hub_genes_35))  # ← Use t30_vs_t35
B6_corrected <- length(intersect(not_switching_genes_30_vs_35_corrected, hub_genes_35))
C6_corrected <- length(intersect(switching_genes_t30_vs_t35, not_hub_genes_35))  # ← Use t30_vs_t35
D6_corrected <- length(intersect(not_switching_genes_30_vs_35_corrected, not_hub_genes_35))

fisher_matrix_6_corrected <- matrix(c(A6_corrected, B6_corrected, C6_corrected, D6_corrected), nrow = 2, byrow = TRUE,
                        dimnames = list(
                          Switching = c("Switching", "Not Switching"),
                          Hub = c("35°C Hub", "Not 35°C Hub")
                        ))

cat("\nCORRECTED Fisher's test for 30°C vs 35°C with 35°C hub genes:\n")
print(fisher_matrix_6_corrected)
fisher_result_6_corrected <- fisher.test(fisher_matrix_6_corrected)
print(fisher_result_6_corrected)
# data:  fisher_matrix_6_corrected
# p-value = 0.1179
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.8796596 2.0911365
# sample estimates:
# odds ratio 
#   1.384856 



###### COMBINED COMPARISON ANALYSIS: (30°C vs 35°C) + (Control vs 35°C)

# 1. CREATE COMBINED SWITCHING GENE LIST
# Union of switching genes from both comparisons
combined_switching_genes <- union(switching_genes_t30_vs_t35, switching_genes_control_vs_t35)

# 2. CREATE COMBINED NON-SWITCHING GENE LIST  
# Genes that DON'T switch in EITHER comparison
combined_not_switching_genes <- setdiff(all_WGCNA_genes, combined_switching_genes)

# 3. FISHER'S EXACT TESTS FOR EACH HUB GENE SET

# A. Combined comparison with Control hub genes

A_combined_control <- length(intersect(combined_switching_genes, hub_genes_control))
B_combined_control <- length(intersect(combined_not_switching_genes, hub_genes_control))
C_combined_control <- length(intersect(combined_switching_genes, not_hub_genes_control))
D_combined_control <- length(intersect(combined_not_switching_genes, not_hub_genes_control))

fisher_matrix_combined_control <- matrix(c(A_combined_control, B_combined_control, C_combined_control, D_combined_control), 
                                        nrow = 2, byrow = TRUE,
                                        dimnames = list(
                                          Switching = c("Switching", "Not Switching"),
                                          Hub = c("Control Hub", "Not Control Hub")
                                        ))

cat("Contingency table:\n")
print(fisher_matrix_combined_control)
fisher_result_combined_control <- fisher.test(fisher_matrix_combined_control)
print(fisher_result_combined_control)
# data:  fisher_matrix_combined_control
# p-value = 0.03996
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9920534 2.1860305
# sample estimates:
# odds ratio 
#   1.497848 


# Calculate percentages for interpretation
control_hub_switch_pct <- round((A_combined_control / (A_combined_control + B_combined_control)) * 100, 1)
non_control_hub_switch_pct <- round((C_combined_control / (C_combined_control + D_combined_control)) * 100, 1)
cat(sprintf("Control hub genes switching: %d/%d (%s%%)\n", 
           A_combined_control, A_combined_control + B_combined_control, control_hub_switch_pct))
           #Control hub genes switching: 30/1529 (2%)
cat(sprintf("Non-control hub genes switching: %d/%d (%s%%)\n", 
           C_combined_control, C_combined_control + D_combined_control, non_control_hub_switch_pct))
           #Non-control hub genes switching: 354/26849 (1.3%)




# B. Combined comparison with 30°C hub genes
A_combined_30 <- length(intersect(combined_switching_genes, hub_genes_30))
B_combined_30 <- length(intersect(combined_not_switching_genes, hub_genes_30))
C_combined_30 <- length(intersect(combined_switching_genes, not_hub_genes_30))
D_combined_30 <- length(intersect(combined_not_switching_genes, not_hub_genes_30))

fisher_matrix_combined_30 <- matrix(c(A_combined_30, B_combined_30, C_combined_30, D_combined_30), 
                                   nrow = 2, byrow = TRUE,
                                   dimnames = list(
                                     Switching = c("Switching", "Not Switching"),
                                     Hub = c("30°C Hub", "Not 30°C Hub")
                                   ))

cat("Contingency table:\n")
print(fisher_matrix_combined_30)
fisher_result_combined_30 <- fisher.test(fisher_matrix_combined_30)
print(fisher_result_combined_30)
# data:  fisher_matrix_combined_30
# p-value = 0.03152
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.014245 2.207744
# sample estimates:
# odds ratio 
#   1.521034 

# Calculate percentages for interpretation
hub_30_switch_pct <- round((A_combined_30 / (A_combined_30 + B_combined_30)) * 100, 1)
non_hub_30_switch_pct <- round((C_combined_30 / (C_combined_30 + D_combined_30)) * 100, 1)
cat(sprintf("30°C hub genes switching: %d/%d (%s%%)\n", 
           A_combined_30, A_combined_30 + B_combined_30, hub_30_switch_pct))
           #30°C hub genes switching: 31/1559 (2%)
cat(sprintf("Non-30°C hub genes switching: %d/%d (%s%%)\n", 
           C_combined_30, C_combined_30 + D_combined_30, non_hub_30_switch_pct))
           #Non-30°C hub genes switching: 353/26819 (1.3%)



# C. Combined comparison with 35°C hub genes
A_combined_35 <- length(intersect(combined_switching_genes, hub_genes_35))
B_combined_35 <- length(intersect(combined_not_switching_genes, hub_genes_35))
C_combined_35 <- length(intersect(combined_switching_genes, not_hub_genes_35))
D_combined_35 <- length(intersect(combined_not_switching_genes, not_hub_genes_35))

fisher_matrix_combined_35 <- matrix(c(A_combined_35, B_combined_35, C_combined_35, D_combined_35), 
                                   nrow = 2, byrow = TRUE,
                                   dimnames = list(
                                     Switching = c("Switching", "Not Switching"),
                                     Hub = c("35°C Hub", "Not 35°C Hub")
                                   ))

cat("Contingency table:\n")
print(fisher_matrix_combined_35)
fisher_result_combined_35 <- fisher.test(fisher_matrix_combined_35)
print(fisher_result_combined_35)
# data:  fisher_matrix_combined_35
# p-value = 0.1141
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.8905049 2.0140466
# sample estimates:
# odds ratio 
#   1.364162 



# Calculate percentages for interpretation
hub_35_switch_pct <- round((A_combined_35 / (A_combined_35 + B_combined_35)) * 100, 1)
non_hub_35_switch_pct <- round((C_combined_35 / (C_combined_35 + D_combined_35)) * 100, 1)
cat(sprintf("35°C hub genes switching: %d/%d (%s%%)\n", 
           A_combined_35, A_combined_35 + B_combined_35, hub_35_switch_pct))
           #35°C hub genes switching: 28/1554 (1.8%)
cat(sprintf("Non-35°C hub genes switching: %d/%d (%s%%)\n", 
           C_combined_35, C_combined_35 + D_combined_35, non_hub_35_switch_pct))
           #Non-35°C hub genes switching: 356/26824 (1.3%)
           )

# 4. SUMMARY OF ALL COMBINED RESULTS
results_summary <- data.frame(
  Hub_Gene_Type = c("Control Hub", "30°C Hub", "35°C Hub"),
  P_Value = c(fisher_result_combined_control$p.value, 
              fisher_result_combined_30$p.value, 
              fisher_result_combined_35$p.value),
  Odds_Ratio = c(fisher_result_combined_control$estimate, 
                 fisher_result_combined_30$estimate, 
                 fisher_result_combined_35$estimate),
  Hub_Switching_Pct = c(control_hub_switch_pct, hub_30_switch_pct, hub_35_switch_pct),
  NonHub_Switching_Pct = c(non_control_hub_switch_pct, non_hub_30_switch_pct, non_hub_35_switch_pct),
  Significant = c(fisher_result_combined_control$p.value < 0.05,
                  fisher_result_combined_30$p.value < 0.05,
                  fisher_result_combined_35$p.value < 0.05)
)

print(results_summary)

# DETAILED BREAKDOWN OF COMBINED SWITCHING GENES
cat("Genes switching in 30°C vs 35°C only:", length(setdiff(switching_genes_t30_vs_t35, switching_genes_control_vs_t35)), "\n")
cat("Genes switching in Control vs 35°C only:", length(setdiff(switching_genes_control_vs_t35, switching_genes_t30_vs_t35)), "\n")
cat("Genes switching in BOTH comparisons:", length(intersect(switching_genes_t30_vs_t35, switching_genes_control_vs_t35)), "\n")
cat("Total combined switching genes:", length(combined_switching_genes), "\n")

# Save combined results
write.csv(results_summary, "/scratch3/workspace/federica_scucchia_uri_edu-altSplice/20250424_ENCORE_HawaiiTPC_Federica/output/Rstudio/Isoform_switch/Pcom/combined_35C_fisher_results.csv", row.names = FALSE)