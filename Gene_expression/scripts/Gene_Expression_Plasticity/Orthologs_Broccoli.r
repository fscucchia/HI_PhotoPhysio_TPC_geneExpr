
############ Orthologs search #############

## load libraries
library(ggplot2)
library(ggtree)
library(ape)
library(ggtreeExtra)

## Load raw count matrix
counts_mcap <- read.csv("gene_count_matrix_noIso_Mcap2.csv", row.names = 1)
counts_pacu <- read.csv("gene_count_matrix_noIso_Pacu.csv", row.names = 1)
counts_pcom <- read.csv("gene_count_matrix_noIso_Pcom.csv", row.names = 1)

## Load gene lengths
gene_length_mcap <- read.delim("gene_lengths_Mcap.tsv", row.names = 1)
gene_length_pacu <- read.delim("gene_lengths_Pacu.tsv", row.names = 1) 
gene_length_pcom <- read.delim("gene_lengths_Pcom.tsv", row.names = 1)

#For each species individually:
#Compute TPMs using its raw count matrix and gene lengths.
#Result: TPM matrix per species (genes × samples).
#TPM normalization already includes:
#•Gene length correction
#•Sequencing depth correction (through the per-sample scaling to 1 million)

# Make sure gene_length is named with gene IDs
# Reorder gene_length to match the rownames of counts

gene_length_vec_mcap <- gene_length_mcap$length
names(gene_length_vec_mcap) <- rownames(gene_length_mcap)
gene_length_vec_mcap <- gene_length_vec_mcap[rownames(counts_mcap)]
#remove NAs
gene_length_vec_mcap <- gene_length_vec_mcap[!is.na(gene_length_vec_mcap)]

gene_length_vec_pacu <- gene_length_pacu$length
names(gene_length_vec_pacu) <- rownames(gene_length_pacu)
gene_length_vec_pacu <- gene_length_vec_pacu[rownames(counts_pacu)]
#remove NAs
gene_length_vec_pacu <- gene_length_vec_pacu[!is.na(gene_length_vec_pacu)]

gene_length_vec_pcom <- gene_length_pcom$length
names(gene_length_vec_pcom) <- rownames(gene_length_pcom)
gene_length_vec_pcom <- gene_length_vec_pcom[rownames(counts_pcom)]
#remove NAs
gene_length_vec_pcom <- gene_length_vec_pcom[!is.na(gene_length_vec_pcom)]

counts_to_tpm <- function(counts, gene_length) {
  # gene_length should be in base pairs
  gene_length_kb <- gene_length / 1000
  rpk <- counts / gene_length_kb
  scaling_factors <- colSums(rpk)
  tpm <- t(t(rpk) / scaling_factors) * 1e6
  return(tpm)
}

tpm_mcap <- counts_to_tpm(counts_mcap, gene_length_vec_mcap)
tpm_pacu <- counts_to_tpm(counts_pacu, gene_length_vec_pacu)
tpm_pcom <- counts_to_tpm(counts_pcom, gene_length_vec_pcom)

#Because TPMs are continuous and right-skewed, apply a log transformation
logTPM_mcap <- log2(tpm_mcap + 1)
logTPM_pacu <- log2(tpm_pacu + 1)
logTPM_pcom <- log2(tpm_pcom + 1)

#Filter out lowly expressed genes. TPM includes a lot of noise at the low end (e.g., TPM < 1). Filtering improves statistical power and avoids modeling noise. A common filter:
keep <- rowMeans(tpm_mcap) > 1
logTPM_filtered_mcap <- logTPM_mcap[keep, ]

keep <- rowMeans(tpm_pacu) > 1
logTPM_filtered_pacu <- logTPM_pacu[keep, ]

keep <- rowMeans(tpm_pcom) > 1
logTPM_filtered_pcom <- logTPM_pcom[keep, ]

#Use normalizeBetweenArrays() from limma to further adjust for distribution differences.
library(limma)
# Normalize across samples (recommended for TPMs)
logTPM_norm_mcap <- normalizeBetweenArrays(logTPM_filtered_mcap, method = "quantile")
logTPM_norm_pacu <- normalizeBetweenArrays(logTPM_filtered_pacu, method = "quantile")
logTPM_norm_pcom <- normalizeBetweenArrays(logTPM_filtered_pcom, method = "quantile")


## Load orthologs
pairs <- read.table("orthologous_pairs.txt", stringsAsFactors = FALSE)
colnames(pairs) <- c("gene1", "gene2")

## Build ortholog groups (connected components)
# get 1:1:1 orthologs (one gene from each of the three species)

library(igraph)

# Build graph
g <- graph_from_data_frame(pairs, directed = FALSE)

# Find connected components (ortholog groups)
groups <- components(g)$membership

## Extract groups with exactly one gene from each species

#  get the species for each gene:
gene_species <- sub("^(.*?)___.*", "\\1", names(groups))
group_df <- data.frame(gene = names(groups), group = groups, species = gene_species)

# Count number of unique species per group
library(dplyr)
one_to_one_to_one <- group_df %>%
  group_by(group) %>%
  filter(n_distinct(species) == 3, n() == 3) %>%  # 3 species, 3 genes
  ungroup()


# Spread to wide format: one row per group, one column per species
library(tidyr)
ortholog_table <- one_to_one_to_one %>%
  dplyr::select(group, species, gene) %>%
  tidyr::pivot_wider(names_from = species, values_from = gene)
#9010 orthologous groups with one gene from each species


# For each species, get the gene IDs from ortholog_table
genes_mcap <- ortholog_table$Montipora_capitata_HIv3   
genes_pacu <- ortholog_table$Pocillopora_acuta_HIv2   
genes_pcom <- ortholog_table$Porites_compressa_HIv1 

# Create keep_rows BEFORE subsetting any gene vectors
keep_rows <- genes_mcap %in% rownames(logTPM_norm_mcap) &
             genes_pacu %in% rownames(logTPM_norm_pacu) &
             genes_pcom %in% rownames(logTPM_norm_pcom)

ortholog_table_filt <- ortholog_table[keep_rows, ]
genes_mcap <- ortholog_table_filt$Montipora_capitata_HIv3
genes_pacu <- ortholog_table_filt$Pocillopora_acuta_HIv2
genes_pcom <- ortholog_table_filt$Porites_compressa_HIv1

## merge the three logTPM tables by matching the order of genes to the ortholog_table and then using cbind() to combine them
# Ensure the rows are in the same order as in ortholog_table for each species
logTPM_norm_mcap_ortho <- logTPM_norm_mcap[genes_mcap, ]
logTPM_norm_pacu_ortho <- logTPM_norm_pacu[genes_pacu, ]
logTPM_norm_pcom_ortho <- logTPM_norm_pcom[genes_pcom, ]

#rename the columns by pasting the species prefix to each column name
colnames(logTPM_norm_mcap_ortho) <- paste0("mcap_", colnames(logTPM_norm_mcap_ortho))
colnames(logTPM_norm_pacu_ortho) <- paste0("pacu_", colnames(logTPM_norm_pacu_ortho))
colnames(logTPM_norm_pcom_ortho) <- paste0("pcom_", colnames(logTPM_norm_pcom_ortho))

# Combine the three matrices column-wise
logTPM_merged <- cbind(logTPM_norm_mcap_ortho, logTPM_norm_pacu_ortho, logTPM_norm_pcom_ortho)

# Assign arbitrary ortholog group names as rownames
rownames(logTPM_merged) <- paste0("ortho", seq_len(nrow(logTPM_merged)))
### 8092 final orthologous genes with 1-to-1-to-1 mapping across the three species

# keep a mapping of orthoID to gene IDs, create a mapping table:
ortho_map <- data.frame(
  orthoID = rownames(logTPM_merged),
  mcap_gene = genes_mcap,
  pacu_gene = genes_pacu,
  pcom_gene = genes_pcom
)

# save ortho_map as a CSV file
write.csv(ortho_map, "ortholog_mapping.csv", row.names = FALSE)