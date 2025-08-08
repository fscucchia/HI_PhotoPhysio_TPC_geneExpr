

########## Get Gene Lengths

# Use the GenomicFeatures R package to extract gene lengths from a GTF file

## Install if needed
# BiocManager::install("GenomicFeatures")
# BiocManager::install("txdbmaker")
library(GenomicFeatures)
library(txdbmaker)

# Build a TxDb object from GTF files - run each species separately
#txdb <- makeTxDbFromGFF("Pocillopora_acuta_HIv2_modified.gtf", format="gtf")
txdb <- makeTxDbFromGFF("stringtie_merged_noIso_Pacu.gtf", format="gtf")
txdb <- makeTxDbFromGFF("stringtie_merged_noIso_Mcap.gtf", format="gtf")
txdb <- makeTxDbFromGFF("stringtie_merged_noIso_Pcom.gtf", format="gtf")

# Get exonic lengths per gene - Get gene lengths as sum of non-overlapping exons per gene
exons.list.per.gene <- exonsBy(txdb, by="gene")

# Reduce overlapping exons
#reduce() ensures that overlapping exons are merged, so each base is counted only once.
reduced_exons <- reduce(exons.list.per.gene)

# Sum exon lengths
gene_lengths <- sum(width(reduced_exons)) 
# gene.lengths is a named vector: names are gene IDs, values are lengths

# Convert to data frame
gene_lengths_df <- data.frame(
    gene_id = names(gene_lengths),
    length = as.integer(gene_lengths)
)

write.table(gene_lengths_df, file = "gene_lengths_Pacu.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_lengths_df, file = "gene_lengths_Mcap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_lengths_df, file = "gene_lengths_Pcom.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


