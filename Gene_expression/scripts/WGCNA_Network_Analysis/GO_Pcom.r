
####### GO Enrichment Analysis of WGCNA top 10% hub genes per each cluster - Pcom

setwd("...................")
load(".RData") #load WCGNA data

#Load libraries
library("ViSEAGO") 
library(topGO)
library(tidyverse)
library(GSEABase)               #BiocManager::install("GSEABase")
library(data.table)
library(ggplot2)
library(cowplot)                #install.packages("cowplot")
library(dplyr)
library(tidyr)
library(readr)

##treatment information
treatmentinfo <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/RNAseq_Pcom_data.csv", header = TRUE, sep = ";")

#keep only control, 30 and 35°C
treatmentinfo <- treatmentinfo[treatmentinfo$temp %in% c(26.8, 30, 35), ]

##load gene clusters

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

# Find shared genes among all three temperatures for each cluster
shared_genes_c1 <- Reduce(intersect, list(cluster1_hub_genes_control, cluster1_hub_genes_30, cluster1_hub_genes_35))
shared_genes_c2 <- Reduce(intersect, list(cluster2_hub_genes_control, cluster2_hub_genes_30, cluster2_hub_genes_35))
shared_genes_c3 <- Reduce(intersect, list(cluster3_hub_genes_control, cluster3_hub_genes_30, cluster3_hub_genes_35))
shared_genes_c9 <- Reduce(intersect, list(cluster9_hub_genes_control, cluster9_hub_genes_30, cluster9_hub_genes_35))

# Save as CSV
write.csv(data.frame(Gene = shared_genes_c1), "shared_genes_cluster1_Pcom.csv", row.names = FALSE)
write.csv(data.frame(Gene = shared_genes_c2), "shared_genes_cluster2_Pcom.csv", row.names = FALSE)
write.csv(data.frame(Gene = shared_genes_c3), "shared_genes_cluster3_Pcom.csv", row.names = FALSE)
write.csv(data.frame(Gene = shared_genes_c9), "shared_genes_cluster9_Pcom.csv", row.names = FALSE)

# ##Load Pcom annotations
# final_GO_table_Pcom <- read.csv("final_GO_table_Pcom.csv", header = TRUE, sep = ",")

# ### filter clusters data data for just expressed genes with GO annotations
# cluster1_hub_genes_shared_GO <- final_GO_table_Pcom %>%
#   filter(protein_accession %in% cluster1_hub_genes_shared) ##no GO

# cluster2_hub_genes_shared_GO <- final_GO_table_Pcom %>%
#   filter(protein_accession %in% cluster2_hub_genes_shared) ##no GO

# cluster3_hub_genes_shared_GO <- final_GO_table_Pcom %>%
#   filter(protein_accession %in% cluster3_hub_genes_shared)  ##no GO

# cluster9_hub_genes_shared_GO <- final_GO_table_Pcom %>%
#   filter(protein_accession %in% cluster9_hub_genes_shared)  ##no GO


#Load diamond results
blast <- read_tsv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Pcom/Pcom.diamondBlastpNCBInr", col_names = FALSE)
colnames(blast) <- c("Gene", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(blast)
dim(blast)
#[1] 40740    12

### find top_hit from blast output corresponding to the shared hub genes per each cluster
# Create a data frame from your gene list
shared_genes_cluster1 <- data.frame(Gene = shared_genes_c1, stringsAsFactors = FALSE)
shared_genes_cluster2 <- data.frame(Gene = shared_genes_c2, stringsAsFactors = FALSE)
shared_genes_cluster3 <- data.frame(Gene = shared_genes_c3, stringsAsFactors = FALSE)
shared_genes_cluster9 <- data.frame(Gene = shared_genes_c9, stringsAsFactors = FALSE)

# Inner join to get matching rows
shared_genes_blast_cluster1  <- inner_join(shared_genes_cluster1, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster2  <- inner_join(shared_genes_cluster2, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster3  <- inner_join(shared_genes_cluster3, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster9  <- inner_join(shared_genes_cluster9, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

#search proteins in NCBI
library(rentrez)

# Get the vector of accession numbers
accessions_cluster1 <- shared_genes_blast_cluster1$top_hit
accessions_cluster2 <- shared_genes_blast_cluster2$top_hit
accessions_cluster3 <- shared_genes_blast_cluster3$top_hit
accessions_cluster9 <- shared_genes_blast_cluster9$top_hit

# Function to fetch protein name and organism
fetch_protein_info <- function(acc) {
  rec <- tryCatch(
    entrez_fetch(db = "protein", id = acc, rettype = "gb", retmode = "text"),
    error = function(e) NA
  )
  if (is.na(rec)) return(data.frame(accessions = acc, protein_name = NA, organism = NA))
  protein_name <- sub("DEFINITION  (.*)", "\\1", regmatches(rec, regexpr("DEFINITION  .*", rec)))
  organism <- sub("  ORGANISM  (.*)", "\\1", regmatches(rec, regexpr("  ORGANISM  .*", rec)))
  data.frame(accessions = acc, protein_name = protein_name, organism = organism)
}

# Loop through all accessions
results_cluster1  <- lapply(accessions_cluster1, fetch_protein_info)
final_cluster1 <- do.call(rbind, results_cluster1)

results_cluster2  <- lapply(accessions_cluster2, fetch_protein_info)
final_cluster2 <- do.call(rbind, results_cluster2)

results_cluster3  <- lapply(accessions_cluster3, fetch_protein_info)
final_cluster3 <- do.call(rbind, results_cluster3)

results_cluster9  <- lapply(accessions_cluster9, fetch_protein_info)
final_cluster9 <- do.call(rbind, results_cluster9)

## Extract only the protein name from the long GenBank record

extract_protein_name <- function(x) {
  # Try /product="..."
  prod <- sub('.*?/product="([^"]+)".*', '\\1', x)
  if (!is.na(prod) && prod != x) return(prod)
  # Try /name="..."
  name <- sub('.*?/name="([^"]+)".*', '\\1', x)
  if (!is.na(name) && name != x) return(name)
  # Try DEFINITION line
  def <- sub('.*DEFINITION\\s+([^\\[]+?)\\s*\\[.*', '\\1', x)
  if (!is.na(def) && def != x) return(trimws(def))
  # If all fail, return NA
  return(NA)
}

final_cluster1$protein_name_clean <- vapply(final_cluster1$protein_name, extract_protein_name, character(1))
final_cluster2$protein_name_clean <- sub('.*?/product="([^"]+)".*', '\\1', final_cluster2$protein_name)
final_cluster3$protein_name_clean <- vapply(final_cluster3$protein_name, extract_protein_name, character(1))
final_cluster9$protein_name_clean <- vapply(final_cluster9$protein_name, extract_protein_name, character(1))

# Save only accessions and cleaned protein name
write.csv(final_cluster1[, c("accessions", "protein_name_clean")], "shared_hub_genes_cluster1_accession_protein.csv", row.names = FALSE)
write.csv(final_cluster2[, c("accessions", "protein_name_clean")], "shared_hub_genes_cluster2_accession_protein.csv", row.names = FALSE)
write.csv(final_cluster3[, c("accessions", "protein_name_clean")], "shared_hub_genes_cluster3_accession_protein.csv", row.names = FALSE)
write.csv(final_cluster9[, c("accessions", "protein_name_clean")], "shared_hub_genes_cluster9_accession_protein.csv", row.names = FALSE)







### keep only unique hub genes per each temperature and cluster

# create vectors of gene IDs
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

# Get unique hub genes for each temp 30 (non-shared genes)
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

##Load Pcom annotations
final_GO_table_Pcom <- read.csv("final_GO_table_Pcom.csv", header = TRUE, sep = ",")

### filter clusters data data for just expressed genes with GO annotations
cluster1_hub_genes_control_GO <- cluster1_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster1_hub_genes_30_GO <- cluster1_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster1_hub_genes_35_GO <- cluster1_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster2_hub_genes_control_GO <- cluster2_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster2_hub_genes_30_GO <- cluster2_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster2_hub_genes_35_GO <- cluster2_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster3_hub_genes_control_GO <- cluster3_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster3_hub_genes_30_GO <- cluster3_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster3_hub_genes_35_GO <- cluster3_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster9_hub_genes_control_GO <- cluster9_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster9_hub_genes_30_GO <- cluster9_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

cluster9_hub_genes_35_GO <- cluster9_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Pcom$protein_accession)

#% genes with GO annotations
percentage_retained_c1_control <- (nrow(cluster1_hub_genes_control_GO) / nrow(cluster1_hub_gene_counts_control)) * 100
## [1] 71.54

percentage_retained_c1_30 <- (nrow(cluster1_hub_genes_30_GO) / nrow(cluster1_hub_gene_counts_30)) * 100
## [1] 68.37

percentage_retained_c1_35 <- (nrow(cluster1_hub_genes_35_GO) / nrow(cluster1_hub_gene_counts_35)) * 100
## [1] 70.52

percentage_retained_c2_control <- (nrow(cluster2_hub_genes_control_GO) / nrow(cluster2_hub_gene_counts_control)) * 100
## [1] 46.22

percentage_retained_c2_30 <- (nrow(cluster2_hub_genes_30_GO) / nrow(cluster2_hub_gene_counts_30)) * 100
## [1] 35.01

percentage_retained_c2_35 <- (nrow(cluster2_hub_genes_35_GO) / nrow(cluster2_hub_gene_counts_35)) * 100
## [1] 46.22

percentage_retained_c3_control <- (nrow(cluster3_hub_genes_control_GO) / nrow(cluster3_hub_gene_counts_control)) * 100
## [1] 65.71

percentage_retained_c3_30 <- (nrow(cluster3_hub_genes_30_GO) / nrow(cluster3_hub_gene_counts_30)) * 100
## [1] 60.95

percentage_retained_c3_35 <- (nrow(cluster3_hub_genes_35_GO) / nrow(cluster3_hub_gene_counts_35)) * 100
## [1] 64.99

percentage_retained_c9_control <- (nrow(cluster9_hub_genes_control_GO) / nrow(cluster9_hub_gene_counts_control)) * 100
## [1] 68.52

percentage_retained_c9_30 <- (nrow(cluster9_hub_genes_30_GO) / nrow(cluster9_hub_gene_counts_30)) * 100
## [1] 57.41

percentage_retained_c9_35 <- (nrow(cluster9_hub_genes_35_GO) / nrow(cluster9_hub_gene_counts_35)) * 100
## [1] 64.81


## Create custom GO annotation file for ViSEAGO

##Get a list of GO Terms for all clusters
go_terms_list <- final_GO_table_Pcom %>%
  filter(protein_accession %in% cluster1_hub_genes_control_GO$Gene | 
         protein_accession %in% cluster1_hub_genes_30_GO$Gene | 
         protein_accession %in% cluster1_hub_genes_35_GO$Gene |
            protein_accession %in% cluster2_hub_genes_control_GO$Gene |
            protein_accession %in% cluster2_hub_genes_30_GO$Gene |
            protein_accession %in% cluster2_hub_genes_35_GO$Gene |
            protein_accession %in% cluster3_hub_genes_control_GO$Gene |
            protein_accession %in% cluster3_hub_genes_30_GO$Gene |
            protein_accession %in% cluster3_hub_genes_35_GO$Gene |
            protein_accession %in% cluster9_hub_genes_control_GO$Gene |
            protein_accession %in% cluster9_hub_genes_30_GO$Gene |
            protein_accession %in% cluster9_hub_genes_35_GO$Gene) %>%
 dplyr::select(protein_accession, go_ids) %>% dplyr::rename(GO.terms = go_ids) %>% dplyr::rename(query = protein_accession)

##3341


# format into the format required by ViSEAGO for custom mappings
Custom_list_GOs <- go_terms_list %>%
  # Separate GO terms into individual rows
  separate_rows(GO.terms, sep = ";") %>%
  # Add necessary columns
  mutate(
    taxid = "46720",
    gene_symbol = query,
    evidence = "Custom"
  ) %>%
  # Rename columns
  dplyr::rename(
    gene_id = query,
    GOID = GO.terms
  ) %>%
  dplyr::select(taxid, gene_id, gene_symbol, GOID, evidence)

go_count <- sum(stringr::str_count(go_terms_list$GO.terms, "GO"))
go_count
#[1] 299245

Custom_GOs_validated <- Custom_list_GOs %>% filter(GOID %in% keys(GO.db))
# 269,690

#find lost GO terms
lost_gos <- setdiff(Custom_list_GOs, Custom_GOs_validated)
length(lost_gos) 

length(unique(Custom_list_GOs$gene_id))
#[1] 3341
length(unique(Custom_GOs_validated$gene_id))
#[1] 3341

### I didn't loose any gene

write.table(Custom_GOs_validated, "Viseago_custom_GOs_hubGenes_Pcom.txt", row.names = FALSE, sep = "\t", quote = FALSE,col.names=TRUE)


### load into ViSEAGO

Custom_Pcom <- ViSEAGO::Custom2GO("Viseago_custom_GOs_hubGenes_Pcom.txt")

myGENE2GO_Pcom <- ViSEAGO::annotate(
    id="46720",
    Custom_Pcom
)
##46720 is Pcom in NCBI taxonomy


## Create gene lists for enrichment

# load genes background, all expressed genes (from datExpr in WGCNA)
background<-scan(
    "background.txt",
    quiet=TRUE,
    what=""
)

# Only keep expressed genes that have a GO annotation
background <- intersect(
  scan("background.txt", quiet = TRUE, what = ""),
  unique(Custom_GOs_validated$gene_id)
)

## load hub genes for each cluster
# extract all gene IDs
hub_gene_ids_cluster1_control <- cluster1_hub_genes_control_GO$Gene
hub_gene_ids_cluster1_30 <- cluster1_hub_genes_30_GO$Gene
hub_gene_ids_cluster1_35 <- cluster1_hub_genes_35_GO$Gene
hub_gene_ids_cluster2_control <- cluster2_hub_genes_control_GO$Gene
hub_gene_ids_cluster2_30 <- cluster2_hub_genes_30_GO$Gene
hub_gene_ids_cluster2_35 <- cluster2_hub_genes_35_GO$Gene
hub_gene_ids_cluster3_control <- cluster3_hub_genes_control_GO$Gene
hub_gene_ids_cluster3_30 <- cluster3_hub_genes_30_GO$Gene
hub_gene_ids_cluster3_35 <- cluster3_hub_genes_35_GO$Gene
hub_gene_ids_cluster9_control <- cluster9_hub_genes_control_GO$Gene
hub_gene_ids_cluster9_30 <- cluster9_hub_genes_30_GO$Gene
hub_gene_ids_cluster9_35 <- cluster9_hub_genes_35_GO$Gene

# Write to selection.txt
writeLines(hub_gene_ids_cluster1_control, "selection_cluster1_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster1_30, "selection_cluster1_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster1_35, "selection_cluster1_hub_genes_35_GO.txt")
writeLines(hub_gene_ids_cluster2_control, "selection_cluster2_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster2_30, "selection_cluster2_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster2_35, "selection_cluster2_hub_genes_35_GO.txt")
writeLines(hub_gene_ids_cluster3_control, "selection_cluster3_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster3_30, "selection_cluster3_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster3_35, "selection_cluster3_hub_genes_35_GO.txt")
writeLines(hub_gene_ids_cluster9_control, "selection_cluster9_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster9_30, "selection_cluster9_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster9_35, "selection_cluster9_hub_genes_35_GO.txt")

#load for Viseago
selection_cluster1_control<-scan(
    "selection_cluster1_hub_genes_control_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster1_30<-scan(
    "selection_cluster1_hub_genes_30_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster1_35<-scan(
    "selection_cluster1_hub_genes_35_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster2_control<-scan(
    "selection_cluster2_hub_genes_control_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster2_30<-scan(
    "selection_cluster2_hub_genes_30_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster2_35<-scan(
    "selection_cluster2_hub_genes_35_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster3_control<-scan(
    "selection_cluster3_hub_genes_control_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster3_30<-scan(
    "selection_cluster3_hub_genes_30_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster3_35<-scan(
    "selection_cluster3_hub_genes_35_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster9_control<-scan(
    "selection_cluster9_hub_genes_control_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster9_30<-scan(
    "selection_cluster9_hub_genes_30_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster9_35<-scan(
    "selection_cluster9_hub_genes_35_GO.txt",
    quiet=TRUE,
    what=""
)


# Define your gene lists for each cluster/temperature
# gene_lists <- list(
#   cluster1_control = selection_cluster1_control,
#   cluster1_30      = selection_cluster1_30,
#   cluster1_35      = selection_cluster1_35,
#   cluster2_control = selection_cluster2_control,
#   cluster2_30      = selection_cluster2_30,
#   cluster2_35      = selection_cluster2_35,
#   cluster3_control = selection_cluster3_control,
#   cluster3_30      = selection_cluster3_30,
#   cluster3_35      = selection_cluster3_35,
#   cluster9_control = selection_cluster9_control,
#   cluster9_30      = selection_cluster9_30,
#   cluster9_35      = selection_cluster9_35
# )  ### I'll run the enrichment separately for cluster 9, as it creates problem in the loop (likely due to low enriched GO terms)

gene_lists <- list(
  cluster1_control = selection_cluster1_control,
  cluster1_30      = selection_cluster1_30,
  cluster1_35      = selection_cluster1_35,
  cluster2_control = selection_cluster2_control,
  cluster2_30      = selection_cluster2_30,
  cluster2_35      = selection_cluster2_35,
  cluster3_control = selection_cluster3_control,
  cluster3_30      = selection_cluster3_30,
  cluster3_35      = selection_cluster3_35
)

library(ViSEAGO)

results_list <- list()

for (name in names(gene_lists)) {
  for (ont in c("BP", "MF")) {
    topgo <- ViSEAGO::create_topGOdata(
      geneSel = gene_lists[[name]],
      allGenes = background,
      gene2GO = myGENE2GO_Pcom,
      ont = ont,
      nodeSize = 5
    )
    res <- topGO::runTest(topgo, algorithm = "classic", statistic = "fisher")
    pvals <- topGO::score(res)
    # pvals_adj <- p.adjust(pvals, method = "BH")
    # sig_terms <- names(pvals_adj)[pvals_adj < 0.05]
    sig_terms <- names(pvals)[pvals < 0.01] 
    if (length(sig_terms) == 0) {
      message(paste("No enriched GO terms for", name, ont, "- skipping."))
      next
    }
    enrich <- ViSEAGO::merge_enrich_terms(Input = list(name = c("topgo", "res")))
    if (is.null(enrich@data) || nrow(enrich@data) == 0) {
      message(paste("No enriched GO terms for", name, ont, "- skipping."))
      next
    }
    ss <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = enrich)
    ss <- ViSEAGO::compute_SS_distances(ss, distance = "Wang")
    heatmap <- ViSEAGO::GOterms_heatmap(ss, showIC = TRUE, showGOlabels = TRUE)
    # Only assign if everything succeeded
    results_list[[paste0(name, "_", ont)]] <- list(
      topgo = topgo,
      enrich = enrich,
      ss = ss,
      heatmap = heatmap
#     # Only run if heatmap object is valid for table extraction
# if ("x" %in% slotNames(heatmap) && !is.null(heatmap@x$data) && nrow(heatmap@x$data) > 0) {
#   ViSEAGO::show_table(heatmap, paste0(ont, "_", name, "_table.xls"))
#   pdf(paste0(ont, "_", name, "_heatmap.pdf"), width = 10, height = 8)
#   ViSEAGO::show_heatmap(heatmap, "GOterms")
#   dev.off()
# }
#   }
# }
    )
    ViSEAGO::show_table(heatmap, paste0(ont, "_", name, "_table.xls"))
    pdf(paste0(ont, "_", name, "_heatmap.pdf"), width = 10, height = 8)
    ViSEAGO::show_heatmap(heatmap, "GOterms")
    dev.off()
  }
}

#After running the loop, check which entries are non-empty
valid_results <- results_list[!sapply(results_list, is.null)]
length(valid_results)
names(valid_results)



##### Cluster 9 run enrichment separately
# create viseago object
cluster9_control_BP <- ViSEAGO::create_topGOdata(
    geneSel=selection_cluster9_control,
    allGenes=background,
    gene2GO=myGENE2GO_Pcom, 
    ont="BP",
    nodeSize=5
)

cluster9_control_MF <- ViSEAGO::create_topGOdata(
    geneSel=selection_cluster9_control,
    allGenes=background,
    gene2GO=myGENE2GO_Pcom, 
    ont="MF",
    nodeSize=5
)

cluster9_30_BP <- ViSEAGO::create_topGOdata(
    geneSel=selection_cluster9_30,
    allGenes=background,
    gene2GO=myGENE2GO_Pcom, 
    ont="BP",
    nodeSize=5
)

cluster9_30_MF <- ViSEAGO::create_topGOdata(
    geneSel=selection_cluster9_30,
    allGenes=background,
    gene2GO=myGENE2GO_Pcom, 
    ont="MF",
    nodeSize=5
)

cluster9_35_BP <- ViSEAGO::create_topGOdata(
    geneSel=selection_cluster9_35,
    allGenes=background,
    gene2GO=myGENE2GO_Pcom, 
    ont="BP",
    nodeSize=5
)

cluster9_35_MF <- ViSEAGO::create_topGOdata(
    geneSel=selection_cluster9_35,
    allGenes=background,
    gene2GO=myGENE2GO_Pcom, 
    ont="MF",
    nodeSize=5
)

# perform TopGO test using classic algorithm
classic_cluster9_control_BP <- topGO::runTest(
    cluster9_control_BP,
    algorithm ="classic",
    statistic = "fisher"
)

pvals <- topGO::score(classic_cluster9_control_BP)
# pvals_adj <- p.adjust(pvals, method = "BH")
# sig_terms <- names(pvals_adj)[pvals_adj < 0.05]
sig_terms <- names(pvals)[pvals < 0.01] 

BP_Results_control_cluster9 <- ViSEAGO::merge_enrich_terms(
    Input = list(cluster9 = c("cluster9_control_BP", "classic_cluster9_control_BP"))
)

BP_Results_control_cluster9
#  - enrichment pvalue cutoff:
#         cluster9 : 0.01
# - enrich GOs (in at least one list): 65 GO terms of 1 conditions.
#         cluster9 : 65 terms

classic_cluster9_control_MF <- topGO::runTest(
    cluster9_control_MF,
    algorithm ="classic",
    statistic = "fisher"
)

MF_Results_control_cluster9 <- ViSEAGO::merge_enrich_terms(
    Input = list(cluster9 = c("cluster9_control_MF", "classic_cluster9_control_MF"))
)

MF_Results_control_cluster9
# - enrich GOs (in at least one list): 15 GO terms of 1 conditions.
#         cluster9 : 15 terms

classic_cluster9_30_BP <- topGO::runTest(
    cluster9_30_BP,
    algorithm ="classic",
    statistic = "fisher"
)

BP_Results_30_cluster9 <- ViSEAGO::merge_enrich_terms(
    Input = list(cluster9 = c("cluster9_30_BP", "classic_cluster9_30_BP"))
)

BP_Results_30_cluster9
# - enrich GOs (in at least one list): 55 GO terms of 1 conditions.
#         cluster9 : 55 terms

classic_cluster9_30_MF <- topGO::runTest(
    cluster9_30_MF,
    algorithm ="classic",
    statistic = "fisher"
)

MF_Results_30_cluster9 <- ViSEAGO::merge_enrich_terms(
    Input = list(cluster9 = c("cluster9_30_MF", "classic_cluster9_30_MF"))
)

MF_Results_30_cluster9
# - enrich GOs (in at least one list): 10 GO terms of 1 conditions.
#         cluster9 : 10 terms

classic_cluster9_35_BP <- topGO::runTest(
    cluster9_35_BP,
    algorithm ="classic",
    statistic = "fisher"
)

BP_Results_35_cluster9 <- ViSEAGO::merge_enrich_terms(
    Input = list(cluster9 = c("cluster9_35_BP", "classic_cluster9_35_BP"))
)

BP_Results_35_cluster9
# - enrich GOs (in at least one list): 16 GO terms of 1 conditions.
#         cluster9 : 16 terms

classic_cluster9_35_MF <- topGO::runTest(
    cluster9_35_MF,
    algorithm ="classic",
    statistic = "fisher"
)

MF_Results_35_cluster9 <- ViSEAGO::merge_enrich_terms(
    Input = list(cluster9 = c("cluster9_35_MF", "classic_cluster9_35_MF"))
)

MF_Results_35_cluster9
# - enrich GOs (in at least one list): 21 GO terms of 1 conditions.
#         cluster9 : 21 terms


#### add cluster 9 to valid results
# For each cluster9 condition, build the heatmap and add to valid_results
# BP control
ss_control_BP <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = BP_Results_control_cluster9)
ss_control_BP <- ViSEAGO::compute_SS_distances(ss_control_BP, distance = "Wang")
heatmap_control_BP <- ViSEAGO::GOterms_heatmap(ss_control_BP, showIC = TRUE, showGOlabels = TRUE)

valid_results[["cluster9_control_BP"]] <- list(
  topgo = cluster9_control_BP,
  enrich = BP_Results_control_cluster9,
  ss = ss_control_BP,
  heatmap = heatmap_control_BP
) ##this one gives issues
# Even though there are 65 enriched GO terms, ViSEAGO could not compute semantic similarity distances between them.

heatmap_control_BP <- tryCatch(
  ViSEAGO::GOterms_heatmap(ss_control_BP, showIC = TRUE, showGOlabels = TRUE),
  error = function(e) {
    message("Heatmap failed: ", e$message)
    NULL
  }
)

valid_results[["cluster9_control_BP"]] <- list(
  topgo = cluster9_control_BP,
  enrich = BP_Results_control_cluster9,
  ss = ss_control_BP,
  heatmap = heatmap_control_BP # will be NULL if failed
)

# MF control
ss_control_MF <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = MF_Results_control_cluster9)
ss_control_MF <- ViSEAGO::compute_SS_distances(ss_control_MF, distance = "Wang")
heatmap_control_MF <- ViSEAGO::GOterms_heatmap(ss_control_MF, showIC = TRUE, showGOlabels = TRUE)

valid_results[["cluster9_control_MF"]] <- list(
  topgo = cluster9_control_MF,
  enrich = MF_Results_control_cluster9,
  ss = ss_control_MF,
  heatmap = heatmap_control_MF
)

# BP 30
ss_30_BP <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = BP_Results_30_cluster9)
ss_30_BP <- ViSEAGO::compute_SS_distances(ss_30_BP, distance = "Wang")
heatmap_30_BP <- ViSEAGO::GOterms_heatmap(ss_30_BP, showIC = TRUE, showGOlabels = TRUE)

valid_results[["cluster9_30_BP"]] <- list(
  topgo = cluster9_30_BP,
  enrich = BP_Results_30_cluster9,
  ss = ss_30_BP,
  heatmap = heatmap_30_BP
)

# MF 30
ss_30_MF <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = MF_Results_30_cluster9)
ss_30_MF <- ViSEAGO::compute_SS_distances(ss_30_MF, distance = "Wang")
heatmap_30_MF <- ViSEAGO::GOterms_heatmap(ss_30_MF, showIC = TRUE, showGOlabels = TRUE)

valid_results[["cluster9_30_MF"]] <- list(
  topgo = cluster9_30_MF,
  enrich = MF_Results_30_cluster9,
  ss = ss_30_MF,
  heatmap = heatmap_30_MF
)

# BP 35
ss_35_BP <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = BP_Results_35_cluster9)
ss_35_BP <- ViSEAGO::compute_SS_distances(ss_35_BP, distance = "Wang")
heatmap_35_BP <- ViSEAGO::GOterms_heatmap(ss_35_BP, showIC = TRUE, showGOlabels = TRUE)

valid_results[["cluster9_35_BP"]] <- list(
  topgo = cluster9_35_BP,
  enrich = BP_Results_35_cluster9,
  ss = ss_35_BP,
  heatmap = heatmap_35_BP
)

# MF 35
ss_35_MF <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Pcom, enrich_GO_terms = MF_Results_35_cluster9)
ss_35_MF <- ViSEAGO::compute_SS_distances(ss_35_MF, distance = "Wang")
heatmap_35_MF <- ViSEAGO::GOterms_heatmap(ss_35_MF, showIC = TRUE, showGOlabels = TRUE)

valid_results[["cluster9_35_MF"]] <- list(
  topgo = cluster9_35_MF,
  enrich = MF_Results_35_cluster9,
  ss = ss_35_MF,
  heatmap = heatmap_35_MF
)



##### plot top 10 enriched GO terms, facet by temperature
library(dplyr)
library(ggplot2)
library(forcats)

# Extract enrichment info from ViSEAGO results
compare_df <- bind_rows(
  lapply(valid_results, function(res) {
    if (!is.null(res$heatmap)) {
      dt_obj <- ViSEAGO::show_table(res$heatmap)
      df <- as.data.frame(dt_obj$x$data)
      df$ID <- sub(".*GO:(\\d+).*", "GO:\\1", df$`GO ID`)
      df
    } else {
      NULL
    }
  }),
  .id = "cluster_condition"
)

#put pvalues and -log10_pvalues in the same columns for cluster 9
compare_df <- compare_df %>%
  mutate(
    `name pvalue` = coalesce(`name pvalue`, `cluster9 pvalue`),
    `name -log10_pvalue` = coalesce(`name -log10_pvalue`, `cluster9 -log10_pvalue`)
  )

# # Count number of genes per GO term 
# For each cluster, filter Custom_GOs_validated to include only genes in that cluster's gene list.
# Count the number of genes per GO term for that cluster.
# Join this cluster-specific gene count to the enrichment results for plotting.

cluster_gene_counts_list <- lapply(names(gene_lists), function(clust_name) {
  cluster_genes <- gene_lists[[clust_name]]
  cluster_GO_annot <- Custom_GOs_validated %>%
    filter(gene_id %in% cluster_genes)
  cluster_gene_counts <- cluster_GO_annot %>%
    group_by(GOID) %>%
    summarise(gene_count = n())
  cluster_gene_counts$cluster_condition <- clust_name
  cluster_gene_counts
})

all_cluster_gene_counts <- bind_rows(cluster_gene_counts_list)

# Extract gene counts for cluster9 (control, 30, 35)
for (clust_name in c("cluster9_control", "cluster9_30", "cluster9_35")) {
  cluster_genes <- get(paste0("selection_", clust_name))
  cluster_GO_annot <- Custom_GOs_validated %>%
    filter(gene_id %in% cluster_genes)
  cluster_gene_counts <- cluster_GO_annot %>%
    group_by(GOID) %>%
    summarise(gene_count = n())
  cluster_gene_counts$cluster_condition <- clust_name
  cluster_gene_counts_list[[clust_name]] <- cluster_gene_counts
}

# Combine all clusters
all_cluster_gene_counts <- bind_rows(cluster_gene_counts_list)

# Add ontology column and clean cluster_condition in compare_df
compare_df <- compare_df %>%
  mutate(
    temperature = case_when(
      grepl("control", cluster_condition) ~ "control",
      grepl("30", cluster_condition) ~ "30",
      grepl("35", cluster_condition) ~ "35",
      TRUE ~ NA_character_
    ),
    temperature = factor(temperature, levels = c("control", "30", "35")),
    ontology = ifelse(grepl("_BP$", cluster_condition), "BP",
                      ifelse(grepl("_MF$", cluster_condition), "MF", NA)),
    cluster_condition = sub("(_BP|_MF)$", "", cluster_condition)
  )

# Join to compare_df
compare_df <- compare_df %>%
  left_join(all_cluster_gene_counts, by = c("ID" = "GOID", "cluster_condition" = "cluster_condition"))

# Get top 10 GO terms per cluster
N <- 10
top_terms <- compare_df %>%
  group_by(cluster_condition, temperature) %>%
  arrange(`name pvalue`) %>%
  slice_head(n = N) %>%
  ungroup()

# Split cluster_condition into cluster and temperature
top_terms <- top_terms %>%
  mutate(
    cluster = sub("(_control|_30|_35)$", "", cluster_condition),
    temperature = case_when(
      grepl("_control$", cluster_condition) ~ "control",
      grepl("_30$", cluster_condition) ~ "30",
      grepl("_35$", cluster_condition) ~ "35",
      TRUE ~ NA_character_
    ),
    temperature = factor(temperature, levels = c("control", "30", "35"))
  )

ggplot(top_terms, aes(
    x = cluster,
    y = forcats::fct_reorder2(term, cluster, `name pvalue`),
    size = gene_count,
    color = -log10(`name pvalue`)
  )) +
  geom_point(alpha = 0.8) +
scale_size_continuous(
  range = c(1, 8),
  breaks = c(10, 20, 50, 100, 200, 400)
) +
  scale_color_gradient(low = "moccasin", high = "darkorange") +
  labs(
    size = "Gene Count",
    color = "-log10(p-value)"
  ) +
  scale_x_discrete(expand = c(0.2, 0.2)) + # less space on x-axis
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1) # less margin
  ) +
  facet_wrap(~ temperature, ncol = 1, scales = "free_y")

ggsave("ViSEAGO_top10_GOterms_faceted_dotplot.pdf", width = 8, height = 25)


# Save compare_df
# Split cluster_condition into cluster and temperature
compare_df <- compare_df %>%
  mutate(
    cluster = sub("(_control|_30|_35)$", "", cluster_condition),
    temperature = case_when(
      grepl("_control$", cluster_condition) ~ "control",
      grepl("_30$", cluster_condition) ~ "30",
      grepl("_35$", cluster_condition) ~ "35",
      TRUE ~ NA_character_
    ),
    temperature = factor(temperature, levels = c("control", "30", "35"))
  )

# Add direction column to compare_df
compare_df <- compare_df %>%
  mutate(
    direction = case_when(
      cluster %in% c("cluster3", "cluster9") ~ "up",
      cluster %in% c("cluster1", "cluster2") ~ "down",
      TRUE ~ NA_character_
    )
  )

write.csv(compare_df, "Sign_GO_p01Fisher_Pcom.csv", row.names = FALSE)

# Save top_terms
# Add direction column to top_terms
top_terms <- top_terms %>%
  mutate(
    direction = case_when(
      cluster %in% c("cluster3", "cluster9") ~ "up",
      cluster %in% c("cluster1", "cluster2") ~ "down",
      TRUE ~ NA_character_
    )
  )

#remove empty columns
top_terms <- top_terms %>%
  select(-`cluster9 pvalue`, -`cluster9 -log10_pvalue`)

write.csv(top_terms, "Top10_GO_p01Fisher_Pcom.csv", row.names = FALSE)




###### GO Slim Analysis
library(GSEABase)

# Download GO slim OBO file and read as a GeneSetCollection
goslim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo")
#format-version: 1.2

# Extract all enriched GO terms from your results
enriched_go_terms <- unique(unlist(
  lapply(valid_results, function(res) {
    if (!is.null(res$heatmap)) {
      dt_obj <- tryCatch(
        ViSEAGO::show_table(res$heatmap),
        error = function(e) NULL
      )
      if (!is.null(dt_obj)) {    ## skip clusters with NULL or invalid heatmap
        df <- as.data.frame(dt_obj$x$data)
        sub(".*GO:(\\d+).*", "GO:\\1", df$`GO ID`)
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
))

# Create a GOCollection object from your enriched GO terms
GO_collection <- GOCollection(enriched_go_terms)

# Map to BP slim categories
slims_BP <- data.frame(goSlim(GO_collection, goslim, "BP"))
slims_BP$category <- row.names(slims_BP)

slims_MF <- data.frame(goSlim(GO_collection, goslim, "MF"))
slims_MF$category <- row.names(slims_MF)

#Get mapped terms, using function mappedIds to get the query terms that mapped to the slim categories
mappedIds <-
  function(df, collection, OFFSPRING) #the command to run requires a dataframe of slim terms, like slims_MF above, your list of query terms, and the offspring from the GOCollection by goSlim
  {
    map <- as.list(OFFSPRING)[rownames(df)] # Subset GOcollection offspring by the rownames of your dataframe
    mapped <- lapply(map, intersect, ids(collection)) #Find the terms that intersect between the subset made above of your query terms and the GOids from the GO collection
    df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L)) #Add column "go_terms" with matching terms 
    df #show resulting dataframe
  }

#Run function for MF and BP terms
BPslim <- mappedIds(slims_BP, GO_collection, GOBPOFFSPRING)
MFslim <- mappedIds(slims_MF, GO_collection, GOMFOFFSPRING)

#Remove duplicate matches, keeping the broader umbrella term
#BP
BPslim <- filter(BPslim, Count>0 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(BPslim$go_terms), ";") #split into multiple GO ids
BPslimX <- data.frame(Term = rep.int(BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
BPslimX <- merge(BPslimX, BPslim[,c(1,3:4)], by="Term") #Add back counts, term, and category info
BPslimX <- unique(setDT(BPslimX)[order(go_term, -Count)], by = "go_term") #remove duplicate offspring terms, keeping only those that appear in the larger umbrella term (larger Count number)
BPslim <- data.frame(slim_term=BPslimX$Term, slim_cat=BPslimX$category, category=BPslimX$go_term) #rename columns
head(BPslim)

#MF
MFslim <- filter(MFslim, Count>0 & Term!="molecular_function") #filter out empty slims and term "molecular function"
MFsplitted <- strsplit(as.character(MFslim$go_terms), ";") #split into multiple GO ids
MFslimX <- data.frame(Term = rep.int(MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
MFslimX <- merge(MFslimX, MFslim[,c(1,3:4)], by="Term")  #Add back counts, term, and category info
MFslimX <- unique(setDT(MFslimX)[order(go_term, -Count)], by = "go_term")  #remove duplicate offspring terms, keeping only
MFslim <- data.frame(slim_term=MFslimX$Term, slim_cat=MFslimX$category, category=MFslimX$go_term) #rename columns
head(MFslim)

###Save slim info with GO enrichment info for heatmap dataframes
#rename category column to ID
colnames(BPslim)[3] <- "ID"
colnames(MFslim)[3] <- "ID"

#add column with BP or MF named ontology
BPslim$ontology <- "BP"
MFslim$ontology <- "MF"

# Extract enrichment info from ViSEAGO results
compare_df <- bind_rows(
  lapply(valid_results, function(res) {
    if (!is.null(res$heatmap)) {
      dt_obj <- tryCatch(
        ViSEAGO::show_table(res$heatmap),
        error = function(e) NULL
      )
      if (!is.null(dt_obj)) {
        df <- as.data.frame(dt_obj$x$data)
        # Clean GO IDs
        df$ID <- sub(".*GO:(\\d+).*", "GO:\\1", df$`GO ID`)
        df
      } else {
        NULL
      }
    } else {
      NULL
    }
  }),
  .id = "cluster_condition"
)

#remove empty columns
compare_df <- compare_df %>%
  select(-`cluster9 pvalue`, -`cluster9 -log10_pvalue`)

GO.BP <- right_join(BPslim, compare_df, by="ID") #add back GO enrichment info for each offspring term
GO.MF <- right_join(MFslim, compare_df, by="ID") #add back GO enrichment info for each offspring term

#remove empty columns
GO.BP <- GO.BP %>%
  select(-`name pvalue`, -`name -log10_pvalue`)

GO.MF <- GO.MF %>%
  select(-`name pvalue`, -`name -log10_pvalue`)

GO.BP <- na.omit(GO.BP) 
GO.MF <- na.omit(GO.MF) 

### Dotplot: GO slim term enrichment across clusters
library(reshape2)

GO.BP %>%
  group_by(cluster_condition, slim_term) %>%
  summarise(go_term_count = n_distinct(ID)) %>%
  arrange(desc(go_term_count)) %>%
  ggplot(aes(x = cluster_condition, y = reorder(slim_term, go_term_count), size = go_term_count, color = go_term_count)) +
    geom_point() +
    labs(title = "GO Slim BP Dotplot (Number of GO Terms per Cluster)",
         x = "Cluster",
         y = "GO Slim Term",
         size = "GO Term Count",
         color = "GO Term Count") +
    theme_bw()

GO.MF %>%
  group_by(cluster_condition, slim_term) %>%
  summarise(go_term_count = n_distinct(ID)) %>%
  arrange(desc(go_term_count)) %>%
  ggplot(aes(x = cluster_condition, y = reorder(slim_term, go_term_count), size = go_term_count, color = go_term_count)) +
    geom_point() +
    labs(title = "GO Slim MF Dotplot (Number of GO Terms per Cluster)",
         x = "Cluster",
         y = "GO Slim Term",
         size = "GO Term Count",
         color = "GO Term Count") +
    theme_bw()

###combine GO.BP and GO.MF into a single data frame and plot
# Combine BP and MF slim results
GO.slim <- bind_rows(GO.BP, GO.MF)

# Summarize by number of GO terms per slim term, cluster, and ontology
GO.slim_summary <- GO.slim %>%
  group_by(cluster_condition, slim_term, ontology) %>%
  summarise(go_term_count = n_distinct(ID), .groups = "drop")

# Dotplot: show BP and MF together, colored by ontology
GO.slim_summary <- GO.slim %>%
  mutate(cluster_condition = sub("(_BP|_MF)$", "", cluster_condition)) %>%  # Remove _BP or _MF at end
  group_by(cluster_condition, slim_term, ontology) %>%
  summarise(go_term_count = n_distinct(ID), .groups = "drop")

p <- ggplot(GO.slim_summary, aes(x = cluster_condition, y = reorder(slim_term, go_term_count), 
                            size = go_term_count, color = ontology)) +
  geom_point(alpha = 0.8) +
  labs(title = "GO Slim Dotplot (BP & MF)",
       x = "Cluster",
       y = "GO Slim Term",
       size = "GO Term Count",
       color = "Ontology") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
ggsave("GO_slim_dotplot.pdf", plot = p, width = 10, height = 8)

# Create a custom order for cluster_condition
GO.slim_summary <- GO.slim %>%
  mutate(cluster_condition = sub("(_BP|_MF)$", "", cluster_condition)) %>%
  group_by(cluster_condition, slim_term, ontology) %>%
  summarise(go_term_count = n_distinct(ID), .groups = "drop") %>%
  mutate(
    cluster_condition = factor(
      cluster_condition,
      levels = c(
        grep("control", unique(cluster_condition), value = TRUE),
        grep("30", unique(cluster_condition), value = TRUE),
        grep("35", unique(cluster_condition), value = TRUE)
      )
    )
  )

p <- ggplot(GO.slim_summary, aes(x = cluster_condition, y = reorder(slim_term, go_term_count), 
                            size = go_term_count, color = ontology)) +
  geom_point(alpha = 0.8) +
  labs(title = "GO Slim Dotplot (BP & MF)",
       x = "Cluster",
       y = "GO Slim Term",
       size = "GO Term Count",
       color = "Ontology") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
ggsave("GO_slim_dotplot.pdf", plot = p, width = 10, height = 8)


### facets plots
GO.slim_summary <- GO.slim_summary %>%
  mutate(
    temperature = case_when(
      grepl("control", cluster_condition) ~ "control",
      grepl("30", cluster_condition) ~ "30",
      grepl("35", cluster_condition) ~ "35",
      TRUE ~ NA_character_
    ),
    temperature = factor(temperature, levels = c("control", "30", "35"))
  )

#Faceted dotplot by temperature
library(ggplot2)
library(dplyr)
library(cowplot)

# Ensure temperature is a factor with correct order
GO.slim_summary$temperature <- factor(
  GO.slim_summary$temperature,
  levels = c("control", "30", "35")
)

# For each temperature, plot only the top N GO slim terms per cluster
#N <- 20
plots <- lapply(levels(GO.slim_summary$temperature), function(temp) {
  df_temp <- GO.slim_summary %>%
    filter(temperature == temp) %>%
    group_by(cluster_condition) %>%
    arrange(desc(go_term_count)) %>%
    #slice_head(n = N) %>%
    ungroup()
  if (nrow(df_temp) > 0) {
ggplot(df_temp, aes(x = cluster_condition, y = reorder(slim_term, go_term_count))) +
  geom_point(aes(size = go_term_count, color = ontology)) +
  scale_x_discrete(expand = c(0.2, 0.2)) + # less space on x-axis
  scale_size_continuous(range = c(2, 8)) + # smaller dots if needed
  labs(
    x = NULL,
    y = paste("GO Slim Term (", temp, ")", sep = ""),
    size = "GO Term Count",
    color = "Ontology"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1), # align y labels
    plot.title = element_blank(),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2) # less margin
  )
  } else {
    NULL
  }
})

# Remove NULL plots
plots <- Filter(Negate(is.null), plots)

# Extract legend from the last plot
legend <- cowplot::get_legend(plots[[length(plots)]])

# Remove legends from all plots
plots_nolegend <- lapply(plots, function(p) p + theme(legend.position = "none"))

# Combine plots vertically and add legend on the right
cowplot::plot_grid(
  cowplot::plot_grid(plotlist = plots_nolegend, ncol = 1, align = "v"),
  legend,
  rel_widths = c(4, 1)
)

pdf("GO_slim_faceted_dotplot_Pcom.pdf", width = 8, height = 15)
print(
  cowplot::plot_grid(
    cowplot::plot_grid(plotlist = plots_nolegend, ncol = 1, align = "v"),
    legend,
    rel_widths = c(4, 1)
  )
)
dev.off()


# Save GO slims
# Split cluster_condition into cluster and temperature
GO.slim_summary <- GO.slim_summary %>%
  mutate(
    cluster = sub("(_control|_30|_35)$", "", cluster_condition),
    temperature = case_when(
      grepl("_control$", cluster_condition) ~ "control",
      grepl("_30$", cluster_condition) ~ "30",
      grepl("_35$", cluster_condition) ~ "35",
      TRUE ~ NA_character_
    ),
    temperature = factor(temperature, levels = c("control", "30", "35"))
  )

# Add direction column to GO.slim_summary
GO.slim_summary <- GO.slim_summary %>%
  mutate(
    direction = case_when(
      cluster %in% c("cluster3", "cluster9") ~ "up",
      cluster %in% c("cluster1", "cluster2") ~ "down",
      TRUE ~ NA_character_
    )
  )

write.csv(GO.slim_summary, "GOslim_All_p01Fisher_Pcom.csv", row.names = FALSE)
