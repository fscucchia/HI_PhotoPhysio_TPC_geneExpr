
####### GO Enrichment Analysis of WGCNA top 10% hub genes per each cluster - Mcap

setwd("...................")
load(".RData") #load WCGNA data

#Load libraries
BiocManager::install("ViSEAGO") # Install ViSEAGO package from Bioconductor
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
treatmentinfo <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/RNAseq_Mcap_data.csv", header = TRUE, sep = ";")

#keep only control, 30 and 35°C
treatmentinfo <- treatmentinfo[treatmentinfo$temp %in% c(26.8, 30, 35), ]

##load gene clusters

#cluster5
cluster5_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_control_cluster5_2.csv")
cluster5_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat30_cluster5_2.csv")
cluster5_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat35_cluster5_2.csv")

#cluster8
cluster8_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_control_cluster8_2.csv")
cluster8_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat30_cluster8_2.csv")
cluster8_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat35_cluster8_2.csv")

#cluster9
cluster9_hub_gene_counts_control <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_control_cluster9_2.csv")
cluster9_hub_gene_counts_30 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat30_cluster9_2.csv")
cluster9_hub_gene_counts_35 <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat35_cluster9_2.csv")

### keep only unique hub genes per each temperature and cluster

# create vectors of gene IDs
cluster5_hub_genes_control <- cluster5_hub_gene_counts_control$Gene
cluster5_hub_genes_30 <- cluster5_hub_gene_counts_30$Gene
cluster5_hub_genes_35 <- cluster5_hub_gene_counts_35$Gene

cluster8_hub_genes_control <- cluster8_hub_gene_counts_control$Gene
cluster8_hub_genes_30 <- cluster8_hub_gene_counts_30$Gene
cluster8_hub_genes_35 <- cluster8_hub_gene_counts_35$Gene

cluster9_hub_genes_control <- cluster9_hub_gene_counts_control$Gene
cluster9_hub_genes_30 <- cluster9_hub_gene_counts_30$Gene
cluster9_hub_genes_35 <- cluster9_hub_gene_counts_35$Gene

# Get unique hub genes for each temp 30 (non-shared genes)
cluster5_hub_genes_30_unique <- setdiff(
  cluster5_hub_genes_30,
  union(cluster5_hub_genes_control, cluster5_hub_genes_35)
)

cluster5_hub_genes_35_unique <- setdiff(
  cluster5_hub_genes_35,
  union(cluster5_hub_genes_control, cluster5_hub_genes_30)
)

cluster5_hub_genes_control_unique <- setdiff(
  cluster5_hub_genes_control,
  union(cluster5_hub_genes_30, cluster5_hub_genes_35)
)

cluster8_hub_genes_30_unique <- setdiff(
  cluster8_hub_genes_30,
  union(cluster8_hub_genes_control, cluster8_hub_genes_35)
)

cluster8_hub_genes_35_unique <- setdiff(
  cluster8_hub_genes_35,
  union(cluster8_hub_genes_control, cluster8_hub_genes_30)
)

cluster8_hub_genes_control_unique <- setdiff(
  cluster8_hub_genes_control,
  union(cluster8_hub_genes_30, cluster8_hub_genes_35)
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


##Load Mcap annotations
final_GO_table_Mcap <- read.csv("final_GO_table_Mcap.csv", header = TRUE, sep = ",")

### filter clusters data data for just expressed genes with GO annotations
cluster5_hub_genes_control_GO <- cluster5_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster5_hub_genes_30_GO <- cluster5_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster5_hub_genes_35_GO <- cluster5_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster8_hub_genes_control_GO <- cluster8_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster8_hub_genes_30_GO <- cluster8_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster8_hub_genes_35_GO <- cluster8_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster9_hub_genes_control_GO <- cluster9_hub_gene_counts_control %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster9_hub_genes_30_GO <- cluster9_hub_gene_counts_30 %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)

cluster9_hub_genes_35_GO <- cluster9_hub_gene_counts_35 %>%
  filter(Gene %in% final_GO_table_Mcap$protein_accession)


#% genes with GO annotations
percentage_retained_c5_control <- (nrow(cluster5_hub_genes_control_GO) / nrow(cluster5_hub_gene_counts_control)) * 100
## [1] 62.25%

percentage_retained_c5_30 <- (nrow(cluster5_hub_genes_30_GO) / nrow(cluster5_hub_gene_counts_30)) * 100
## [1] 59.54%

percentage_retained_c5_35 <- (nrow(cluster5_hub_genes_35_GO) / nrow(cluster5_hub_gene_counts_35)) * 100
## [1] 59.24%

percentage_retained_c8_control <- (nrow(cluster8_hub_genes_control_GO) / nrow(cluster8_hub_gene_counts_control)) * 100
## [1] 55.47%

percentage_retained_c8_30 <- (nrow(cluster8_hub_genes_30_GO) / nrow(cluster8_hub_gene_counts_30)) * 100
## [1] 53.57%

percentage_retained_c8_35 <- (nrow(cluster8_hub_genes_35_GO) / nrow(cluster8_hub_gene_counts_35)) * 100
## [1] 57.84%

percentage_retained_c9_control <- (nrow(cluster9_hub_genes_control_GO) / nrow(cluster9_hub_gene_counts_control)) * 100
## [1] 68.45%

percentage_retained_c9_30 <- (nrow(cluster9_hub_genes_30_GO) / nrow(cluster9_hub_gene_counts_30)) * 100
## [1] 64.56%

percentage_retained_c9_35 <- (nrow(cluster9_hub_genes_35_GO) / nrow(cluster9_hub_gene_counts_35)) * 100
## [1] 62.14%


## Create custom GO annotation file for ViSEAGO

##Get a list of GO Terms for all clusters
go_terms_list <- final_GO_table_Mcap %>%
  filter(protein_accession %in% cluster5_hub_genes_control_GO$Gene | 
         protein_accession %in% cluster5_hub_genes_30_GO$Gene | 
         protein_accession %in% cluster5_hub_genes_35_GO$Gene |
            protein_accession %in% cluster8_hub_genes_control_GO$Gene |
            protein_accession %in% cluster8_hub_genes_30_GO$Gene |
            protein_accession %in% cluster8_hub_genes_35_GO$Gene |
            protein_accession %in% cluster9_hub_genes_control_GO$Gene |
            protein_accession %in% cluster9_hub_genes_30_GO$Gene |
            protein_accession %in% cluster9_hub_genes_35_GO$Gene) %>%
 dplyr::select(protein_accession, go_ids) %>% dplyr::rename(GO.terms = go_ids) %>% dplyr::rename(query = protein_accession)

##2950


# format into the format required by ViSEAGO for custom mappings
Custom_list_GOs <- go_terms_list %>%
  # Separate GO terms into individual rows
  separate_rows(GO.terms, sep = ";") %>%
  # Add necessary columns
  mutate(
    taxid = "46704",
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
#[1] 330782

Custom_GOs_validated <- Custom_list_GOs %>% filter(GOID %in% keys(GO.db))
# 299,576

#find lost GO terms
lost_gos <- setdiff(Custom_list_GOs, Custom_GOs_validated)
length(lost_gos) 

length(unique(Custom_list_GOs$gene_id))
#[1] 2950
length(unique(Custom_GOs_validated$gene_id))
#[1] 2950

### I didn't loose any gene

write.table(Custom_GOs_validated, "Viseago_custom_GOs_hubGenes.txt", row.names = FALSE, sep = "\t", quote = FALSE,col.names=TRUE)


### load into ViSEAGO

Custom_Mcap <- ViSEAGO::Custom2GO("Viseago_custom_GOs_hubGenes.txt")

myGENE2GO_Mcap <- ViSEAGO::annotate(
    id="46704",
    Custom_Mcap
)
##46704 is Mcap in NCBI taxonomy


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
hub_gene_ids_cluster5_control <- cluster5_hub_genes_control_GO$Gene
hub_gene_ids_cluster5_30 <- cluster5_hub_genes_30_GO$Gene
hub_gene_ids_cluster5_35 <- cluster5_hub_genes_35_GO$Gene
hub_gene_ids_cluster8_control <- cluster8_hub_genes_control_GO$Gene
hub_gene_ids_cluster8_30 <- cluster8_hub_genes_30_GO$Gene
hub_gene_ids_cluster8_35 <- cluster8_hub_genes_35_GO$Gene
hub_gene_ids_cluster9_control <- cluster9_hub_genes_control_GO$Gene
hub_gene_ids_cluster9_30 <- cluster9_hub_genes_30_GO$Gene
hub_gene_ids_cluster9_35 <- cluster9_hub_genes_35_GO$Gene


# Write to selection.txt
writeLines(hub_gene_ids_cluster5_control, "selection_cluster5_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster5_30, "selection_cluster5_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster5_35, "selection_cluster5_hub_genes_35_GO.txt")
writeLines(hub_gene_ids_cluster8_control, "selection_cluster8_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster8_30, "selection_cluster8_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster8_35, "selection_cluster8_hub_genes_35_GO.txt")
writeLines(hub_gene_ids_cluster9_control, "selection_cluster9_hub_genes_control_GO.txt")
writeLines(hub_gene_ids_cluster9_30, "selection_cluster9_hub_genes_30_GO.txt")
writeLines(hub_gene_ids_cluster9_35, "selection_cluster9_hub_genes_35_GO.txt")

#load for Viseago
selection_cluster5_control<-scan(
    "selection_cluster5_hub_genes_control_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster5_30<-scan(
    "selection_cluster5_hub_genes_30_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster5_35<-scan(
    "selection_cluster5_hub_genes_35_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster8_control<-scan(
    "selection_cluster8_hub_genes_control_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster8_30<-scan(
    "selection_cluster8_hub_genes_30_GO.txt",
    quiet=TRUE,
    what=""
)

selection_cluster8_35<-scan(
    "selection_cluster8_hub_genes_35_GO.txt",
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


# Define gene lists for each cluster/temperature
gene_lists <- list(
  cluster5_control = selection_cluster5_control,
  cluster5_30      = selection_cluster5_30,
  cluster5_35      = selection_cluster5_35,
  cluster8_control = selection_cluster8_control,
  cluster8_30      = selection_cluster8_30,
  cluster8_35      = selection_cluster8_35,
  cluster9_control = selection_cluster9_control,
  cluster9_30      = selection_cluster9_30,
  cluster9_35      = selection_cluster9_35
)

library(ViSEAGO)

results_list <- list()

for (name in names(gene_lists)) {
  for (ont in c("BP", "MF")) {
    topgo <- ViSEAGO::create_topGOdata(
      geneSel = gene_lists[[name]],
      allGenes = background,
      gene2GO = myGENE2GO_Mcap,
      ont = ont,
      nodeSize = 5
    )
    res <- topGO::runTest(topgo, algorithm = "classic", statistic = "fisher")
    pvals <- topGO::score(res)
    # pvals_adj <- p.adjust(pvals, method = "BH")
    # sig_terms <- names(pvals_adj)[pvals_adj < 0.1]
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
    ss <- ViSEAGO::build_GO_SS(gene2GO = myGENE2GO_Mcap, enrich_GO_terms = enrich)
    ss <- ViSEAGO::compute_SS_distances(ss, distance = "Wang")
    heatmap <- ViSEAGO::GOterms_heatmap(ss, showIC = TRUE, showGOlabels = TRUE)
    # Only assign if everything succeeded
    results_list[[paste0(name, "_", ont)]] <- list(
      topgo = topgo,
      enrich = enrich,
      ss = ss,
      heatmap = heatmap
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


##### plot top 20 enriched GO terms, facet by temperature
library(dplyr)
library(ggplot2)
library(forcats)

# Extract enrichment info from ViSEAGO results
compare_df <- bind_rows(
  lapply(valid_results, function(res) {
    dt_obj <- ViSEAGO::show_table(res$heatmap)
    df <- as.data.frame(dt_obj$x$data)
    df$ID <- sub(".*GO:(\\d+).*", "GO:\\1", df$`GO ID`)
    df
  }),
  .id = "cluster_condition"
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
  scale_color_gradient(low = "lightgreen", high = "darkgreen") +
  labs(
    title = "Top 20 Enriched GO Terms per Cluster",
    size = "Gene Count",
    color = "-log10(p-value)"
  ) +
  scale_x_discrete(expand = c(0, 1)) + # less space on x-axis 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1)
  ) +
  facet_wrap(~ temperature, ncol = 1, scales = "free_y")

ggsave("ViSEAGO_top10_GOterms_faceted_dotplot_FDR.pdf", width = 10, height = 19)


# Save compare_df
# Add cluster column and clean cluster_condition in compare_df
compare_df <- compare_df %>%
  mutate(
    cluster = sub("(_control|_30|_35)$", "", cluster_condition)
  )

# Add direction column to compare_df
compare_df <- compare_df %>%
  mutate(
    direction = case_when(
      cluster %in% c("cluster8", "cluster9") ~ "up",
      cluster %in% c("cluster5") ~ "down",
      TRUE ~ NA_character_
    )
  )

write.csv(compare_df, "Sign_GO_p01Fisher_Mcap.csv", row.names = FALSE)

# Save top_terms
# Add direction column to top_terms
top_terms <- top_terms %>%
  mutate(
    direction = case_when(
      cluster %in% c("cluster8", "cluster9") ~ "up",
      cluster %in% c("cluster5") ~ "down",
      TRUE ~ NA_character_
    )
  )

write.csv(top_terms, "Top10_GO_p01Fisher_Mcap.csv", row.names = FALSE)





###### GO Slim Analysis
library(GSEABase)

# Download GO slim OBO file and read as a GeneSetCollection
goslim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo")
#format-version: 1.2

# Extract all enriched GO terms from your results
enriched_go_terms <- unique(unlist(
  lapply(valid_results, function(res) {
    dt_obj <- ViSEAGO::show_table(res$heatmap)
    df <- as.data.frame(dt_obj$x$data)
    sub(".*GO:(\\d+).*", "GO:\\1", df$`GO ID`)
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
    dt_obj <- ViSEAGO::show_table(res$heatmap)
    df <- as.data.frame(dt_obj$x$data)
    # Clean GO IDs
    df$ID <- sub(".*GO:(\\d+).*", "GO:\\1", df$`GO ID`)
    df
  }),
  .id = "cluster_condition"
)

GO.BP <- right_join(BPslim, compare_df, by="ID") #add back GO enrichment info for each offspring term
GO.MF <- right_join(MFslim, compare_df, by="ID") #add back GO enrichment info for each offspring term

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
  scale_x_discrete(expand = c(0, 0.9)) + # less space on x-axis
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
    plot.margin = margin(2, 2, 2, 2) # less margin
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

pdf("GO_slim_faceted_dotplot_Mcap.pdf", width = 10, height = 11)
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
      cluster %in% c("cluster8", "cluster9") ~ "up",
      cluster %in% c("cluster5") ~ "down",
      TRUE ~ NA_character_
    )
  )

write.csv(GO.slim_summary, "GOslim_All_p01Fisher_Mcap.csv", row.names = FALSE)
