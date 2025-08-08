
############## Functional Annotation - Mcap
## Compilation of the output of different methods: Uniprot, Interproscan and EggNog

# Load libraries
library("tidyverse") 
library("dplyr")
library("readr")
library("tidyr")

#Diamond results

blast <- read_tsv("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Mcap/Mcap.diamondBlastpNCBInr", col_names = FALSE)
colnames(blast) <- c("seqName", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(blast)
dim(blast)
#[1] 49746    12


#Unitprot results

#8 out of 29,916 UniProtKB AC/ID identifiers were successfully mapped to UniProtKB IDs
uni1 <- read_tsv("Uniprot_UniProtKB_to_UniProtKB.tsv", col_names = TRUE)
uni1 <- uni1[,c(1,2,4,5:8,11,13)]
colnames(uni1) <- c("from", "top_hit", "uniprotkb_entry", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology")
head(uni1)

#281 out of 29,916 RefSeq Protein identifiers were successfully mapped to UniProtKB IDs 
uni2 <- read_tsv("Uniprot_RefSeqProtein_to_UniprotKB.tsv", col_names = TRUE)
uni2 <- uni2[,c(1,2,4,5:8,11,13)]
colnames(uni2) <- c("from", "top_hit", "uniprotkb_entry", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology")
head(uni2)

#7,922 out of 29,916 EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to UniProtKB
uni3 <- read_tsv("Uniprot_EMBLGenBankDDBJ_to_UniProtKB.tsv", col_names = TRUE)
uni3 <- uni3[,c(1,2,4,5:8,11,13)]
colnames(uni3) <- c("from", "top_hit", "uniprotkb_entry", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology")
head(uni3)

#506 out of 29,916 Ensembl Genomes Protein identifiers were successfully mapped to UniProtKB IDs 
uni4 <- read_tsv("Uniprot_EnsemblGenomesProtein_to_UniProtKB.tsv", col_names = TRUE)
uni4 <- uni4[,c(1,2,4,5:8,11,13)]
colnames(uni4) <- c("from", "top_hit", "uniprotkb_entry", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology")
head(uni4)

#compile results
Uniprot_results2 <- bind_rows(uni1, uni2, uni3, uni4)
# Uniprot_results2 <- unique(Uniprot_results)
# Uniprot_results$go_ids <- gsub(" ", "", Uniprot_results$go_ids)

Uniprot_results2 <- Uniprot_results2 %>% filter(grepl("GO",go_ids)) 

## Add gene name from Diamond results
# Clean top_hit to remove version suffix (e.g., .1, .2)
blast <- blast %>%
  mutate(top_hit_clean = sub("\\..*", "", top_hit))

Uniprot_results2 <- Uniprot_results2 %>%
  mutate(top_hit_clean = sub("\\..*", "", from))

# Now join using the cleaned top_hit
Uniprot_final2 <- Uniprot_results2 %>%
  left_join(blast, by = "top_hit_clean")
  
#This keeps the first row for each top_hit_clean, along with all columns.
Uniprot_results_unique2 <- Uniprot_final2  %>%
  distinct(top_hit_clean, .keep_all = TRUE)

Uniprot_results_unique2 <- Uniprot_results_unique2[,c(1:4,8:11)]  #4951 Genes with GO terms 
colnames(Uniprot_results_unique2)[8] <- "protein_accession"

##### InterpoScan
# Load InterProScan results
interpro <- read_tsv("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Mcap/mcap.interpro.tsv", col_names = FALSE)

#remove 15th column
interpro <- interpro[, -15]

# Step 2: Assign correct column names manually
colnames(interpro) <- c(
  "protein_accession", "sequence_md5", "sequence_length", "analysis",
  "signature_accession", "signature_description", "start_location", "stop_location",
  "score", "status", "date", "interpro_accession", "interpro_description", "go_ids"
)

#retain only rows with GO terms
interpro <- interpro %>% filter(grepl("GO", go_ids))

#remove rows with "score" above 0.00001
interpro <- interpro %>%
  mutate(score = as.numeric(score))

interpro <- interpro %>% filter(score < 0.00001)

#remove rows were "score" is - and were the scire is exaclty 0
interpro <- interpro %>% filter(score != "0")
interpro <- interpro %>% filter(score != "-")  

#for rows with the same protein_accession keep the one with the lowest score
interpro_final <- interpro %>%
  group_by(protein_accession) %>%
  slice(which.min(score)) %>%
  ungroup()     ### 8881 genes with Go annotations

#keep neeeded columns
interpro_final <- interpro_final[,c(1,6,13,14)]


## EggNog

# Load InterProScan results
eggnog <- read_tsv("/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes.EggNog_results.txt", col_names = TRUE)

#rename needed columns
colnames(eggnog)[1] <- "protein_accession"
colnames(eggnog)[10] <- "go_ids"

#retain only rows with GO terms
eggnog <- eggnog %>% filter(grepl("GO", go_ids))  ### 9588 genes with Go annotations

#keep neeeded columns
eggnog_final  <- eggnog [,c(1,5,8,10,12,13,17)]


### Find unique and overlapping GO terms

# add a column to each dataset indicate the source of each entry
Uniprot_results_unique2 <- Uniprot_results_unique2 %>%
  mutate(source = "Uniprot")

interpro_final <- interpro_final %>%
    mutate(source = "InterProScan")

eggnog_final <- eggnog_final %>%
    mutate(source = "EggNog")

#Generate lists of GO terms for each method
Uniprot_GO <- Uniprot_results_unique2 %>% dplyr::select(protein_accession, go_ids, source)
Interpro_GO <- interpro_final %>% dplyr::select(protein_accession, go_ids, source)
Eggnog_GO <- eggnog_final %>% dplyr::select(protein_accession, go_ids, source)

# Standardize separators
Uniprot_GO <- Uniprot_GO %>%
  mutate(go_ids = gsub(";", ";", go_ids)) # already ;
Interpro_GO <- Interpro_GO %>%
  mutate(go_ids = gsub("\\|", ";", go_ids)) # replace | with ;
Eggnog_GO <- Eggnog_GO %>%
  mutate(go_ids = gsub(",", ";", go_ids)) # replace , with ;

# Split GO terms into individual rows
Uniprot_long <- Uniprot_GO %>%
  separate_rows(go_ids, sep = ";") %>%
  mutate(go_ids = trimws(go_ids))
Interpro_long <- Interpro_GO %>%
  separate_rows(go_ids, sep = ";") %>%
  mutate(go_ids = trimws(go_ids))
Eggnog_long <- Eggnog_GO %>%
  separate_rows(go_ids, sep = ";") %>%
  mutate(go_ids = trimws(go_ids))

# Combine all datasets
all_GO <- bind_rows(Uniprot_long, Interpro_long, Eggnog_long)

# Summarize: for each protein_accession, get unique and overlapping GOs
GO_summary <- all_GO %>%
  group_by(protein_accession, go_ids) %>%
  summarise(sources = paste(sort(unique(source)), collapse = ","), .groups = "drop")

# Find GOs present in all three sources for each protein
GO_in_all <- GO_summary %>%
  filter(grepl("Uniprot", sources) & grepl("InterProScan", sources) & grepl("EggNog", sources))
## 487

# Find GOs unique to one source for each protein
GO_unique <- GO_summary %>%
  group_by(protein_accession, go_ids) %>%
  filter(length(strsplit(sources, ",")[[1]]) == 1) %>%
  ungroup()
### 1361923

# Collapse all unique GO terms for each protein_accession
final_GO_table_Mcap <- all_GO %>%
  group_by(protein_accession) %>%
  summarise(go_ids = paste(sort(unique(go_ids)), collapse = ";")) %>%
  ungroup()
### 17220

