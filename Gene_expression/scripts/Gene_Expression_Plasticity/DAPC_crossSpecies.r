

############# Gene expression plasticity - Discriminant Analysis of Principal Components #############

## load libraries
library("DESeq2")               #BiocManager::install("DESeq2")
library("factoextra")           #install.packages("factoextra")
library("tidyverse")            
library("RColorBrewer")
library("ggplot2")              
library("patchwork")            #install.packages("patchwork")
library("dplyr")
library(adegenet)       #install.packages("adegenet")

#treatment information
treatmentinfo <- read.csv("RNAseq_data_all.csv", header = TRUE, sep = ";")
#keep only control, 30 and 35°C
treatmentinfo <- treatmentinfo[treatmentinfo$temp %in% c(26.8, 30, 35), ]

#gene count matrix - Broccoli orthologs
logTPM_merged <- as.data.frame(read.csv("logTPM_merged.csv"))
#keep only control, 30 and 35°C
logTPM_merged_filt <- logTPM_merged[, grep("10|7|9", colnames(logTPM_merged))]

#remove samples below 5 million reads
logTPM_merged_filt <- as.data.frame(logTPM_merged_filt)
logTPM_merged <- logTPM_merged_filt %>% dplyr::select(-pacu_F10)
treatmentinfo <- treatmentinfo[!(treatmentinfo$sample_id %in% c("pacu_B1", "pacu_F10")), ]
logTPM_merged <- logTPM_merged_filt %>% dplyr::select(-pcom_B7)
treatmentinfo <- treatmentinfo[!(treatmentinfo$sample_id %in% c("pcom_B7")), ]

# Reorder treatmentinfo to match the columns of logTPM_merged_filt
treatmentinfo <- treatmentinfo[match(colnames(logTPM_merged_filt), treatmentinfo$sample_id), ]


######### Subset for 26.8 vs 30°C
# Subset metadata for 26.8 and 30°C
treatmentinfo_30 <- treatmentinfo[treatmentinfo$temp %in% c(26.8, 30), ]

# Subset expression matrix columns to match samples in treatmentinfo_30
logTPM_merged_filt_30 <- logTPM_merged[, colnames(logTPM_merged) %in% treatmentinfo_30$sample_id]

# Reorder metadata to match expression matrix columns
treatmentinfo_30 <- treatmentinfo_30[match(colnames(logTPM_merged_filt_30), treatmentinfo_30$sample_id), ]

# Transpose for DAPC
datExpr_30 <- t(logTPM_merged_filt_30)

#Assigning groups from themetadata, ensure treatmentinfo is ordered to match datExpr rownames
treatmentinfo_30 <- treatmentinfo_30[match(rownames(datExpr_30), treatmentinfo_30$sample_id), ]

# Create group variable
group_30 <- as.factor(paste(treatmentinfo_30$temp, treatmentinfo_30$species, sep = "_"))

# Create temp_status variable
treatmentinfo_30$temp_status <- ifelse(treatmentinfo_30$temp == 26.8, "control", "mid")
treatmentinfo_30$species <- as.factor(treatmentinfo_30$species)
treatmentinfo_30$temp_status <- as.factor(treatmentinfo_30$temp_status)

### Run DAPC

# This will suggest the optimal number of PCs to retain:
set.seed(123)
tmp <- dapc(datExpr_30, group_30, n.pca = 50, n.da = 1) # 1 is 2 temp - 1
a_score <- optim.a.score(tmp) #8

# Run DAPC with the optimal number of PCs
dapc_res <- dapc(datExpr_30, group_30)
#The number of discriminant functions (DFs) you can retain is at most (number of groups - 1).
#So, with 2 temperature groups, the maximum is 1.

## Plot DAPC results
# Basic scatter plot
scatter(dapc_res, posi.da="bottomright", bg="white", pch=20, cstar=0, col=c("blue", "orange", "green"))

## Extract DAPC scores for further modeling
dapc_scores_30 <- dapc_res$ind.coord  # Each row = sample, columns = discriminant functions
# You can now use dapc_scores in MCMCglmm or other modeling

# Combine DAPC scores with metadata
dapc_df <- cbind(treatmentinfo_30, dapc_scores_30)

# Plot DAPC results
# Assign colors to each group (make sure the order matches the group levels)
group_colors <- c("blue", "orange", "green")

# Plot with legend
scatter(
  dapc_res,
  grp = group_30,
  col = group_colors,
  posi.leg = "topright",
  bg = "white",
  pch = 20,
  cstar = 0
)

## plotting the DAPC scores as density plots for each group
library(ggplot2)

# Assign colors by species
species_colors <- c("mcap" = "green", "pacu" = "cyan", "pcom" = "orange")
dapc_df$species <- factor(dapc_df$species, levels = c("mcap", "pacu", "pcom"))

# Add a shading variable: TRUE for control, FALSE for high temp
dapc_df$shaded <- dapc_df$temp_status == "control"

# Make sure 'shaded' is a factor with clear labels
dapc_df$Condition <- ifelse(dapc_df$temp_status == "control", "Control (shaded)", "30°C (solid)")
dapc_df$Condition <- factor(dapc_df$Condition, levels = c("Control (shaded)", "30°C (solid)"))

ggplot(dapc_df, aes(x = LD1, fill = species, color = species, alpha = Condition)) +
  geom_density(position = "identity", adjust = 1.2) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  scale_alpha_manual(
    values = c("Control (shaded)" = 0.1, "30°C (solid)" = 1),
    name = "Condition"
  ) +
  theme_minimal() +
  labs(
    title = "DAPC LD1 Density by Species and Temperature",
    x = "DAPC LD1",
    y = "Density",
    fill = "Species",
    color = "Species"
  )


# Calculate mean LD1 for each species and temp_status
library(dplyr)

mean_LD1 <- dapc_df %>%
  group_by(species, temp_status) %>%
  summarise(mean_LD1 = mean(LD1), .groups = "drop")

# Spread to wide format for easy subtraction
mean_LD1_wide <- tidyr::pivot_wider(mean_LD1, names_from = temp_status, values_from = mean_LD1)

# Calculate the difference: mid temp - control
mean_LD1_wide$diff_mid_vs_control <- mean_LD1_wide$mid - mean_LD1_wide$control

print(mean_LD1_wide)
# # A tibble: 3 × 4
#   species control    mid diff_mid_vs_control
#   <fct>     <dbl>  <dbl>               <dbl>
# 1 mcap       23.0  21.7                -1.33
# 2 pacu      -35.9 -38.7                -2.87
# 3 pcom       12.3   8.62               -3.68


## Add arrows and labels to the ggplot:

# the density plot
p <- ggplot(dapc_df, aes(x = LD1, fill = species, color = species, alpha = Condition)) +
  geom_density(position = "identity", adjust = 1.5) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  scale_alpha_manual(
    values = c("Control (shaded)" = 0.1, "30°C (solid)" = 1),
    name = "Condition"
  ) +
  theme_minimal() +
  labs(
    title = "DAPC LD1 Density by Species and Temperature",
    x = "DAPC LD1",
    y = "Density",
    fill = "Species",
    color = "Species"
  )

# Add arrows and labels for each species
p +
  xlim(min(dapc_df$LD1) - 2, max(dapc_df$LD1) + 2) +
  geom_segment(
    data = mean_LD1_wide,
    aes(
      x = control, xend = mid,
      y = 0.05, yend = 0.05
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    size = 0.7,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = mean_LD1_wide,
    aes(
      x = (control + mid) / 2,
      y = 0.06,
      label = round(diff_mid_vs_control, 2)
    ),
    size = 4,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  theme(panel.grid = element_blank())




### MCMCglmm analysis

install.packages("MCMCglmm")
library(MCMCglmm)

#Fit the MCMCglmm model
#Model the first discriminant function (DF1):

# Fit the model: DAPC score ~ tmain effect of temp + main effect of species + interaction
model <- MCMCglmm(
  LD1 ~ temp_status * species,
  data = dapc_df,
  family = "gaussian",
  nitt = 13000, burnin = 3000, thin = 10
)

#Extract posterior samples and calculate differences
post_samples <- model$Sol

# extracts the species-specific effects of being at higher temperature
effect_Mcap_mid <- post_samples[, "temp_statusmid"] #"mcap" is the reference species (the baseline in the model).
effect_Pacu_mid <- post_samples[, "temp_statusmid"] + post_samples[, "temp_statusmid:speciespacu"]
effect_Pcom_mid <- post_samples[, "temp_statusmid"] + post_samples[, "temp_statusmid:speciespcom"]

# Calculate differences in absolute values (plasticity difference)
diff_Mcap_vs_Pacu <- abs(effect_Mcap_mid) - abs(effect_Pacu_mid)
diff_Mcap_vs_Pcom <- abs(effect_Mcap_mid) - abs(effect_Pcom_mid)
diff_Pacu_vs_Pcom <- abs(effect_Pacu_mid) - abs(effect_Pcom_mid)

# Calculate MCMC-based posterior probability values
p_Mcap_vs_Pacu <- mean(diff_Mcap_vs_Pacu > 0) ## 0.09
p_Mcap_vs_Pcom <- mean(diff_Mcap_vs_Pcom > 0) ## 0.02
p_Pacu_vs_Pcom <- mean(diff_Pacu_vs_Pcom > 0) ## 0.25




########### Subset for 26.8 vs 35°C
# Subset metadata for 26.8 and 35°C
treatmentinfo_35 <- treatmentinfo[treatmentinfo$temp %in% c(26.8, 35), ]

# Subset expression matrix columns to match samples in treatmentinfo_35
logTPM_merged_filt_35 <- logTPM_merged[, colnames(logTPM_merged) %in% treatmentinfo_35$sample_id]

# Reorder metadata to match expression matrix columns
treatmentinfo_35 <- treatmentinfo_35[match(colnames(logTPM_merged_filt_35), treatmentinfo_35$sample_id), ]

# Transpose for DAPC
datExpr_35 <- t(logTPM_merged_filt_35)
# grouping factor: create vector group that assigns each sample to its group
#Assigning groups from themetadata, ensure treatmentinfo is ordered to match datExpr rownames
treatmentinfo_35 <- treatmentinfo_35[match(rownames(datExpr_35), treatmentinfo_35$sample_id), ]
group_35 <- as.factor(paste(treatmentinfo_35$temp, treatmentinfo_35$species, sep = "_"))

# Create temp_status variable
treatmentinfo_35$temp_status <- ifelse(treatmentinfo_35$temp == 26.8, "control", "high")
treatmentinfo_35$species <- as.factor(treatmentinfo_35$species)
treatmentinfo_35$temp_status <- as.factor(treatmentinfo_35$temp_status)

### Run DAPC

# This will suggest the optimal number of PCs to retain:
set.seed(123)
tmp <- dapc(datExpr_35, group_35, n.pca = 50, n.da = 1) # 1 is 2 temp - 1
a_score <- optim.a.score(tmp) #8

# Run DAPC with the optimal number of PCs
dapc_res <- dapc(datExpr_35, group_35)
#The number of discriminant functions (DFs) you can retain is at most (number of groups - 1).
#So, with 2 temperature groups, the maximum is 1.

## Plot DAPC results
# Basic scatter plot
scatter(dapc_res, posi.da="bottomright", bg="white", pch=20, cstar=0, col=c("blue", "orange", "green"))

## Extract DAPC scores for further modeling
dapc_scores_35 <- dapc_res$ind.coord  # Each row = sample, columns = discriminant functions
# You can now use dapc_scores in MCMCglmm or other modeling

# Combine DAPC scores with metadata
dapc_df <- cbind(treatmentinfo_35, dapc_scores_35)

# Plot DAPC results
# Assign colors to each group (make sure the order matches the group levels)
group_colors <- c("blue", "orange", "green")

# Plot with legend
scatter(
  dapc_res,
  grp = group_35,
  col = group_colors,
  posi.leg = "topright",
  bg = "white",
  pch = 20,
  cstar = 0
)

## plotting the DAPC scores as density plots for each group
library(ggplot2)

# Assign colors by species
species_colors <- c("mcap" = "green", "pacu" = "cyan", "pcom" = "orange")
dapc_df$species <- factor(dapc_df$species, levels = c("mcap", "pacu", "pcom"))

# Add a shading variable: TRUE for control, FALSE for high temp
dapc_df$shaded <- dapc_df$temp_status == "control"

# Make sure 'shaded' is a factor with clear labels
dapc_df$Condition <- ifelse(dapc_df$temp_status == "control", "Control (shaded)", "35°C (solid)")
dapc_df$Condition <- factor(dapc_df$Condition, levels = c("Control (shaded)", "35°C (solid)"))

ggplot(dapc_df, aes(x = LD1, fill = species, color = species, alpha = Condition)) +
  geom_density(position = "identity", adjust = 1.2) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  scale_alpha_manual(
    values = c("Control (shaded)" = 0.1, "35°C (solid)" = 1),
    name = "Condition"
  ) +
  theme_minimal() +
  labs(
    title = "DAPC LD1 Density by Species and Temperature",
    x = "DAPC LD1",
    y = "Density",
    fill = "Species",
    color = "Species"
  )


# Calculate mean LD1 for each species and temp_status
library(dplyr)

mean_LD1 <- dapc_df %>%
  group_by(species, temp_status) %>%
  summarise(mean_LD1 = mean(LD1), .groups = "drop")

# Spread to wide format for easy subtraction
mean_LD1_wide <- tidyr::pivot_wider(mean_LD1, names_from = temp_status, values_from = mean_LD1)

# Calculate the difference: high temp - control
mean_LD1_wide$diff_high_vs_control <- mean_LD1_wide$high - mean_LD1_wide$control

print(mean_LD1_wide)
# A tibble: 3 × 4
#   species control  high diff_high_vs_control
#   <fct>     <dbl> <dbl>                <dbl>
# 1 mcap       5.17  10.6                 5.39
# 2 pacu     -14.2  -45.0               -30.8 
# 3 pcom      21.7   24.3                 2.56


## Add arrows and labels to the ggplot:

# the density plot
p <- ggplot(dapc_df, aes(x = LD1, fill = species, color = species, alpha = Condition)) +
  geom_density(position = "identity", adjust = 1.5) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  scale_alpha_manual(
    values = c("Control (shaded)" = 0.1, "35°C (solid)" = 1),
    name = "Condition"
  ) +
  theme_minimal() +
  labs(
    title = "DAPC LD1 Density by Species and Temperature",
    x = "DAPC LD1",
    y = "Density",
    fill = "Species",
    color = "Species"
  )


# Add arrows and labels for each species
p +
  xlim(min(dapc_df$LD1) - 2, max(dapc_df$LD1) + 2) +
  geom_segment(
    data = mean_LD1_wide,
    aes(
      x = control, xend = high,
      y = 0.05, yend = 0.05
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    size = 0.7,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = mean_LD1_wide,
    aes(
      x = (control + high) / 2,
      y = 0.06,
      label = round(diff_high_vs_control, 2)
    ),
    size = 4,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  theme(panel.grid = element_blank())




### MCMCglmm analysis

install.packages("MCMCglmm")
library(MCMCglmm)

#Fit the MCMCglmm model
#Model the first discriminant function (DF1):

# Fit the model: DAPC score ~ tmain effect of temp + main effect of species + interaction
model <- MCMCglmm(
  LD1 ~ temp_status * species,
  data = dapc_df,
  family = "gaussian",
  nitt = 13000, burnin = 3000, thin = 10
)

#Extract posterior samples and calculate differences
post_samples <- model$Sol

# extracts the species-specific effects of being at higher temperature
effect_Mcap_high <- post_samples[, "temp_statushigh"] #"mcap" is the reference species (the baseline in the model).
effect_Pacu_high <- post_samples[, "temp_statushigh"] + post_samples[, "temp_statushigh:speciespacu"]
effect_Pcom_high <- post_samples[, "temp_statushigh"] + post_samples[, "temp_statushigh:speciespcom"]

# Calculate differences in absolute values (plasticity difference)
diff_Mcap_vs_Pacu <- abs(effect_Mcap_high) - abs(effect_Pacu_high)
diff_Mcap_vs_Pcom <- abs(effect_Mcap_high) - abs(effect_Pcom_high)
diff_Pacu_vs_Pcom <- abs(effect_Pacu_high) - abs(effect_Pcom_high)

# Calculate MCMC-based posterior probability values
p_Mcap_vs_Pacu <- mean(diff_Mcap_vs_Pacu > 0) ## 0
p_Mcap_vs_Pcom <- mean(diff_Mcap_vs_Pcom > 0) ## 0.9
p_Pacu_vs_Pcom <- mean(diff_Pacu_vs_Pcom > 0) ## 1

#indicate the probability that the absolute plasticity effect (the effect of being at high temperature) is greater in one species than another.
#p_Mcap_vs_Pcom = 0.9 means that in 90% of the posterior samples, the absolute effect for Mcap is greater than for Pcom.
#p_Mcap_vs_Pacu = 0 means that in 0% of the posterior samples, the absolute effect for Mcap is greater than for Pacu (i.e., Pacu almost always has a greater effect than Mcap).
#A value close to 1 (e.g., 0.9) suggests strong evidence that the first species (Mcap) has a greater plasticity effect than the second (Pcom).
#A value close to 0 suggests strong evidence that the second species (Pacu) has a greater effect than the first (Mcap).
#A value near 0.5 would indicate no difference.
#There is strong evidence that Pacu has greater plasticity than Mcap (p_Mcap_vs_Pacu = 0).
#There is moderate evidence that Mcap has greater plasticity than Pcom (p_Mcap_vs_Pcom = 0.9).
#There is strong evidence that Pacu has greater plasticity than Pcom (p_Pacu_vs_Pcom = 1).



############## Test for differences in gene expression at the control temperature using limma

# Subset metadata and expression matrix for control samples only
control_samples <- treatmentinfo[treatmentinfo$temp == 26.8, ]
logTPM_control <- logTPM_merged_filt[, colnames(logTPM_merged_filt) %in% control_samples$sample_id]

# Reorder metadata to match columns of expression matrix
control_samples <- control_samples[match(colnames(logTPM_control), control_samples$sample_id), ]

# Prepare design matrix for limma
library(limma)

# Make sure species is a factor
control_samples$species <- factor(control_samples$species)

## Design matrix for species effect
design <- model.matrix(~ 0 + species, data = control_samples)
# By default, model.matrix(~ species, ...) would include an intercept (baseline) and then coefficients for all but one species (the reference).
# Using ~ 0 + species removes the intercept and creates a separate column for each species in the design matrix.
# This allows you to directly compare each species’ mean expression
colnames(design) <- levels(control_samples$species)

### Run limma
# Fit the linear model
fit <- lmFit(logTPM_control, design)

# Make pairwise contrasts between species (mcap vs pacu, mcap vs pcom, pacu vs pcom)
contrast.matrix <- makeContrasts(
  mcap_vs_pacu = mcap - pacu,
  mcap_vs_pcom = mcap - pcom,
  pacu_vs_pcom = pacu - pcom,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get differentially expressed genes for each contrast
DE_mcap_vs_pacu <- topTable(fit2, coef = "mcap_vs_pacu", number = Inf)
DE_mcap_vs_pcom <- topTable(fit2, coef = "mcap_vs_pcom", number = Inf)
DE_pacu_vs_pcom <- topTable(fit2, coef = "pacu_vs_pcom", number = Inf)
# logFC: The log2 fold change between the two groups being compared (e.g., mcap vs pacu).
# A value of 1 means a 2-fold increase in the first group compared to the second.
# AveExpr: The average (mean) log2 expression level of the gene across all samples in the comparison.
# t: The moderated t-statistic for the test of differential expression.
# P.Value: The raw p-value for the test of differential expression (before multiple testing correction).
# adj.P.Val: The adjusted p-value (Benjamini-Hochberg FDR) for multiple testing correction.
# B: The log-odds that the gene is differentially expressed (higher B means more likely to be truly differentially expressed).


### Run PCA
# Subset metadata and expression matrix for control samples only
control_samples <- treatmentinfo[treatmentinfo$temp == 26.8, ]
expr_for_pca <- logTPM_merged_filt[, colnames(logTPM_merged_filt) %in% control_samples$sample_id]

# Transpose so samples are rows, genes are columns
expr_for_pca_t <- t(expr_for_pca)

# Run PCA
pca_res <- prcomp(expr_for_pca_t, scale. = TRUE)

# Prepare data for plotting
pca_df <- as.data.frame(pca_res$x)
pca_df$species <- control_samples$species
pca_df$temp <- control_samples$temp

# Plot PCA (PC1 vs PC2, colored by species)
library(ggplot2)

species_colors <- c("mcap" = "green", "pacu" = "cyan", "pcom" = "orange")
ggplot(pca_df, aes(x = PC1, y = PC2, color = species, shape = as.factor(temp))) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_df), color = species), vjust = -1, size = 3) +  # Color labels by species
  scale_color_manual(values = species_colors) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "PCA of Gene Expression - control temperature", color = "Species", shape = "Temperature")




### Create a mapping of ortholog IDs to gene names for each species
# join ortho_map with each DE table by first converting the rownames of the DE table to a column (e.g., orthoID), then using dplyr::left_join().

ortho_map <- data.frame(
  orthoID = rownames(logTPM_merged),
  mcap_gene = genes_mcap,
  pacu_gene = genes_pacu,
  pcom_gene = genes_pcom
) #### this is the output from Broccoli

# For mcap_vs_pacu
DE_mcap_vs_pacu$orthoID <- rownames(DE_mcap_vs_pacu)
DE_mcap_vs_pacu_annot <- left_join(ortho_map, DE_mcap_vs_pacu, by = "orthoID")

# For mcap_vs_pcom
DE_mcap_vs_pcom$orthoID <- rownames(DE_mcap_vs_pcom)
DE_mcap_vs_pcom_annot <- left_join(ortho_map, DE_mcap_vs_pcom, by = "orthoID")

# For pacu_vs_pcom
DE_pacu_vs_pcom$orthoID <- rownames(DE_pacu_vs_pcom)
DE_pacu_vs_pcom_annot <- left_join(ortho_map, DE_pacu_vs_pcom, by = "orthoID")


###### search for top 10% hub genes from WGCNA

##load gene clusters - Mcap

#cluster5
cluster5_hub_gene_counts_control_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_control_cluster5_2.csv")
cluster5_hub_gene_counts_30_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat30_cluster5_2.csv")
cluster5_hub_gene_counts_35_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat35_cluster5_2.csv")

#cluster8
cluster8_hub_gene_counts_control_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_control_cluster8_2.csv")
cluster8_hub_gene_counts_30_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat30_cluster8_2.csv")
cluster8_hub_gene_counts_35_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat35_cluster8_2.csv")

#cluster9
cluster9_hub_gene_counts_control_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_control_cluster9_2.csv")
cluster9_hub_gene_counts_30_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat30_cluster9_2.csv")
cluster9_hub_gene_counts_35_Mcap <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/top10pct_hub_genes_treat35_cluster9_2.csv")

### keep only unique hub genes per each temperature and cluster

# create vectors of gene IDs
cluster5_hub_genes_control_Mcap <- cluster5_hub_gene_counts_control_Mcap$Gene
cluster5_hub_genes_30_Mcap <- cluster5_hub_gene_counts_30_Mcap$Gene
cluster5_hub_genes_35_Mcap <- cluster5_hub_gene_counts_35_Mcap$Gene

cluster8_hub_genes_control_Mcap <- cluster8_hub_gene_counts_control_Mcap$Gene
cluster8_hub_genes_30_Mcap <- cluster8_hub_gene_counts_30_Mcap$Gene
cluster8_hub_genes_35_Mcap <- cluster8_hub_gene_counts_35_Mcap$Gene

cluster9_hub_genes_control_Mcap <- cluster9_hub_gene_counts_control_Mcap$Gene
cluster9_hub_genes_30_Mcap <- cluster9_hub_gene_counts_30_Mcap$Gene
cluster9_hub_genes_35_Mcap <- cluster9_hub_gene_counts_35_Mcap$Gene

# Get unique hub genes for each temp 30 (non-shared genes)
cluster5_hub_genes_30_unique_Mcap <- setdiff(
  cluster5_hub_genes_30_Mcap,
  union(cluster5_hub_genes_control_Mcap, cluster5_hub_genes_35_Mcap)
)

cluster5_hub_genes_35_unique_Mcap <- setdiff(
  cluster5_hub_genes_35_Mcap,
  union(cluster5_hub_genes_control_Mcap, cluster5_hub_genes_30_Mcap)
)

cluster5_hub_genes_control_unique_Mcap <- setdiff(
  cluster5_hub_genes_control_Mcap,
  union(cluster5_hub_genes_30_Mcap, cluster5_hub_genes_35_Mcap)
)

cluster8_hub_genes_30_unique_Mcap <- setdiff(
  cluster8_hub_genes_30_Mcap,
  union(cluster8_hub_genes_control_Mcap, cluster8_hub_genes_35_Mcap)
)

cluster8_hub_genes_35_unique_Mcap <- setdiff(
  cluster8_hub_genes_35_Mcap,
  union(cluster8_hub_genes_control_Mcap, cluster8_hub_genes_30_Mcap)
)

cluster8_hub_genes_control_unique_Mcap <- setdiff(
  cluster8_hub_genes_control_Mcap,
  union(cluster8_hub_genes_30_Mcap, cluster8_hub_genes_35_Mcap)
)

cluster9_hub_genes_30_unique_Mcap <- setdiff(
  cluster9_hub_genes_30_Mcap,
  union(cluster9_hub_genes_control_Mcap, cluster9_hub_genes_35_Mcap)
)

cluster9_hub_genes_35_unique_Mcap <- setdiff(
  cluster9_hub_genes_35_Mcap,
  union(cluster9_hub_genes_control_Mcap, cluster9_hub_genes_30_Mcap)
)

cluster9_hub_genes_control_unique_Mcap <- setdiff(
  cluster9_hub_genes_control_Mcap,
  union(cluster9_hub_genes_30_Mcap, cluster9_hub_genes_35_Mcap)
)

# Get shared hub genesacross the 3 temperatures

cluster5_hub_genes_shared_Mcap <- Reduce(intersect, list(
  cluster5_hub_genes_control_Mcap,
  cluster5_hub_genes_30_Mcap,
  cluster5_hub_genes_35_Mcap
))

cluster8_hub_genes_shared_Mcap <- Reduce(intersect, list(
  cluster8_hub_genes_control_Mcap,
  cluster8_hub_genes_30_Mcap,
  cluster8_hub_genes_35_Mcap
))

cluster9_hub_genes_shared_Mcap <- Reduce(intersect, list(
  cluster9_hub_genes_control_Mcap,
  cluster9_hub_genes_30_Mcap,
  cluster9_hub_genes_35_Mcap
))

#Load diamond results
blast <- read_tsv("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Mcap/Mcap.diamondBlastpNCBInr", col_names = FALSE)
colnames(blast) <- c("Gene", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(blast)
dim(blast)
#[1] 49746    12

### find top_hit from blast output corresponding to the shared and unique hub genes per each cluster
# Create a data frame from the gene list
shared_genes_cluster5_Mcap <- data.frame(Gene = cluster5_hub_genes_shared_Mcap, stringsAsFactors = FALSE)
shared_genes_cluster8_Mcap <- data.frame(Gene = cluster8_hub_genes_shared_Mcap, stringsAsFactors = FALSE)
shared_genes_cluster9_Mcap <- data.frame(Gene = cluster9_hub_genes_shared_Mcap, stringsAsFactors = FALSE)

unique_cluster5_hub_genes_30_Mcap <- data.frame(Gene = cluster5_hub_genes_30_unique_Mcap, stringsAsFactors = FALSE)
unique_cluster5_hub_genes_35_Mcap <- data.frame(Gene = cluster5_hub_genes_35_unique_Mcap, stringsAsFactors = FALSE)
unique_cluster5_hub_genes_control_Mcap <- data.frame(Gene = cluster5_hub_genes_control_unique_Mcap, stringsAsFactors = FALSE)

unique_cluster8_hub_genes_30_Mcap <- data.frame(Gene = cluster8_hub_genes_30_unique_Mcap, stringsAsFactors = FALSE)
unique_cluster8_hub_genes_35_Mcap <- data.frame(Gene = cluster8_hub_genes_35_unique_Mcap, stringsAsFactors = FALSE)
unique_cluster8_hub_genes_control_Mcap <- data.frame(Gene = cluster8_hub_genes_control_unique_Mcap, stringsAsFactors = FALSE)

unique_cluster9_hub_genes_30_Mcap <- data.frame(Gene = cluster9_hub_genes_30_unique_Mcap, stringsAsFactors = FALSE)
unique_cluster9_hub_genes_35_Mcap <- data.frame(Gene = cluster9_hub_genes_35_unique_Mcap, stringsAsFactors = FALSE)
unique_cluster9_hub_genes_control_Mcap <- data.frame(Gene = cluster9_hub_genes_control_unique_Mcap, stringsAsFactors = FALSE)

# Inner join to get matching rows
shared_genes_blast_cluster5_Mcap  <- inner_join(shared_genes_cluster5_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster8_Mcap  <- inner_join(shared_genes_cluster8_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster9_Mcap  <- inner_join(shared_genes_cluster9_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster5_30_Mcap  <- inner_join(unique_cluster5_hub_genes_30_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster5_35_Mcap  <- inner_join(unique_cluster5_hub_genes_35_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster5_control_Mcap  <- inner_join(unique_cluster5_hub_genes_control_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster8_30_Mcap  <- inner_join(unique_cluster8_hub_genes_30_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster8_35_Mcap  <- inner_join(unique_cluster8_hub_genes_35_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster8_control_Mcap  <- inner_join(unique_cluster8_hub_genes_control_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster9_30_Mcap  <- inner_join(unique_cluster9_hub_genes_30_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster9_35_Mcap  <- inner_join(unique_cluster9_hub_genes_35_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_cluster9_control_Mcap  <- inner_join(unique_cluster9_hub_genes_control_Mcap, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)


#search proteins in NCBI
library(rentrez)

# Get the vector of accession numbers
accessions_shared_cluster5_Mcap <- shared_genes_blast_cluster5_Mcap$top_hit
accessions_shared_cluster8_Mcap <- shared_genes_blast_cluster8_Mcap$top_hit
accessions_shared_cluster9_Mcap <- shared_genes_blast_cluster9_Mcap$top_hit

accessions_unique_cluster5_30_Mcap <- unique_genes_blast_cluster5_30_Mcap$top_hit
accessions_unique_cluster5_35_Mcap <- unique_genes_blast_cluster5_35_Mcap$top_hit
accessions_unique_cluster5_control_Mcap <- unique_genes_blast_cluster5_control_Mcap$top_hit

accessions_unique_cluster8_30_Mcap <- unique_genes_blast_cluster8_30_Mcap$top_hit
accessions_unique_cluster8_35_Mcap <- unique_genes_blast_cluster8_35_Mcap$top_hit
accessions_unique_cluster8_control_Mcap <- unique_genes_blast_cluster8_control_Mcap$top_hit

accessions_unique_cluster9_30_Mcap <- unique_genes_blast_cluster9_30_Mcap$top_hit
accessions_unique_cluster9_35_Mcap <- unique_genes_blast_cluster9_35_Mcap$top_hit
accessions_unique_cluster9_control_Mcap <- unique_genes_blast_cluster9_control_Mcap$top_hit

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
results_shared_cluster5_Mcap  <- lapply(accessions_shared_cluster5_Mcap, fetch_protein_info)
final_shared_cluster5_Mcap <- do.call(rbind, results_shared_cluster5_Mcap)

results_shared_cluster8_Mcap  <- lapply(accessions_shared_cluster8_Mcap, fetch_protein_info)
final_shared_cluster8_Mcap <- do.call(rbind, results_shared_cluster8_Mcap)

results_shared_cluster9_Mcap  <- lapply(accessions_shared_cluster9_Mcap, fetch_protein_info)
final_shared_cluster9_Mcap <- do.call(rbind, results_shared_cluster9_Mcap)

results_unique_cluster5_30_Mcap  <- lapply(accessions_unique_cluster5_30_Mcap, fetch_protein_info)
final_unique_cluster5_30_Mcap <- do.call(rbind, results_unique_cluster5_30_Mcap)

results_unique_cluster5_35_Mcap  <- lapply(accessions_unique_cluster5_35_Mcap, fetch_protein_info)  
final_unique_cluster5_35_Mcap <- do.call(rbind, results_unique_cluster5_35_Mcap)

results_unique_cluster5_control_Mcap  <- lapply(accessions_unique_cluster5_control_Mcap, fetch_protein_info)
final_unique_cluster5_control_Mcap <- do.call(rbind, results_unique_cluster5_control_Mcap)

results_unique_cluster8_30_Mcap  <- lapply(accessions_unique_cluster8_30_Mcap, fetch_protein_info)
final_unique_cluster8_30_Mcap <- do.call(rbind, results_unique_cluster8_30_Mcap)

results_unique_cluster8_35_Mcap  <- lapply(accessions_unique_cluster8_35_Mcap, fetch_protein_info)
final_unique_cluster8_35_Mcap <- do.call(rbind, results_unique_cluster8_35_Mcap)

results_unique_cluster8_control_Mcap  <- lapply(accessions_unique_cluster8_control_Mcap, fetch_protein_info)
final_unique_cluster8_control_Mcap <- do.call(rbind, results_unique_cluster8_control_Mcap)

results_unique_cluster9_30_Mcap  <- lapply(accessions_unique_cluster9_30_Mcap, fetch_protein_info)
final_unique_cluster9_30_Mcap <- do.call(rbind, results_unique_cluster9_30_Mcap)

results_unique_cluster9_35_Mcap  <- lapply(accessions_unique_cluster9_35_Mcap, fetch_protein_info)
final_unique_cluster9_35_Mcap <- do.call(rbind, results_unique_cluster9_35_Mcap)

results_unique_cluster9_control_Mcap  <- lapply(accessions_unique_cluster9_control_Mcap, fetch_protein_info)
final_unique_cluster9_control_Mcap <- do.call(rbind, results_unique_cluster9_control_Mcap)

# keep only first 2 columns
final_shared_cluster5_Mcap <- final_shared_cluster5_Mcap[, 1:2]
final_shared_cluster8_Mcap <- final_shared_cluster8_Mcap[, 1:2]
final_shared_cluster9_Mcap <- final_shared_cluster9_Mcap[, 1:2]

final_unique_cluster5_30_Mcap <- final_unique_cluster5_30_Mcap[, 1:2]
final_unique_cluster5_35_Mcap <- final_unique_cluster5_35_Mcap[, 1:2]
final_unique_cluster5_control_Mcap <- final_unique_cluster5_control_Mcap[, 1:2]

final_unique_cluster8_30_Mcap <- final_unique_cluster8_30_Mcap[, 1:2]
final_unique_cluster8_35_Mcap <- final_unique_cluster8_35_Mcap[, 1:2]
final_unique_cluster8_control_Mcap <- final_unique_cluster8_control_Mcap[, 1:2]

final_unique_cluster9_30_Mcap <- final_unique_cluster9_30_Mcap[, 1:2]
final_unique_cluster9_35_Mcap <- final_unique_cluster9_35_Mcap[, 1:2]
final_unique_cluster9_control_Mcap <- final_unique_cluster9_control_Mcap[, 1:2]

# Join with gene table
colnames(final_shared_cluster5_Mcap)[1] <- "top_hit"
colnames(final_shared_cluster8_Mcap)[1] <- "top_hit"
colnames(final_shared_cluster9_Mcap)[1] <- "top_hit"

colnames(final_unique_cluster5_30_Mcap)[1] <- "top_hit"
colnames(final_unique_cluster5_35_Mcap)[1] <- "top_hit"
colnames(final_unique_cluster5_control_Mcap)[1] <- "top_hit"

colnames(final_unique_cluster8_30_Mcap)[1] <- "top_hit"
colnames(final_unique_cluster8_35_Mcap)[1] <- "top_hit"
colnames(final_unique_cluster8_control_Mcap)[1] <- "top_hit"

colnames(final_unique_cluster9_30_Mcap)[1] <- "top_hit"
colnames(final_unique_cluster9_35_Mcap)[1] <- "top_hit"
colnames(final_unique_cluster9_control_Mcap)[1] <- "top_hit"

joined_cluster5_Mcap <- left_join(final_shared_cluster5_Mcap, shared_genes_blast_cluster5_Mcap, by = "top_hit")
joined_cluster8_Mcap <- left_join(final_shared_cluster8_Mcap, shared_genes_blast_cluster8_Mcap, by = "top_hit")
joined_cluster9_Mcap <- left_join(final_shared_cluster9_Mcap, shared_genes_blast_cluster9_Mcap, by = "top_hit")

joined_unique_cluster5_30_Mcap <- left_join(final_unique_cluster5_30_Mcap, unique_genes_blast_cluster5_30_Mcap, by = "top_hit")
joined_unique_cluster5_35_Mcap <- left_join(final_unique_cluster5_35_Mcap, unique_genes_blast_cluster5_35_Mcap, by = "top_hit")
joined_unique_cluster5_control_Mcap <- left_join(final_unique_cluster5_control_Mcap, unique_genes_blast_cluster5_control_Mcap, by = "top_hit")

joined_unique_cluster8_30_Mcap <- left_join(final_unique_cluster8_30_Mcap, unique_genes_blast_cluster8_30_Mcap, by = "top_hit")
joined_unique_cluster8_35_Mcap <- left_join(final_unique_cluster8_35_Mcap, unique_genes_blast_cluster8_35_Mcap, by = "top_hit")
joined_unique_cluster8_control_Mcap <- left_join(final_unique_cluster8_control_Mcap, unique_genes_blast_cluster8_control_Mcap, by = "top_hit")

joined_unique_cluster9_30_Mcap <- left_join(final_unique_cluster9_30_Mcap, unique_genes_blast_cluster9_30_Mcap, by = "top_hit")
joined_unique_cluster9_35_Mcap <- left_join(final_unique_cluster9_35_Mcap, unique_genes_blast_cluster9_35_Mcap, by = "top_hit")
joined_unique_cluster9_control_Mcap <- left_join(final_unique_cluster9_control_Mcap, unique_genes_blast_cluster9_control_Mcap, by = "top_hit")


#### keep unique rows from the joined tables

joined_unique_cluster5_30_Mcap <- joined_unique_cluster5_30_Mcap[!duplicated(joined_unique_cluster5_30_Mcap$top_hit), ]
joined_unique_cluster5_35_Mcap <- joined_unique_cluster5_35_Mcap[!duplicated(joined_unique_cluster5_35_Mcap$top_hit), ]
joined_unique_cluster5_control_Mcap <- joined_unique_cluster5_control_Mcap[!duplicated(joined_unique_cluster5_control_Mcap$top_hit), ]

joined_unique_cluster8_30_Mcap <- joined_unique_cluster8_30_Mcap[!duplicated(joined_unique_cluster8_30_Mcap$top_hit), ]
joined_unique_cluster8_35_Mcap <- joined_unique_cluster8_35_Mcap[!duplicated(joined_unique_cluster8_35_Mcap$top_hit), ]
joined_unique_cluster8_control_Mcap <- joined_unique_cluster8_control_Mcap[!duplicated(joined_unique_cluster8_control_Mcap$top_hit), ]

joined_unique_cluster9_30_Mcap <- joined_unique_cluster9_30_Mcap[!duplicated(joined_unique_cluster9_30_Mcap$top_hit), ]
joined_unique_cluster9_35_Mcap <- joined_unique_cluster9_35_Mcap[  !duplicated(joined_unique_cluster9_35_Mcap$top_hit), ]
joined_unique_cluster9_control_Mcap <- joined_unique_cluster9_control_Mcap[!duplicated(joined_unique_cluster9_control_Mcap$top_hit), ] 

#### Join to ortholog map to get gene names for Mcap

joined_cluster5_Mcap_ortho <- left_join(joined_cluster5_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_cluster8_Mcap_ortho <- left_join(joined_cluster8_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_cluster9_Mcap_ortho <- left_join(joined_cluster9_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))

joined_unique_cluster5_30_Mcap_ortho <- left_join(joined_unique_cluster5_30_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_unique_cluster5_35_Mcap_ortho <- left_join(joined_unique_cluster5_35_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_unique_cluster5_control_Mcap_ortho <- left_join(joined_unique_cluster5_control_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))

joined_unique_cluster8_30_Mcap_ortho <- left_join(joined_unique_cluster8_30_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_unique_cluster8_35_Mcap_ortho <- left_join(joined_unique_cluster8_35_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_unique_cluster8_control_Mcap_ortho <- left_join(joined_unique_cluster8_control_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))

joined_unique_cluster9_30_Mcap_ortho <- left_join(joined_unique_cluster9_30_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_unique_cluster9_35_Mcap_ortho <- left_join(joined_unique_cluster9_35_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))
joined_unique_cluster9_control_Mcap_ortho <- left_join(joined_unique_cluster9_control_Mcap, ortho_map, by = c("Gene" = "mcap_gene"))

## remove rows with NAs in the ortho_ID column
joined_cluster5_Mcap_ortho <- joined_cluster5_Mcap_ortho[!is.na(joined_cluster5_Mcap_ortho$orthoID), ]
joined_cluster8_Mcap_ortho <- joined_cluster8_Mcap_ortho[!is.na(joined_cluster8_Mcap_ortho$orthoID), ]
joined_cluster9_Mcap_ortho <- joined_cluster9_Mcap_ortho[!is.na(joined_cluster9_Mcap_ortho$orthoID), ]

joined_unique_cluster5_30_Mcap_ortho <- joined_unique_cluster5_30_Mcap_ortho[!is.na(joined_unique_cluster5_30_Mcap_ortho$orthoID), ]
joined_unique_cluster5_35_Mcap_ortho <- joined_unique_cluster5_35_Mcap_ortho[!is.na(joined_unique_cluster5_35_Mcap_ortho$orthoID), ]
joined_unique_cluster5_control_Mcap_ortho <- joined_unique_cluster5_control_Mcap_ortho[!is.na(joined_unique_cluster5_control_Mcap_ortho$orthoID), ]

joined_unique_cluster8_30_Mcap_ortho <- joined_unique_cluster8_30_Mcap_ortho[!is.na(joined_unique_cluster8_30_Mcap_ortho$orthoID), ]
joined_unique_cluster8_35_Mcap_ortho <- joined_unique_cluster8_35_Mcap_ortho[!is.na(joined_unique_cluster8_35_Mcap_ortho$orthoID), ]
joined_unique_cluster8_control_Mcap_ortho <- joined_unique_cluster8_control_Mcap_ortho[!is.na(joined_unique_cluster8_control_Mcap_ortho$orthoID), ]

joined_unique_cluster9_30_Mcap_ortho <- joined_unique_cluster9_30_Mcap_ortho[!is.na(joined_unique_cluster9_30_Mcap_ortho$orthoID), ]
joined_unique_cluster9_35_Mcap_ortho <- joined_unique_cluster9_35_Mcap_ortho[!is.na(joined_unique_cluster9_35_Mcap_ortho$orthoID), ]
joined_unique_cluster9_control_Mcap_ortho <- joined_unique_cluster9_control_Mcap_ortho[!is.na(joined_unique_cluster9_control_Mcap_ortho$orthoID), ]

## get percentage of genes lost in the ortholog mapping
percentage_shared_cluster5_Mcap <- nrow(joined_cluster5_Mcap_ortho) / nrow(joined_unique_cluster5_30_Mcap) * 100 #0.91
percentage_shared_cluster8_Mcap <- nrow(joined_cluster8_Mcap_ortho) / nrow(joined_unique_cluster8_30_Mcap) * 100 #0.42
percentage_shared_cluster9_Mcap <- nrow(joined_cluster9_Mcap_ortho) / nrow(joined_unique_cluster9_30_Mcap) * 100 #1.30

percentage_unique_cluster5_30_Mcap <- nrow(joined_unique_cluster5_30_Mcap_ortho) / nrow(joined_unique_cluster5_30_Mcap) * 100 #44.24
percentage_unique_cluster5_35_Mcap <- nrow(joined_unique_cluster5_35_Mcap_ortho) / nrow(joined_unique_cluster5_35_Mcap) * 100 #46.60
percentage_unique_cluster5_control_Mcap <- nrow(joined_unique_cluster5_control_Mcap_ortho) / nrow(joined_unique_cluster5_control_Mcap) * 100 #48.63

percentage_unique_cluster8_30_Mcap <- nrow(joined_unique_cluster8_30_Mcap_ortho) / nrow(joined_unique_cluster8_30_Mcap) * 100 #37.58
percentage_unique_cluster8_35_Mcap <- nrow(joined_unique_cluster8_35_Mcap_ortho) / nrow(joined_unique_cluster8_35_Mcap) * 100 #35.65
percentage_unique_cluster8_control_Mcap <- nrow(joined_unique_cluster8_control_Mcap_ortho) / nrow(joined_unique_cluster8_control_Mcap) * 100 #33.63

percentage_unique_cluster9_30_Mcap <- nrow(joined_unique_cluster9_30_Mcap_ortho) / nrow(joined_unique_cluster9_30_Mcap) * 100 #47.40
percentage_unique_cluster9_35_Mcap <- nrow(joined_unique_cluster9_35_Mcap_ortho) / nrow(joined_unique_cluster9_35_Mcap) * 100 #43.64
percentage_unique_cluster9_control_Mcap <- nrow(joined_unique_cluster9_control_Mcap_ortho) / nrow(joined_unique_cluster9_control_Mcap) * 100 #53.53

## set row names of DE table as orthoID
DE_mcap_vs_pacu_complete <- DE_mcap_vs_pacu
DE_mcap_vs_pacu_complete$orthoID <- rownames(DE_mcap_vs_pacu)


# join DE table witg joined tables to get gene names
DE_mcap_vs_pacu_sharedHubs_cluster5 <- left_join(DE_mcap_vs_pacu_complete, joined_cluster5_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_sharedHubs_cluster8 <- left_join(DE_mcap_vs_pacu_complete, joined_cluster8_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_sharedHubs_cluster9 <- left_join(DE_mcap_vs_pacu_complete, joined_cluster9_Mcap_ortho , by = "orthoID") %>% na.omit()

DE_mcap_vs_pacu_uniqueHubs_cluster5_30 <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster5_30_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster5_35 <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster5_35_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster5_control <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster5_control_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster8_30 <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster8_30_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster8_35 <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster8_35_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster8_control <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster8_control_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster9_30 <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster9_30_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster9_35 <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster9_35_Mcap_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pacu_uniqueHubs_cluster9_control <- left_join(DE_mcap_vs_pacu_complete, joined_unique_cluster9_control_Mcap_ortho , by = "orthoID") %>% na.omit()






##load gene clusters - Pacu

#blue
blue_hub_gene_counts_control_pacu <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/top10pct_hub_genes_control_blue_4.csv")
blue_hub_gene_counts_30_pacu <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/top10pct_hub_genes_treat30_blue_4.csv")
blue_hub_gene_counts_35_pacu <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/top10pct_hub_genes_treat35_blue_4.csv")

#turquoise
turquoise_hub_gene_counts_control_pacu <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/top10pct_hub_genes_control_turquoise_3.csv")
turquoise_hub_gene_counts_30_pacu <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/top10pct_hub_genes_treat30_turquoise_3.csv")
turquoise_hub_gene_counts_35_pacu <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/top10pct_hub_genes_treat35_turquoise_3.csv")

### keep only unique hub genes per each temperature and module

# create vectors of gene IDs
blue_hub_genes_control_pacu <- blue_hub_gene_counts_control_pacu$Gene
blue_hub_genes_30_pacu <- blue_hub_gene_counts_30_pacu$Gene
blue_hub_genes_35_pacu <- blue_hub_gene_counts_35_pacu$Gene

turquoise_hub_genes_control_pacu <- turquoise_hub_gene_counts_control_pacu$Gene
turquoise_hub_genes_30_pacu <- turquoise_hub_gene_counts_30_pacu$Gene
turquoise_hub_genes_35_pacu <- turquoise_hub_gene_counts_35_pacu$Gene






##load gene clusters - Pcom

#cluster1
cluster1_hub_gene_counts_control_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster1.csv")
cluster1_hub_gene_counts_30_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster1.csv")
cluster1_hub_gene_counts_35_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster1.csv")

#cluster2
cluster2_hub_gene_counts_control_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster2.csv")
cluster2_hub_gene_counts_30_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster2.csv")
cluster2_hub_gene_counts_35_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster2.csv")

#cluster3
cluster3_hub_gene_counts_control_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster3.csv")
cluster3_hub_gene_counts_30_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster3.csv")
cluster3_hub_gene_counts_35_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster3.csv")

#cluster9
cluster9_hub_gene_counts_control_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_control_cluster9.csv")
cluster9_hub_gene_counts_30_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat30_cluster9.csv")
cluster9_hub_gene_counts_35_Pcom <- read_csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/top10pct_hub_genes_treat35_cluster9.csv")

# Get shared hub genes across the 3 temperatures
cluster1_hub_genes_shared_Pcom <- Reduce(intersect, list(
  cluster1_hub_gene_counts_control_Pcom,
  cluster1_hub_gene_counts_30_Pcom,
  cluster1_hub_gene_counts_35_Pcom
))

cluster2_hub_genes_shared_Pcom <- Reduce(intersect, list(
  cluster2_hub_gene_counts_control_Pcom,
  cluster2_hub_gene_counts_30_Pcom,
  cluster2_hub_gene_counts_35_Pcom
))

cluster3_hub_genes_shared_Pcom <- Reduce(intersect, list(
  cluster3_hub_gene_counts_control_Pcom,
  cluster3_hub_gene_counts_30_Pcom,
  cluster3_hub_gene_counts_35_Pcom
))

cluster9_hub_genes_shared_Pcom <- Reduce(intersect, list(
  cluster9_hub_gene_counts_control_Pcom,
  cluster9_hub_gene_counts_30_Pcom,
  cluster9_hub_gene_counts_35_Pcom
))

# Get all genes at 35 minus the shared genes across the 3 temperatures  
cluster1_hub_genes_35_unique_Pcom <- setdiff(
  cluster1_hub_gene_counts_35_Pcom,
  cluster1_hub_genes_shared_Pcom
)

cluster2_hub_genes_35_unique_Pcom <- setdiff(
  cluster2_hub_gene_counts_35_Pcom,
  cluster2_hub_genes_shared_Pcom
)

cluster3_hub_genes_35_unique_Pcom <- setdiff(
  cluster3_hub_gene_counts_35_Pcom,
  cluster3_hub_genes_shared_Pcom
)

cluster9_hub_genes_35_unique_Pcom <- setdiff(
  cluster9_hub_gene_counts_35_Pcom,
  cluster9_hub_genes_shared_Pcom
)

#Load diamond results
blast <- read_tsv("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pcom/Pcom.diamondBlastpNCBInr", col_names = FALSE)
colnames(blast) <- c("Gene", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(blast)
dim(blast)
#[1] 40740    12

### find top_hit from blast output corresponding to the hub genes
# Create a data frame from the gene list
shared_genes_cluster1_Pcom <- data.frame(Gene = cluster1_hub_genes_shared_Pcom, stringsAsFactors = FALSE)
shared_genes_cluster2_Pcom <- data.frame(Gene = cluster2_hub_genes_shared_Pcom, stringsAsFactors = FALSE)
shared_genes_cluster3_Pcom <- data.frame(Gene = cluster3_hub_genes_shared_Pcom, stringsAsFactors = FALSE)
shared_genes_cluster9_Pcom <- data.frame(Gene = cluster9_hub_genes_shared_Pcom, stringsAsFactors = FALSE)

unique_genes_35_cluster1_Pcom <- data.frame(Gene = cluster1_hub_genes_35_unique_Pcom, stringsAsFactors = FALSE)
unique_genes_35_cluster2_Pcom <- data.frame(Gene = cluster2_hub_genes_35_unique_Pcom, stringsAsFactors = FALSE)
unique_genes_35_cluster3_Pcom <- data.frame(Gene = cluster3_hub_genes_35_unique_Pcom, stringsAsFactors = FALSE)
unique_genes_35_cluster9_Pcom <- data.frame(Gene = cluster9_hub_genes_35_unique_Pcom, stringsAsFactors = FALSE)

# Inner join to get matching rows
shared_genes_blast_cluster1_Pcom  <- inner_join(shared_genes_cluster1_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster2_Pcom  <- inner_join(shared_genes_cluster2_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster3_Pcom  <- inner_join(shared_genes_cluster3_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
shared_genes_blast_cluster9_Pcom  <- inner_join(shared_genes_cluster9_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

unique_genes_blast_35_cluster1_Pcom  <- inner_join(unique_genes_35_cluster1_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
unique_genes_blast_35_cluster2_Pcom  <- inner_join(unique_genes_35_cluster2_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
unique_genes_blast_35_cluster3_Pcom  <- inner_join(unique_genes_35_cluster3_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)
unique_genes_blast_35_cluster9_Pcom  <- inner_join(unique_genes_35_cluster9_Pcom, blast, by = "Gene") %>%
  dplyr::select(Gene, top_hit)

#search proteins in NCBI
library(rentrez)

# Get the vector of accession numbers
accessions_shared_cluster1_Pcom <- shared_genes_blast_cluster1_Pcom$top_hit
accessions_shared_cluster2_Pcom <- shared_genes_blast_cluster2_Pcom$top_hit
accessions_shared_cluster3_Pcom <- shared_genes_blast_cluster3_Pcom$top_hit
accessions_shared_cluster9_Pcom <- shared_genes_blast_cluster9_Pcom$top_hit

accessions_unique_cluster1_35_Pcom <- unique_genes_blast_35_cluster1_Pcom$top_hit
accessions_unique_cluster2_35_Pcom <- unique_genes_blast_35_cluster2_Pcom$top_hit
accessions_unique_cluster3_35_Pcom <- unique_genes_blast_35_cluster3_Pcom$top_hit
accessions_unique_cluster9_35_Pcom <- unique_genes_blast_35_cluster9_Pcom$top_hit

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
results_shared_cluster1_Pcom  <- lapply(accessions_shared_cluster1_Pcom, fetch_protein_info)
final_shared_cluster1_Pcom <- do.call(rbind, results_shared_cluster1_Pcom)

results_shared_cluster2_Pcom  <- lapply(accessions_shared_cluster2_Pcom, fetch_protein_info)
final_shared_cluster2_Pcom <- do.call(rbind, results_shared_cluster2_Pcom)

results_shared_cluster3_Pcom  <- lapply(accessions_shared_cluster3_Pcom, fetch_protein_info)
final_shared_cluster3_Pcom <- do.call(rbind, results_shared_cluster3_Pcom)

results_shared_cluster9_Pcom  <- lapply(accessions_shared_cluster9_Pcom, fetch_protein_info)
final_shared_cluster9_Pcom <- do.call(rbind, results_shared_cluster9_Pcom)

results_unique_cluster1_35_Pcom  <- lapply(accessions_unique_cluster1_35_Pcom, fetch_protein_info)
final_unique_cluster1_35_Pcom <- do.call(rbind, results_unique_cluster1_35_Pcom)

results_unique_cluster2_35_Pcom  <- lapply(accessions_unique_cluster2_35_Pcom, fetch_protein_info)
final_unique_cluster2_35_Pcom <- do.call(rbind, results_unique_cluster2_35_Pcom)

results_unique_cluster3_35_Pcom  <- lapply(accessions_unique_cluster3_35_Pcom, fetch_protein_info)
final_unique_cluster3_35_Pcom <- do.call(rbind, results_unique_cluster3_35_Pcom)

results_unique_cluster9_35_Pcom  <- lapply(accessions_unique_cluster9_35_Pcom, fetch_protein_info)
final_unique_cluster9_35_Pcom <- do.call(rbind, results_unique_cluster9_35_Pcom)

# keep only first 2 columns
final_shared_cluster1_Pcom <- final_shared_cluster1_Pcom[, 1:2]
final_shared_cluster2_Pcom <- final_shared_cluster2_Pcom[, 1:2]
final_shared_cluster3_Pcom <- final_shared_cluster3_Pcom[, 1:2]
final_shared_cluster9_Pcom <- final_shared_cluster9_Pcom[, 1:2]

final_unique_cluster1_35_Pcom <- final_unique_cluster1_35_Pcom[, 1:2]
final_unique_cluster2_35_Pcom <- final_unique_cluster2_35_Pcom[, 1:2]
final_unique_cluster3_35_Pcom <- final_unique_cluster3_35_Pcom[, 1:2]
final_unique_cluster9_35_Pcom <- final_unique_cluster9_35_Pcom[, 1:2]

# Join with gene table
colnames(final_shared_cluster1_Pcom)[1] <- "top_hit"
colnames(final_shared_cluster2_Pcom)[1] <- "top_hit"
colnames(final_shared_cluster3_Pcom)[1] <- "top_hit"
colnames(final_shared_cluster9_Pcom)[1] <- "top_hit"

colnames(final_unique_cluster1_35_Pcom)[1] <- "top_hit"
colnames(final_unique_cluster2_35_Pcom)[1] <- "top_hit"
colnames(final_unique_cluster3_35_Pcom)[1] <- "top_hit"
colnames(final_unique_cluster9_35_Pcom)[1] <- "top_hit"

joined_shared_cluster1_Pcom <- left_join(final_shared_cluster1_Pcom, shared_genes_blast_cluster1_Pcom, by = "top_hit")
joined_shared_cluster2_Pcom <- left_join(final_shared_cluster2_Pcom, shared_genes_blast_cluster2_Pcom, by = "top_hit")
joined_shared_cluster3_Pcom <- left_join(final_shared_cluster3_Pcom, shared_genes_blast_cluster3_Pcom, by = "top_hit")
joined_shared_cluster9_Pcom <- left_join(final_shared_cluster9_Pcom, shared_genes_blast_cluster9_Pcom, by = "top_hit")

joined_unique_cluster1_35_Pcom <- left_join(final_unique_cluster1_35_Pcom, unique_genes_blast_35_cluster1_Pcom, by = "top_hit")
joined_unique_cluster2_35_Pcom <- left_join(final_unique_cluster2_35_Pcom, unique_genes_blast_35_cluster2_Pcom, by = "top_hit")
joined_unique_cluster3_35_Pcom <- left_join(final_unique_cluster3_35_Pcom, unique_genes_blast_35_cluster3_Pcom, by = "top_hit")
joined_unique_cluster9_35_Pcom <- left_join(final_unique_cluster9_35_Pcom, unique_genes_blast_35_cluster9_Pcom, by = "top_hit")

#### keep unique rows from the joined tables
joined_unique_cluster1_35_Pcom <- joined_unique_cluster1_35_Pcom[!duplicated(joined_unique_cluster1_35_Pcom$top_hit), ]
joined_unique_cluster2_35_Pcom <- joined_unique_cluster2_35_Pcom[!duplicated(joined_unique_cluster2_35_Pcom$top_hit), ]
joined_unique_cluster3_35_Pcom <- joined_unique_cluster3_35_Pcom[!duplicated(joined_unique_cluster3_35_Pcom$top_hit), ]
joined_unique_cluster9_35_Pcom <- joined_unique_cluster9_35_Pcom[!duplicated(joined_unique_cluster9_35_Pcom$top_hit), ]

#### Join to ortholog map to get gene names for Pcom
joined_shared_cluster1_Pcom_ortho <- left_join(joined_shared_cluster1_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))
joined_shared_cluster2_Pcom_ortho <- left_join(joined_shared_cluster2_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))
joined_shared_cluster3_Pcom_ortho <- left_join(joined_shared_cluster3_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))
joined_shared_cluster9_Pcom_ortho <- left_join(joined_shared_cluster9_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))

joined_unique_cluster1_35_Pcom_ortho <- left_join(joined_unique_cluster1_35_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))
joined_unique_cluster2_35_Pcom_ortho <- left_join(joined_unique_cluster2_35_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))  
joined_unique_cluster3_35_Pcom_ortho <- left_join(joined_unique_cluster3_35_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))
joined_unique_cluster9_35_Pcom_ortho <- left_join(joined_unique_cluster9_35_Pcom, ortho_map, by = c("Gene" = "pcom_gene"))

## remove rows with NAs in the ortho_ID column
joined_shared_cluster1_Pcom_ortho <- joined_shared_cluster1_Pcom_ortho[!is.na(joined_shared_cluster1_Pcom_ortho$orthoID), ]
joined_shared_cluster2_Pcom_ortho <- joined_shared_cluster2_Pcom_ortho[!is.na(joined_shared_cluster2_Pcom_ortho$orthoID), ]
joined_shared_cluster3_Pcom_ortho <- joined_shared_cluster3_Pcom_ortho[!is.na(joined_shared_cluster3_Pcom_ortho$orthoID), ]
joined_shared_cluster9_Pcom_ortho <- joined_shared_cluster9_Pcom_ortho[!is.na(joined_shared_cluster9_Pcom_ortho$orthoID), ]

joined_unique_cluster1_35_Pcom_ortho <- joined_unique_cluster1_35_Pcom_ortho[!is.na(joined_unique_cluster1_35_Pcom_ortho$orthoID), ]
joined_unique_cluster2_35_Pcom_ortho <- joined_unique_cluster2_35_Pcom_ortho[!is.na(joined_unique_cluster2_35_Pcom_ortho$orthoID), ]
joined_unique_cluster3_35_Pcom_ortho <- joined_unique_cluster3_35_Pcom_ortho[!is.na(joined_unique_cluster3_35_Pcom_ortho$orthoID), ]
joined_unique_cluster9_35_Pcom_ortho <- joined_unique_cluster9_35_Pcom_ortho[!is.na(joined_unique_cluster9_35_Pcom_ortho$orthoID), ]

## get percentage of genes lost in the ortholog mapping
percentage_shared_cluster1_Pcom <- nrow(joined_shared_cluster1_Pcom_ortho) / nrow(joined_shared_cluster1_Pcom) * 100 #66.67
percentage_shared_cluster2_Pcom <- nrow(joined_shared_cluster2_Pcom_ortho) / nrow(joined_shared_cluster2_Pcom) * 100 #0
percentage_shared_cluster3_Pcom <- nrow(joined_shared_cluster3_Pcom_ortho) / nrow(joined_shared_cluster3_Pcom) * 100 #40
percentage_shared_cluster9_Pcom <- nrow(joined_shared_cluster9_Pcom_ortho) / nrow(joined_shared_cluster9_Pcom) * 100 #100

percentage_unique_cluster1_Pcom <- nrow(joined_unique_cluster1_35_Pcom_ortho) / nrow(joined_unique_cluster1_35_Pcom) * 100 #49.35
percentage_unique_cluster2_Pcom <- nrow(joined_unique_cluster2_35_Pcom_ortho) / nrow(joined_unique_cluster2_35_Pcom) * 100 #7.17
percentage_unique_cluster3_Pcom <- nrow(joined_unique_cluster3_35_Pcom_ortho) / nrow(joined_unique_cluster3_35_Pcom) * 100 #37.76
percentage_unique_cluster9_Pcom <- nrow(joined_unique_cluster9_35_Pcom_ortho) / nrow(joined_unique_cluster9_35_Pcom) * 100 #28.85

## set row names of DE table as orthoID
DE_pacu_vs_pcom_complete <- DE_pacu_vs_pcom
DE_pacu_vs_pcom_complete$orthoID <- rownames(DE_pacu_vs_pcom)

DE_mcap_vs_pcom_complete <- DE_mcap_vs_pcom
DE_mcap_vs_pcom_complete$orthoID <- rownames(DE_mcap_vs_pcom)

# join DE table witg joined tables to get gene names
DE_pacu_vs_pcom_sharedHubs_cluster1 <- left_join(DE_pacu_vs_pcom_complete, joined_shared_cluster1_Pcom_ortho , by = "orthoID") %>% na.omit()
## 4 all unnamed proteins [Porites evermanni]
DE_pacu_vs_pcom_sharedHubs_cluster2 <- left_join(DE_pacu_vs_pcom_complete, joined_shared_cluster2_Pcom_ortho , by = "orthoID") %>% na.omit()
## 0
DE_pacu_vs_pcom_sharedHubs_cluster3 <- left_join(DE_pacu_vs_pcom_complete, joined_shared_cluster3_Pcom_ortho , by = "orthoID") %>% na.omit()
## 2 all unnamed
DE_pacu_vs_pcom_sharedHubs_cluster9 <- left_join(DE_pacu_vs_pcom_complete, joined_shared_cluster9_Pcom_ortho , by = "orthoID") %>% na.omit()
## 1 unnamed

DE_pacu_vs_pcom_uniqueHubs_cluster1_35 <- left_join(DE_pacu_vs_pcom_complete, joined_unique_cluster1_35_Pcom_ortho , by = "orthoID") %>% na.omit()
## 406 unnamed, 10 with names
DE_pacu_vs_pcom_uniqueHubs_cluster2_35 <- left_join(DE_pacu_vs_pcom_complete, joined_unique_cluster2_35_Pcom_ortho , by = "orthoID") %>% na.omit()
## 22 all unnamed
DE_pacu_vs_pcom_uniqueHubs_cluster3_35 <- left_join(DE_pacu_vs_pcom_complete, joined_unique_cluster3_35_Pcom_ortho , by = "orthoID") %>% na.omit()
## 232 unnamed, 7 with names
DE_pacu_vs_pcom_uniqueHubs_cluster9_35 <- left_join(DE_pacu_vs_pcom_complete, joined_unique_cluster9_35_Pcom_ortho , by = "orthoID") %>% na.omit()
## 13 unnamed, 2 with names

DE_mcap_vs_pcom_sharedHubs_cluster1 <- left_join(DE_mcap_vs_pcom_complete, joined_shared_cluster1_Pcom_ortho , by = "orthoID") %>% na.omit()
## 4 all unnamed proteins 
DE_mcap_vs_pcom_sharedHubs_cluster2 <- left_join(DE_mcap_vs_pcom_complete, joined_shared_cluster2_Pcom_ortho , by = "orthoID") %>% na.omit()
## 0
DE_mcap_vs_pcom_sharedHubs_cluster3 <- left_join(DE_mcap_vs_pcom_complete, joined_shared_cluster3_Pcom_ortho , by = "orthoID") %>% na.omit()
## 2 all unnamed
DE_mcap_vs_pcom_sharedHubs_cluster9 <- left_join(DE_mcap_vs_pcom_complete, joined_shared_cluster9_Pcom_ortho , by = "orthoID") %>% na.omit()
## 1 unnamed

DE_mcap_vs_pcom_uniqueHubs_cluster1_35 <- left_join(DE_mcap_vs_pcom_complete, joined_unique_cluster1_35_Pcom_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pcom_uniqueHubs_cluster2_35 <- left_join(DE_mcap_vs_pcom_complete, joined_unique_cluster2_35_Pcom_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pcom_uniqueHubs_cluster3_35 <- left_join(DE_mcap_vs_pcom_complete, joined_unique_cluster3_35_Pcom_ortho , by = "orthoID") %>% na.omit()
DE_mcap_vs_pcom_uniqueHubs_cluster9_35 <- left_join(DE_mcap_vs_pcom_complete, joined_unique_cluster9_35_Pcom_ortho , by = "orthoID") %>% na.omit()


## many unnamed proteins form Pcom blast, I will look at the blast hit for Mcap and Pacu to see if they are the same for the orthologs
## the DE tables above also have the gene names for Mcap and Pacu, so I will add the blast hit for those as well

blast_mcap <- read_tsv("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Mcap/Mcap.diamondBlastpNCBInr", col_names = FALSE)
colnames(blast_mcap) <- c("Gene", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blast_pacu <- read_tsv("/scratch/workspace/federica_scucchia_uri_edu-hawaii/20250424_ENCORE_HawaiiTPC_Federica/output/Functional_Annotation/Pacu/Pacu.diamondBlastpNCBInr", col_names = FALSE)
colnames(blast_pacu) <- c("Gene", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

## add protein names to the blast hits
#search proteins in NCBI
library(rentrez)
# Get the vector of accession numbers
accessions_pacu <- blast_pacu$top_hit
accessions_mcap <- blast_mcap$top_hit

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
results_accessions_pacu  <- lapply(accessions_pacu, fetch_protein_info)
final_accessions_pacu <- do.call(rbind, results_accessions_pacu)

results_accessions_mcap  <- lapply(accessions_mcap, fetch_protein_info)
final_accessions_mcap <- do.call(rbind, results_accessions_mcap)
