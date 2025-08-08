
############# "Host ash-free dry weight measurements" #############

## load libraries
library("ggplot2")
library("plotrix")
library('dplyr')
library(tidyverse)
library(cowplot)
library(broom)
library(here)


## Import metadata
Data <- read_delim(here("Samples_ASFDW_raw.csv"), delim = ";")

## Import surface area data
Surf_Data <- read_delim(here("1_surface_area_final.csv"), delim = ",")

## Import homogenate volumes data
HomogVol_Data <- read_delim(here("Homog_volumes.csv"), delim = ";")

## Join datasets based on Sample_ID
Combined_Data <- Data %>%
  left_join(Surf_Data, by = "colony_id") %>%
  left_join(HomogVol_Data, by = "colony_id") %>%
  distinct(colony_id, .keep_all = TRUE)

## Add a new column "Tot DW g" to Combined_Data
Combined_Data <- Combined_Data %>%
  mutate(`Tot DW g` = `AFDW (g/ml)` * `PBS ml`)

## Add a new column "Dry weight per area" to Combined_Data
Combined_Data <- Combined_Data %>%
  mutate(`Dry weight per area_g/cm2` = `Tot DW g` / `surface.area.cm2`)

Combined_Data <- Combined_Data %>%
  mutate(`Dry weight per area_ug/cm2` = `Dry weight per area_g/cm2` * 1000000)

## Add a new column with species names
Combined_Data <- Combined_Data %>%
  mutate(species = case_when(
    grepl("MCAP", colony_id) ~ "Montipora capitata",
    grepl("PACT", colony_id) ~ "Pocillopora acuta",
    grepl("PCOM", colony_id) ~ "Porites compressa",
    TRUE ~ NA_character_  # Default to NA if no match
  ))


# Plot results
Combined_Data %>% 
  ggplot(aes(x = species, y = `Dry weight per area_ug/cm2`, color = species, group =species)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +  # Add jittered points
  #geom_text(aes(label = colony_id), hjust = 1.5, vjust = 1.5, size = 3) +  # Add colony_id labels
  theme_minimal() +  # Set a minimal theme with a white background
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_text(face = "bold"),  # Customize facet labels
    legend.position = "right",     # Add legend for colors
    panel.border = element_rect(color = "black", fill = NA),  # Add black borders
    axis.ticks = element_line(color = "black")  # Show axis ticks
  ) +
  scale_color_manual(values = c(
    "Montipora capitata" = "green",
    "Pocillopora acuta" = "cyan",
    "Porites compressa" = "orange"
  ))

############ Remove outliers 

#set quantile values 
q <- c(0.25, 0.75) 

Quants_AFDW <- Combined_Data %>%
  group_by(species) %>%
  summarise(quant25 = quantile(`Dry weight per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`Dry weight per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`Dry weight per area_ug/cm2`))

#Calculate Quantile upper and lower ranges
Quants_AFDW$upper <-  Quants_AFDW$quant75+1.5*Quants_AFDW$IQRbyGroup # Upper Range
Quants_AFDW$lower <- Quants_AFDW$quant25-1.5*Quants_AFDW$IQRbyGroup # Lower Range

# Join quantile ranges back to the original dataset
Combined_Data_with_bounds <- Combined_Data %>%
  left_join(Quants_AFDW, by = "species")

# Filter out outliers based on the quantile ranges
x1 <- Combined_Data_with_bounds %>%
  filter(`Dry weight per area_ug/cm2` < upper & `Dry weight per area_ug/cm2` > lower)

# Plot results
x1  %>% 
  ggplot(aes(x = species, y = `Dry weight per area_ug/cm2`, color = species, group =species)) + 
  geom_boxplot(outlier.shape = NA) +
  #geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +  # Add jittered points
  #geom_text(aes(label = colony_id), hjust = 1.5, vjust = 1.5, size = 3) +  # Add colony_id labels
  theme_minimal() +  # Set a minimal theme with a white background
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_text(face = "bold"),  # Customize facet labels
    legend.position = "right",     # Add legend for colors
    panel.border = element_rect(color = "black", fill = NA),  # Add black borders
    axis.ticks = element_line(color = "black")  # Show axis ticks
  ) +
  scale_color_manual(values = c(
    "Montipora capitata" = "green",
    "Pocillopora acuta" = "cyan",
    "Porites compressa" = "orange"
  ))


##### change graph type

x1  %>% 
  ggplot(aes(x = species, y = `Dry weight per area_ug/cm2`, color = species, fill = species)) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  #theme(axis.text.x = element_text(vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata", 
                              "Pocillopora acuta", "Porites compressa"))  +
  labs(y = "Host Dry weight per area (ug/cm^2)") +  # Add y-axis label
  geom_boxplot(
    aes(fill = species, fill = after_scale(colorspace::lighten(fill, .7))),
    size = 0.5  ) +
  geom_point(
    position = position_jitter(width = .2, seed = 0),
    size = 2.5, alpha = .5
  ) +
  geom_point(
    position = position_jitter(width = .2, seed = 0),
    size = 2.5, stroke = .5, shape = 1, color = "white"
  )


library(writexl)

# Save x1 as an Excel file
write_xlsx(x1, path = "AFDW_Rdata.xlsx")



################################
  
  ###### **Statistical Analysis**     
  
  ## Load packages
  library(dplyr)
  library(ggpubr)
  library(onewaytests)
  library(purrr)
  library(tidyr)
  
  ## Assess the normality of the data per each variable
  #Shapiro-Wilk’s method is widely recommended for normality test and it provides better power than K-S. 
  #It is based on the correlation between the data and the corresponding normal scores.
  
  # Ensure Temp.Cat is a factor
  x1$species <- as.factor(x1$species)

  # Perform Shapiro-Wilk test for normality for each species
  shapiro_results <- x1 %>%
    group_by(species) %>%
    summarise(
      p_value = shapiro.test(`Dry weight per area_ug/cm2`)$p.value
    )
  
  # View the results
  print(shapiro_results)
  # # A tibble: 3 × 2
  # species            p_value
  # <fct>                <dbl>
  # 1 Montipora capitata  0.542 
  # 2 Pocillopora acuta   0.703 
  # 3 Porites compressa   0.0410
  
  
  # Generate Q-Q plots for each species
  qqplots <- x1 %>%
    group_by(species) %>%
    nest() %>%
    mutate(
      qqplot = map(data, ~ ggqqplot(.x$`Dry weight per area_ug/cm2`, 
                                    title = paste("Q-Q Plot for", unique(.x$species))))
    )
  
  # Save the Q-Q plots to files
  output_dir <- "output/QQplots_species"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (i in seq_len(nrow(qqplots))) {
    plot <- qqplots$qqplot[[i]]
    species <- qqplots$species[i]
    filename <- paste0(output_dir, "/qqplot_", species, ".png")
    ggsave(filename = filename, plot = plot, width = 8, height = 6)
  }
  
  ## perform the Levene test for homogeneity of variances
  
  library(car)
  
  levene_results <- leveneTest(`Dry weight per area_ug/cm2` ~ species, data = x1)
  
  # View the results
  print(levene_results)
  #   Levene's Test for Homogeneity of Variance (center = median)
  # Df F value Pr(>F)
  # group  2  2.3797 0.1211
  # 18  
  
  
  
  
  #### Perform Kruskal-Wallis test (not normally distributed data)
  
  kruskal_results <- kruskal.test(`Dry weight per area_ug/cm2` ~ species, data = x1)
  
  # View the results of the Kruskal-Wallis test
  
  print(kruskal_results)
  # Kruskal-Wallis rank sum test
  # 
  # data:  Dry weight per area_ug/cm2 by species
  # Kruskal-Wallis chi-squared = 14.083, df = 2, p-value = 0.0008748
  
  
  ## If significant, perform Dunn's test
  library(dunn.test)

  dunn_results <- dunn.test(
    x = x1$`Dry weight per area_ug/cm2`,
    g = x1$species,
    method = "bonferroni"
  )
  
  # Extract significant comparisons
  significant_results <- data.frame(
    comparisons = dunn_results$comparisons,
    p_value = dunn_results$P.adjusted
  ) %>%
    filter(p_value < 0.05) %>%
    separate(comparisons, into = c("group1", "group2"), sep = " - ") %>%
    mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*"
      )
    )
  
  # Calculate y-positions for the significance bars
  significant_results <- significant_results %>%
    left_join(
      x1 %>% group_by(species) %>% summarise(max_y = max(`Dry weight per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),
      by = c("group1" = "species")
    ) %>%
    rename(max_y_group1 = max_y) %>%
    left_join(
      x1 %>% group_by(species) %>% summarise(max_y = max(`Dry weight per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),
      by = c("group2" = "species")
    ) %>%
    rename(max_y_group2 = max_y) %>%
    mutate(
      base_y = pmax(max_y_group1, max_y_group2) + 200,  # Adjust distance from data points
      y.position = base_y + 100  # Add extra padding for asterisks
    )
  
  # Ensure group1 and group2 are characters
  significant_results <- significant_results %>%
    mutate(
      group1 = as.character(group1),
      group2 = as.character(group2)
    )
  
  # Adjust y-positions to avoid overlapping bars
  significant_results <- significant_results %>%
    arrange(group1, group2) %>%  # Ensure consistent ordering
    mutate(
      y.position = base_y + row_number() * 700 # Increment y.position for each comparison
    )
  
  # Create the plot with significance bars
  library(ggsignif)
  library(ggpubr)
  
  # Create the plot with statistical annotations
  x1 %>%
    ggplot(aes(x = species, y = `Dry weight per area_ug/cm2`, color = species, fill = species)) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(y = "Host Dry weight per area (ug/cm^2)") +
    geom_boxplot(
      aes(fill = species, fill = after_scale(colorspace::lighten(fill, .7))),
      size = 0.5
    ) +
    geom_point(
      position = position_jitter(width = .2, seed = 0),
      size = 2.5, alpha = .5
    ) +
    geom_point(
      position = position_jitter(width = .2, seed = 0),
      size = 2.5, stroke = .5, shape = 1, color = "white"
    ) +
    stat_pvalue_manual(
      significant_results,
      label = "label",  # Use the label column for significance levels
      xmin = "group1",  # Start of the comparison
      xmax = "group2",  # End of the comparison
      y.position = "y.position",  # Position of the significance bars
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )

  
 
  
#### save ANOVA and Tukey test results in Excel file

library(writexl)
library(broom)
library(dplyr)
library(FSA)  

# Create publication-ready Kruskal-Wallis table
kruskal_publication_table <- data.frame(
  Test = "Kruskal-Wallis",
  Statistic = round(kruskal_results$statistic, 3),
  df = kruskal_results$parameter,
  p_value = case_when(
    kruskal_results$p.value < 0.001 ~ "< 0.001",
    kruskal_results$p.value < 0.01 ~ as.character(round(kruskal_results$p.value, 3)),
    TRUE ~ as.character(round(kruskal_results$p.value, 3))
  ),
  significance = case_when(
    kruskal_results$p.value <= 0.001 ~ "***",
    kruskal_results$p.value <= 0.01 ~ "**",
    kruskal_results$p.value <= 0.05 ~ "*",
    TRUE ~ ""
  ),
  Method = kruskal_results$method
) %>%
  rename(
    "Chi-squared" = Statistic,
    "p-value" = p_value,
    "Significance" = significance
  )

# Create publication-ready Dunn's test table using FSA package
dunn_results_fsa <- dunnTest(`Dry weight per area_ug/cm2` ~ species, data = x1, method = "bonferroni")

dunn_publication_table <- dunn_results_fsa$res %>%
  mutate(
    # Round values appropriately
    Z = round(Z, 3),
    P.unadj = case_when(
      P.unadj < 0.001 ~ "< 0.001",
      P.unadj < 0.01 ~ as.character(round(P.unadj, 3)),
      TRUE ~ as.character(round(P.unadj, 3))
    ),
    P.adj = case_when(
      P.adj < 0.001 ~ "< 0.001",
      P.adj < 0.01 ~ as.character(round(P.adj, 3)),
      TRUE ~ as.character(round(P.adj, 3))
    ),
    # Add significance stars for adjusted p-values
    significance = case_when(
      as.numeric(gsub("< ", "", P.adj)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", P.adj)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", P.adj)) <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rename(
    "Comparison" = Comparison,
    "Z-statistic" = Z,
    "p-value (unadj)" = P.unadj,
    "p-value (adj)" = P.adj,
    "Significance" = significance
  )

# Create summary statistics table (using medians and IQR for non-parametric data)
summary_stats <- x1 %>%
  group_by(species) %>%
  summarise(
    N = n(),
    Median = round(median(`Dry weight per area_ug/cm2`, na.rm = TRUE), 3),
    IQR = round(IQR(`Dry weight per area_ug/cm2`, na.rm = TRUE), 3),
    Min = round(min(`Dry weight per area_ug/cm2`, na.rm = TRUE), 3),
    Max = round(max(`Dry weight per area_ug/cm2`, na.rm = TRUE), 3),
    Q1 = round(quantile(`Dry weight per area_ug/cm2`, 0.25, na.rm = TRUE), 3),
    Q3 = round(quantile(`Dry weight per area_ug/cm2`, 0.75, na.rm = TRUE), 3),
    Mean = round(mean(`Dry weight per area_ug/cm2`, na.rm = TRUE), 3),
    SD = round(sd(`Dry weight per area_ug/cm2`, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  rename("Species" = species)

# Create assumption testing results table
assumption_tests <- data.frame(
  Test = c("Shapiro-Wilk (Montipora capitata)", 
           "Shapiro-Wilk (Pocillopora acuta)", 
           "Shapiro-Wilk (Porites compressa)",
           "Levene's Test"),
  Purpose = c("Normality", "Normality", "Normality", "Homogeneity of Variance"),
  Statistic = c(
    round(shapiro.test(x1$`Dry weight per area_ug/cm2`[x1$species == "Montipora capitata"])$statistic, 4),
    round(shapiro.test(x1$`Dry weight per area_ug/cm2`[x1$species == "Pocillopora acuta"])$statistic, 4),
    round(shapiro.test(x1$`Dry weight per area_ug/cm2`[x1$species == "Porites compressa"])$statistic, 4),
    round(levene_results$`F value`[1], 4)
  ),
  p_value = c(
    round(shapiro_results$p_value[1], 4),
    round(shapiro_results$p_value[2], 4),
    round(shapiro_results$p_value[3], 4),
    round(levene_results$`Pr(>F)`[1], 4)
  ),
  Interpretation = c(
    ifelse(shapiro_results$p_value[1] > 0.05, "Normal", "Non-normal"),
    ifelse(shapiro_results$p_value[2] > 0.05, "Normal", "Non-normal"),
    ifelse(shapiro_results$p_value[3] > 0.05, "Normal", "Non-normal"),
    ifelse(levene_results$`Pr(>F)`[1] > 0.05, "Homogeneous", "Heterogeneous")
  ),
  Test_Used = c(
    ifelse(shapiro_results$p_value[1] > 0.05, "Parametric (ANOVA)", "Non-parametric (Kruskal-Wallis)"),
    ifelse(shapiro_results$p_value[2] > 0.05, "Parametric (ANOVA)", "Non-parametric (Kruskal-Wallis)"),
    ifelse(shapiro_results$p_value[3] > 0.05, "Parametric (ANOVA)", "Non-parametric (Kruskal-Wallis)"),
    "Overall: Non-parametric (Kruskal-Wallis)"
  )
) %>%
  rename("p-value" = p_value)

# Create comprehensive Excel workbook
afdw_stats_tables <- list(
  "Summary_Statistics" = summary_stats,
  "Kruskal_Wallis_Results" = kruskal_publication_table,
  "Dunn_PostHoc_Results" = dunn_publication_table,
  "Assumption_Tests" = assumption_tests,
  "Raw_Data" = x1 %>% 
    select(colony_id, species, `Dry weight per area_ug/cm2`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "DryWeight_per_area_ug_cm2" = `Dry weight per area_ug/cm2`
    )
)

# Save to Excel
write_xlsx(afdw_stats_tables, "output/AFDW_Statistical_Analysis.xlsx")

