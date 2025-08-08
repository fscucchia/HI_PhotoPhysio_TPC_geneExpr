
############# "Host carbohydrates concentration measurements" #############

## load libraries
library("ggplot2")
library("plotrix")
library('dplyr')
library(tidyverse)
library(cowplot)
library(broom)
library(here)


## Import metadata
Data <- read_delim(here("carbs_raw.csv"), delim = ";")

## Import surface area data
Surf_Data <- read_delim(here("1_surface_area_final.csv"), delim = ",")

## Import homogenate volumes data
HomogVol_Data <- read_delim(here("Homog_volumes.csv"), delim = ";")

## Join datasets based on Sample_ID
Combined_Data <- Data %>%
  left_join(Surf_Data, by = "colony_id") %>%
  left_join(HomogVol_Data, by = "colony_id") %>%
  distinct(colony_id, .keep_all = TRUE)

## Add a new column "Tot conc mg" to Combined_Data, multiply by dilution factor
Combined_Data <- Combined_Data %>%
  mutate(`Tot conc mg` = `carbs_mg/ml` * `PBS ml` * 10)

Combined_Data <- Combined_Data %>%
  mutate(`Tot conc mg` = `carbs_mg/ml` * `PBS ml`)

## Add a new column "Carbs per area" to Combined_Data
Combined_Data <- Combined_Data %>%
  mutate(`Carbs per area_ug/cm2` = `Tot conc mg` / `surface.area.cm2`* 1000)

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
  ggplot(aes(x = species, y = `Carbs per area_ug/cm2`, color = species, group =species)) + 
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

Quants_Carbs <- Combined_Data %>%
  group_by(species) %>%
  summarise(quant25 = quantile(`Carbs per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`Carbs per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`Carbs per area_ug/cm2`))

#Calculate Quantile upper and lower ranges
Quants_Carbs$upper <-  Quants_Carbs$quant75+1.5*Quants_Carbs$IQRbyGroup # Upper Range
Quants_Carbs$lower <- Quants_Carbs$quant25-1.5*Quants_Carbs$IQRbyGroup # Lower Range

# Join quantile ranges back to the original dataset
Combined_Data_with_bounds <- Combined_Data %>%
  left_join(Quants_Carbs, by = "species")

# Filter out outliers based on the quantile ranges
x1 <- Combined_Data_with_bounds %>%
  filter(`Carbs per area_ug/cm2` < upper & `Carbs per area_ug/cm2` > lower)

# Plot results
x1  %>% 
  ggplot(aes(x = species, y = `Carbs per area_ug/cm2`, color = species, group =species)) + 
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



############ Remove outliers 2nd round, Pacu


#set quantile values 
q <- c(0.25, 0.75) 

# Create a new table x2 with the first 8 columns of x1
x2 <- x1 %>% 
  select(1:8)

Quants_Carbs <- x2 %>%
  group_by(species) %>%
  summarise(quant25 = quantile(`Carbs per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`Carbs per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`Carbs per area_ug/cm2`))

#Calculate Quantile upper and lower ranges
Quants_Carbs$upper <-  Quants_Carbs$quant75+1.5*Quants_Carbs$IQRbyGroup # Upper Range
Quants_Carbs$lower <- Quants_Carbs$quant25-1.5*Quants_Carbs$IQRbyGroup # Lower Range

# Join quantile ranges back to the original dataset
Combined_Data_with_bounds <- x2 %>%
  left_join(Quants_Carbs, by = "species")

# Filter out outliers based on the quantile ranges
x2 <- Combined_Data_with_bounds %>%
  filter(`Carbs per area_ug/cm2` < upper & `Carbs per area_ug/cm2` > lower)

# Plot results
x2  %>% 
  ggplot(aes(x = species, y = `Carbs per area_ug/cm2`, color = species, group =species)) + 
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


##### change graph type

x2  %>% 
  ggplot(aes(x = species, y = `Carbs per area_ug/cm2`, color = species, fill = species)) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  #theme(axis.text.x = element_text(vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata", 
                              "Pocillopora acuta", "Porites compressa"))  +
  labs(y = "Host Carbsein (ug/cm^2)") +  # Add y-axis label
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
  x2$species <- as.factor(x2$species)

  # Perform Shapiro-Wilk test for normality for each species
  shapiro_results <- x2 %>%
    group_by(species) %>%
    summarise(
      p_value = shapiro.test(`Carbs per area_ug/cm2`)$p.value
    )
  
  # View the results
  print(shapiro_results)
  # # A tibble: 3 × 2
  # species            p_value
  # <fct>                <dbl>
  #   1 Montipora capitata   0.376
  # 2 Pocillopora acuta    0.549
  # 3 Porites compressa    0.556
  
  # Generate Q-Q plots for each species
  qqplots <- x2 %>%
    group_by(species) %>%
    nest() %>%
    mutate(
      qqplot = map(data, ~ ggqqplot(.x$`Carbs per area_ug/cm2`, 
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
  
  levene_results <- leveneTest(`Carbs per area_ug/cm2` ~ species, data = x2)
  
  # View the results
  print(levene_results)
#   Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  2  2.1732 0.1483
#       15 
  
  
  ####### ANOVA (for normally-distributed data and with homogeneous variance)
  # Perform ANOVA to test for differences between species
  anova_results <- aov(`Carbs per area_ug/cm2` ~ species, data = x2)
  
  # View the summary of the ANOVA results
  # summary(anova_results)
  # Df Sum Sq Mean Sq F value   Pr(>F)    
  # species      2 0.2567 0.12834   11.64 0.000889 ***
  #   Residuals   15 0.1654 0.01103                     
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1    

  # Perform Tukey's HSD test
  tukey_results <- TukeyHSD(anova_results)
  
  # View the results
  print(tukey_results)
  # $species
  # diff        lwr         upr     p adj
  # Pocillopora acuta-Montipora capitata -0.23634748 -0.4073123 -0.06538268 0.0070965
  # Porites compressa-Montipora capitata  0.07765587 -0.0681434  0.22345514 0.3739316
  # Porites compressa-Pocillopora acuta   0.31400335  0.1430386  0.48496815 0.0006823
  
  
  # Extract significant pairwise comparisons and transform Tukey results into the required format for 
  #ggpubr
  tukey_df <- as.data.frame(tukey_results$species)
  tukey_df$comparison <- rownames(tukey_df)
  
  significant_pairs <- tukey_df %>%
    filter(`p adj` < 0.05) %>%
    mutate(
      group1 = sub("-.*", "", comparison),  # Extract the first group
      group2 = sub(".*-", "", comparison),  # Extract the second group
      label = "***"  # Use an asterisk for significant comparisons
    ) %>%
    select(group1, group2, label)
  
  # Add significance annotations to the plot
  plot <- x2 %>%
    ggplot(aes(x = species, y = `Carbs per area_ug/cm2`, color = species, fill = species)) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("Montipora capitata", 
                                "Pocillopora acuta", "Porites compressa")) +
    labs(y = "Host Carbohydrates (ug/cm^2)") +
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
    ) 
  
  plot +
    ggpubr::stat_pvalue_manual(
      data = significant_pairs,
      label = "label",
      y.position = 0.6, 
      step.increase = 0.1
    ) 
  
  
#### save ANOVA and Tukey test results in Excel file
library(writexl)
library(broom)
library(dplyr)

# Create publication-ready ANOVA table
anova_publication_table <- tidy(anova_results) %>%
  mutate(
    # Round numeric values appropriately
    sumsq = format(sumsq, scientific = TRUE, digits = 3),
    meansq = format(meansq, scientific = TRUE, digits = 3),
    statistic = round(statistic, 3),
    p.value = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ as.character(round(p.value, 3)),
      TRUE ~ as.character(round(p.value, 3))
    ),
    # Add significance stars
    significance = case_when(
      as.numeric(gsub("< ", "", p.value)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", p.value)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", p.value)) <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rename(
    "Source" = term,
    "df" = df,
    "Sum of Squares" = sumsq,
    "Mean Square" = meansq,
    "F value" = statistic,
    "p-value" = p.value,
    "Significance" = significance
  )

# Create publication-ready Tukey HSD table
tukey_publication_table <- tukey_df %>%
  mutate(
    # Round values appropriately
    diff = round(diff, 4),
    lwr = round(lwr, 4),
    upr = round(upr, 4),
    `p adj` = case_when(
      `p adj` < 0.001 ~ "< 0.001",
      `p adj` < 0.01 ~ as.character(round(`p adj`, 4)),
      TRUE ~ as.character(round(`p adj`, 4))
    ),
    # Add significance stars
    significance = case_when(
      as.numeric(gsub("< ", "", `p adj`)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", `p adj`)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", `p adj`)) <= 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Create confidence interval column
    conf.interval = paste0("[", lwr, ", ", upr, "]")
  ) %>%
  select(-lwr, -upr) %>%
  rename(
    "Comparison" = comparison,
    "Mean Difference" = diff,
    "95% CI" = conf.interval,
    "p-value (adj)" = `p adj`,
    "Significance" = significance
  )

# Create summary statistics table
summary_stats <- x2 %>%
  group_by(species) %>%
  summarise(
    N = n(),
    Mean = round(mean(`Carbs per area_ug/cm2`, na.rm = TRUE), 4),
    SD = round(sd(`Carbs per area_ug/cm2`, na.rm = TRUE), 4),
    SE = round(sd(`Carbs per area_ug/cm2`, na.rm = TRUE) / sqrt(n()), 4),
    Min = round(min(`Carbs per area_ug/cm2`, na.rm = TRUE), 4),
    Max = round(max(`Carbs per area_ug/cm2`, na.rm = TRUE), 4),
    Median = round(median(`Carbs per area_ug/cm2`, na.rm = TRUE), 4),
    Q1 = round(quantile(`Carbs per area_ug/cm2`, 0.25, na.rm = TRUE), 4),
    Q3 = round(quantile(`Carbs per area_ug/cm2`, 0.75, na.rm = TRUE), 4),
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
    round(shapiro.test(x2$`Carbs per area_ug/cm2`[x2$species == "Montipora capitata"])$statistic, 4),
    round(shapiro.test(x2$`Carbs per area_ug/cm2`[x2$species == "Pocillopora acuta"])$statistic, 4),
    round(shapiro.test(x2$`Carbs per area_ug/cm2`[x2$species == "Porites compressa"])$statistic, 4),
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
  )
) %>%
  rename("p-value" = p_value)

# Create comprehensive Excel workbook
carbs_stats_tables <- list(
  "Summary_Statistics" = summary_stats,
  "ANOVA_Results" = anova_publication_table,
  "Tukey_HSD_Results" = tukey_publication_table,
  "Assumption_Tests" = assumption_tests,
  "Raw_Data" = x2 %>% 
    select(colony_id, species, `Carbs per area_ug/cm2`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "Carbohydrates_per_area_ug_cm2" = `Carbs per area_ug/cm2`
    )
)

# Save to Excel
write_xlsx(carbs_stats_tables, "output/Carbohydrate_Concentration_Statistical_Analysis.xlsx")

