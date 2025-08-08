

############# "GUAVA-hemocytometer comparison" #############

## load libraries
library("ggplot2")
library("plotrix")
library('dplyr')
library(tidyverse)
library(cowplot)
library(broom)
library(here)
library(readxl)
library(tidyr)
if (!require(ggpmisc)) install.packages("ggpmisc")
library(ggpmisc)

## Import metadata
guava <- read_delim(here("Hemo_Guava_trendile.csv"), delim = ";")

## Add a new column with species names
guava <- guava  %>%
  mutate(species = case_when(
    grepl("Mcap", sample_ID) ~ "Montipora capitata",
    grepl("Pacu", sample_ID) ~ "Pocillopora acuta",
    grepl("Pcom", sample_ID) ~ "Porites compressa",
    TRUE ~ NA_character_  # Default to NA if no match
  ))

# Scatter plot: Hemoc_cell/ml vs Guava_cell/ml, colored by species, with trendline and R²
ggplot(guava, aes(x = `Hemoc_cell/ml`, y = `Guava_cell/ml`)) +
  geom_point(aes(color = species), size = 3, alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    size = 0.05,
    linetype = "dashed",
    alpha = 0.2   # Lower alpha for SE shading
  ) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE
  ) +
  theme_classic() +
  labs(
    x = "Hemocytometer cells/ml",
    y = "Guava cells/ml",
    color = "Species",
    title = "Hemocytometer vs Guava Cell Counts"
  )+
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  scale_x_continuous(labels = scales::scientific) +
  scale_y_continuous(labels = scales::scientific)


############# "Algae density measurements" #############

## load libraries
library("ggplot2")
library("plotrix")
library('dplyr')
library(tidyverse)
library(cowplot)
library(broom)
library(here)
library(readxl)
library(tidyr)

## Import metadata
Data <- read_delim(here("Algae_count_guava_raw.csv"), delim = ";")

## Add a new column with species names
Data <- Data %>%
  mutate(species = case_when(
    grepl("MCAP", colony_id) ~ "Montipora capitata",
    grepl("PACT", colony_id) ~ "Pocillopora acuta",
    grepl("PCOM", colony_id) ~ "Porites compressa",
    TRUE ~ NA_character_  # Default to NA if no match
  ))


## Import surface area data
Surf_Data <- read_delim(here("1_surface_area_final.csv"), delim = ",")

## Import homogenate volumes data
HomogVol_Data <- read_delim(here("Homog_volumes.csv"), delim = ";")

## Join datasets based on Sample_ID
Combined_Data <- Data %>%
  left_join(Surf_Data, by = "colony_id") %>%
  left_join(HomogVol_Data, by = "colony_id") %>%
  distinct(colony_id, .keep_all = TRUE)

## Multiply by dilution factor
Combined_Data <- Combined_Data %>%
  mutate(`Algae_dil` = `cells/ml` * `dilution_fact`)

## Add a new column "Tot algae" to Combined_Data - normalize by homog volume
Combined_Data <- Combined_Data %>%
  mutate(`Tot_algae` = `Algae_dil` * `PBS ml`)

## Add a new column "Algae per area" to Combined_Data - normalize by skeleton surface area
Combined_Data <- Combined_Data %>%
  mutate(`Cells/cm2` = `Tot_algae` / `surface.area.cm2`)

# Select the first and 11th columns from Combined_Data
Algae_Data <- Combined_Data %>% 
  select(1, 11)

# Save the new table as a .csv file
write_csv(Algae_Data, "Algae_Data.csv")

# Plot results - algae
Combined_Data %>% 
  ggplot(aes(x = species, y = `Cells/cm2`, color = species, fill =species)) + 
  facet_wrap(~`Temp Cat`, scales = "free_y", nrow = 1) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                              "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) +
  geom_boxplot(
    aes(fill = species, fill = after_scale(colorspace::lighten(fill, .7))),
    size = 0.5, outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(width = .2, seed = 0),
    size = 2.5, alpha = .5
  ) +
  geom_point(
    position = position_jitter(width = .2, seed = 0),
    size = 2.5, stroke = .5, shape = 1, color = "white"
  )



############ Remove outliers 

#set quantile values 
q <- c(0.25, 0.75) 

Quants_algae <- Combined_Data %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`Cells/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`Cells/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`Cells/cm2`))


#Calculate Quantile upper and lower ranges
Quants_algae$upper_algae <-  Quants_algae$quant75+1.5*Quants_algae$IQRbyGroup # Upper Range
Quants_algae$lower_algae <- Quants_algae$quant25-1.5*Quants_algae$IQRbyGroup # Lower Range


# Join quantile ranges back to the original dataset
Combined_Data_with_bounds <- Combined_Data %>%
  left_join(Quants_algae, by = c("species", "Temp Cat"))


# Filter out outliers based on the quantile ranges
x1 <- Combined_Data_with_bounds %>%
  filter(`Cells/cm2` < upper_algae & `Cells/cm2` > lower_algae)

# Plot results by temp - algae
x1 %>% 
  ggplot(aes(x = species, y = `Cells/cm2`, color = species, fill =species)) + 
  facet_wrap(~`Temp Cat`, scales = "free_y", nrow = 1) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                              "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) +
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


### Plot results by species
x1$`Temp Cat`= as.factor(x1$`Temp Cat`)

x1 %>% 
  ggplot(aes(x = `Temp Cat`, y = `Cells/cm2`, color = species, fill =species)) + 
  facet_wrap(~species, scales = "free_y", nrow = 3) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                              "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) +
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
  
# Rename the column Temp Cat to Temp_Cat
x1 <- x1 %>%
  rename(Temp_Cat = `Temp Cat`)

# Ensure Temp.Cat is a factor
x1$Temp_Cat <- as.factor(x1$Temp_Cat)
x1$species <- as.factor(x1$species)


  # Perform Shapiro-Wilk test for normality for each species
  shapiro_results_algae <- x1 %>% #algae
    group_by(species, Temp_Cat) %>%
    summarise(
      p_value = shapiro.test(`Cells/cm2`)$p.value
    )
  
  # View the results
  print(shapiro_results_algae) #all normally distributed

  
  # Generate Q-Q plots for each species
  qqplots <- x1 %>%
    group_by(species, Temp_Cat) %>%
    nest() %>%
    mutate(
      qqplot = map(data, ~ ggqqplot(.x$`Cells/cm2`, 
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
  
  
  ### Perform the Levene test for homogeneity of variances
  
  library(car)
  
  # Perform Levene's test for homogeneity of variances with species as the factor and Temp_Cat as the levels
  levene_results_algae <- x1 %>%
    group_by(Temp_Cat) %>%  # Group by Temp_Cat
    summarise(
      levene_p_value = leveneTest(`Cells/cm2` ~ species, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  
  # Set the number format to display decimals instead of scientific notation
  options(scipen = 999, digits = 4)
    # View the results
  print(levene_results_algae) #12, 25, 26.8 and 35 do not pass the test

  
  levene_results_algae <- x1 %>%
    group_by(species) %>%  # Group by species
    summarise(
      levene_p_value = leveneTest(`Cells/cm2` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  
  # View the results
  print(levene_results_algae)
  # # A tibble: 3 × 2
  # species            levene_p_value
  # <fct>                       <dbl>
  #   1 Montipora capitata         0.0354
  # 2 Pocillopora acuta          0.889 
  # 3 Porites compressa          0.299 
  
  
  
  # Split into those that passed and failed the Levene test  - temp_cat
  passed_temp_cat_algae <- levene_results_algae %>%
    filter(levene_p_value > 0.05) %>%
    pull(Temp_Cat)
  
  failed_temp_cat_algae <- levene_results_algae %>%
    filter(levene_p_value <= 0.05) %>%
    pull(Temp_Cat)
  
  
  # Split species into those that passed and failed the Levene test - species
  passed_temp_cat_algae <- levene_results_algae %>%
    filter(levene_p_value > 0.05) %>%
    pull(species)
  
  failed_temp_cat_algae <- levene_results_algae %>%
    filter(levene_p_value <= 0.05) %>%
    pull(species)
  
  
  
  
  ####### ANOVA (for normally-distributed data and with homogeneous variance)
  
  # Perform ANOVA for Temp_Cat that passed the Levene test
  anova_results_algae <- x1 %>%
    filter(Temp_Cat %in% passed_temp_cat_algae) %>%
    group_by(Temp_Cat) %>%
    summarise(
      anova_p_value = summary(aov(`Cells/cm2` ~ species, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_algae)
  # # A tibble: 3 × 2
  # Temp_Cat anova_p_value
  # <fct>            <dbl>
  # 1 18          0.00000146
  # 2 30          0.00162
  

  # Perform ANOVA for species that passed the Levene test
  anova_results_algae <- x1 %>%
    filter(species %in% passed_temp_cat_algae) %>%
    group_by(species) %>%
    summarise(
      anova_p_value = summary(aov(`Cells/cm2` ~ Temp_Cat, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_algae)
  # # A tibble: 2 × 2
  # species           anova_p_value
  # <fct>                     <dbl>
  #   1 Pocillopora acuta         0.774
  # 2 Porites compressa         0.161
  
  
  # Perform Tukey's HSD test for Temp_Cat that passed the Levene test
  tukey_results_algae <- x1 %>%
    filter(Temp_Cat %in% passed_temp_cat_algae) %>%
    group_by(Temp_Cat) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`Cells/cm2` ~ species, data = cur_data()))),
      .groups = "drop"
    ) %>%
    mutate(
      tukey_table = map(tukey, ~ as.data.frame(.x$species) %>%
                          mutate(Comparison = rownames(.x$species)))
    ) %>%
    select(Temp_Cat, tukey_table) %>%
    unnest(tukey_table)
  
  # View Tukey's HSD results
  print(tukey_results_algae)
  

  # Perform Tukey's HSD test for species that passed the Levene test
  tukey_results_algae <- x1 %>%
    filter(species %in% passed_temp_cat_algae) %>%
    group_by(species) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`Cells/cm2` ~ Temp_Cat, data = cur_data()))),
      .groups = "drop"
    ) %>%
    mutate(
      tukey_table = map(tukey, ~ as.data.frame(.x$Temp_Cat) %>%
                          mutate(Comparison = rownames(.x$Temp_Cat)))
    ) %>%
    select(species, tukey_table) %>%
    unnest(tukey_table)
  
  # View Tukey's HSD results
  print(tukey_results_algae)
  # # A tibble: 30 × 6
  # species               diff      lwr     upr `p adj` Comparison
  # <fct>                <dbl>    <dbl>   <dbl>   <dbl> <chr>     
  #   1 Pocillopora acuta -121179. -450222. 207864.   0.875 18-12     
  # 2 Pocillopora acuta -134890. -423480. 153700.   0.723 25-12     
  # 3 Pocillopora acuta  -45313. -357026. 266399.   0.998 26.8-12   
  # 4 Pocillopora acuta  -82504. -381223. 216215.   0.960 30-12     
  # 5 Pocillopora acuta  -92230. -380819. 196360.   0.927 35-12     
  # 6 Pocillopora acuta  -13711. -342754. 315332.   1.00  25-18     
  # 7 Pocillopora acuta   75866. -273634. 425365.   0.986 26.8-18   
  # 8 Pocillopora acuta   38675. -299287. 376637.   0.999 30-18     
  # 9 Pocillopora acuta   28949. -300094. 357992.   1.00  35-18     
  # 10 Pocillopora acuta   89577. -222136. 401289.   0.952 26.8-25  
  
  
  
  ########## Perform Kruskal-Wallis test for Temp_Cat that failed the Levene test
  kruskal_results_algae <- x1 %>%
    filter(Temp_Cat %in% failed_temp_cat_algae) %>%
    group_by(Temp_Cat) %>%
    summarise(
      kruskal_p_value = kruskal.test(`Cells/cm2` ~ species, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_algae)
  # # A tibble: 3 × 2
  # Temp_Cat kruskal_p_value
  # <fct>              <dbl>
  # 1 12              0.000670
  # 2 25              0.000539
  # 3 26.8            0.00144 
  # 4 35              0.000377
  
  
  # Perform Kruskal-Wallis test for species that failed the Levene test
  kruskal_results_algae <- x1 %>%
    filter(species %in% failed_temp_cat_algae) %>%
    group_by(species) %>%
    summarise(
      kruskal_p_value = kruskal.test(`Cells/cm2` ~ Temp_Cat, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_algae)
  # # A tibble: 1 × 2
  # species            kruskal_p_value
  # <fct>                        <dbl>
  #   1 Montipora capitata           0.147
  
  
 
  
  ####### Perform Dunn's test for Temp_Cat that failed the Levene test
  
  # Load the FSA package
  library(FSA)
  
  dunn_results_algae <- x1 %>%
    filter(Temp_Cat %in% failed_temp_cat_algae) %>%
    group_by(Temp_Cat) %>%
    summarise(
      dunn_test = list(dunnTest(`Cells/cm2` ~ species, data = cur_data(), method = "bonferroni")),
      .groups = "drop"
    ) %>%
    mutate(
      dunn_table = map(dunn_test, ~ as.data.frame(.x$res) %>%
                         mutate(comparison = rownames(.x$res)))
    ) %>%
    select(Temp_Cat, dunn_table) %>%
    unnest(dunn_table)
  
  # View Dunn's test results
  print(dunn_results_algae)
  
  
  # Perform Dunn's test for species that failed the Levene test
  
  dunn_results_algae <- x1 %>%
    filter(species %in% failed_temp_cat_algae) %>%
    group_by(species) %>%
    summarise(
      dunn_test = list(dunnTest(`Cells/cm2` ~ Temp_Cat, data = cur_data(), method = "bonferroni")),
      .groups = "drop"
    ) %>%
    mutate(
      dunn_table = map(dunn_test, ~ as.data.frame(.x$res) %>%
                         mutate(comparison = rownames(.x$res)))
    ) %>%
    select(species, dunn_table) %>%
    unnest(dunn_table)
  
  # View Dunn's test results
  print(dunn_results_algae)
  # # A tibble: 15 × 6
  # species            Comparison       Z P.unadj P.adj comparison
  # <fct>              <chr>        <dbl>   <dbl> <dbl> <chr>     
  #   1 Montipora capitata 12 - 18     0.487  0.627   1     1         
  # 2 Montipora capitata 12 - 25    -0.0725 0.942   1     2         
  # 3 Montipora capitata 18 - 25    -0.540  0.589   1     3         
  # 4 Montipora capitata 12 - 26.8  -0.653  0.514   1     4         
  # 5 Montipora capitata 18 - 26.8  -1.08   0.280   1     5         
  # 6 Montipora capitata 25 - 26.8  -0.562  0.574   1     6         
  # 7 Montipora capitata 12 - 30    -1.65   0.0981  1     7         
  # 8 Montipora capitata 18 - 30    -2.02   0.0436  0.654 8         
  # 9 Montipora capitata 25 - 30    -1.53   0.127   1     9         
  # 10 Montipora capitata 26.8 - 30  -0.945  0.344   1     10        
  # 11 Montipora capitata 12 - 35     0.954  0.340   1     11        
  # 12 Montipora capitata 18 - 35     0.396  0.692   1     12        
  # 13 Montipora capitata 25 - 35     0.994  0.320   1     13        
  # 14 Montipora capitata 26.8 - 35   1.57   0.115   1     14        
  # 15 Montipora capitata 30 - 35     2.61   0.00911 0.137 15 
  
  
 
  # Combine Dunn and Tukey results into a single data frame
  significant_results_algae <- bind_rows(
    dunn_results_algae %>%
      filter(P.adj <= 0.05) %>%  # Filter significant results
      mutate(
        group1 = sub(" - .*", "", Comparison),  # Extract the first group
        group2 = sub(".* - ", "", Comparison),  # Extract the second group
        label = case_when(  # Assign labels based on p-value ranges
          P.adj <= 0.001 ~ "***",
          P.adj <= 0.01 ~ "**",
          P.adj <= 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      select(Temp_Cat, group1, group2, label),
    tukey_results_algae %>%
      filter(`p adj` <= 0.05) %>%  # Filter significant results
      mutate(
        group1 = sub("-.*", "", Comparison),  # Extract the first group
        group2 = sub(".*-", "", Comparison),  # Extract the second group
        label = case_when(  # Assign labels based on p-value ranges
          `p adj` <= 0.001 ~ "***",
          `p adj` <= 0.01 ~ "**",
          `p adj` <= 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      select(Temp_Cat, group1, group2, label)
  )
  

  
  
########### Add y.position and bar_y columns to significant_results_chl
  significant_results_algae <- significant_results_algae %>%
    left_join(
      x1 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`Cells/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group1" = "species")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x1 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`Cells/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group2" = "species")
    ) %>%
    rename(max_y_group2 = max_y) %>%  # Rename for clarity
    mutate(
      base_y = pmax(max_y_group1, max_y_group2) + 170000  # distance from bars and data points
    ) %>%
    group_by(Temp_Cat) %>%
    mutate(
      bar_y = base_y + (row_number() - 1) * 170000,  # Increment bar_y for each comparison within the same Temp_Cat
      y.position = bar_y + 1  # Add extra padding for asterisks
    ) %>%
    ungroup()
  
  
  
  # Ensure group1 and group2 are factors
  significant_results_algae$group1 <- as.factor(significant_results_algae$group1)
  significant_results_algae$group2 <- as.factor(significant_results_algae$group2)
  
  
  # Ensure the ggpubr package is loaded
  library(ggpubr)
  
  # Prepare the annotations for statistical significance
  annotations <- significant_results_algae %>%
    mutate(
      group1 = as.character(group1),  # Ensure group1 is a character
      group2 = as.character(group2)   # Ensure group2 is a character
    )
  
  
  plot_x1_stat <- x1 %>%
    ggplot(aes(x = species, y = `Cells/cm2`, color = species, fill = species)) +
    facet_wrap(~Temp_Cat, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                                "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) +
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
    # Add statistical annotations
    stat_pvalue_manual(
      significant_results_algae,
      label = "label",  # Use the label column for significance levels
      xmin = "group1",  # Start of the comparison
      xmax = "group2",  # End of the comparison
      y.position = significant_results_algae$y.position - 1,  # Decrease space between bars and data points
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.05)),  # Remove extra space below and slightly above
      labels = function(x) {
        sapply(x, function(y) {
          if (is.na(y) || y == 0) {  # Handle NA and 0 values
            return("0")
          }
          exponent <- floor(log10(y))
          base <- round(y / (10^exponent), 1)  # Round the base to 1 decimal place
          paste0(base, "*10^", exponent)
        })
      }  # Custom function to format y-axis labels as "base * 10^exponent"
    )
  
  # Display the updated plot
  plot_x1_stat
  
  
###by species
  
plot_x1_species_stat <- x1 %>%
    ggplot(aes(x = Temp_Cat, y = `Cells/cm2`, color = species, fill = species)) +
    facet_wrap(~species, scales = "free_y", nrow = 3) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text = element_blank()  # Remove facet titles
    ) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                                "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) +
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
  
  plot_x1_species_stat
  


######### Save statistical tests resutls for supplementary 
  
library(writexl)
library(broom)
library(dplyr)
library(purrr)
library(car) 
library(FSA)

# Separation of TEMPERATURE vs SPECIES analysis

# 1. TEMPERATURE ANALYSIS: Run Levene's test by temperature (comparing species within each temp)
levene_results_by_temp_algae <- x1 %>%
  group_by(Temp_Cat) %>%
  summarise(
    levene_p_value = leveneTest(`Cells/cm2` ~ species, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

# Extract temperature categories that passed/failed Levene's test
passed_temp_algae <- levene_results_by_temp_algae %>%
  filter(levene_p_value > 0.05) %>%
  pull(Temp_Cat)

failed_temp_algae <- levene_results_by_temp_algae %>%
  filter(levene_p_value <= 0.05) %>%
  pull(Temp_Cat)

# 2. SPECIES ANALYSIS: Run Levene's test by species (comparing temperatures within each species)
levene_results_by_species_algae <- x1 %>%
  group_by(species) %>%
  summarise(
    levene_p_value = leveneTest(`Cells/cm2` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

# Extract species that passed/failed Levene's test
passed_species_algae <- levene_results_by_species_algae %>%
  filter(levene_p_value > 0.05) %>%
  pull(species)

failed_species_algae <- levene_results_by_species_algae %>%
  filter(levene_p_value <= 0.05) %>%
  pull(species)

# 1. ANOVA by Temperature (comparing species within each temperature)
anova_results_algae_by_temp <- x1 %>%
  filter(Temp_Cat %in% passed_temp_algae) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$species)) < 2) {
      return(data.frame(
        term = "species", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`Cells/cm2` ~ species, data = .x)
      tidy_result <- tidy(anova_result)
      return(tidy_result)
    }, error = function(e) {
      return(data.frame(
        term = "species", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = as.character(e)
      ))
    })
  }) %>%
  filter(term != "Residuals") %>%
  mutate(
    variable = "Algae_density",
    analysis_by = "Temperature",
    test_type = "ANOVA"
  )

# 2. ANOVA by Species (comparing temperatures within each species)
anova_results_algae_by_species <- x1 %>%
  filter(species %in% passed_species_algae) %>%
  group_by(species) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$Temp_Cat)) < 2) {
      return(data.frame(
        term = "Temp_Cat", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`Cells/cm2` ~ Temp_Cat, data = .x)
      tidy_result <- tidy(anova_result)
      return(tidy_result)
    }, error = function(e) {
      return(data.frame(
        term = "Temp_Cat", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = as.character(e)
      ))
    })
  }) %>%
  filter(term != "Residuals") %>%
  mutate(
    variable = "Algae_density",
    analysis_by = "Species",
    test_type = "ANOVA",
    Temp_Cat = NA
  )

# COMBINE ALL ANOVA RESULTS
all_anova_results_algae <- bind_rows(
  anova_results_algae_by_temp,
  anova_results_algae_by_species
) %>%
  mutate(
    # Format results for publication
    statistic = round(statistic, 3),
    p.value = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ as.character(round(p.value, 3)),
      TRUE ~ as.character(round(p.value, 3))
    ),
    significance = case_when(
      as.numeric(gsub("< ", "", p.value)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", p.value)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", p.value)) <= 0.05 ~ "*",
      TRUE ~ ""
    ),
    sumsq = round(sumsq, 4),
    meansq = round(meansq, 4)
  ) %>%
  select(variable, analysis_by, Temp_Cat, species, term, df, sumsq, meansq, statistic, p.value, significance, test_type)

# EXTRACT TUKEY HSD RESULTS FOR PARAMETRIC ANALYSES
# Tukey for Algae by Temperature
tukey_algae_by_temp <- x1 %>%
  filter(Temp_Cat %in% passed_temp_algae) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`Cells/cm2` ~ species, data = .x)
      tukey_result <- TukeyHSD(anova_result)
      tukey_df <- as.data.frame(tukey_result$species)
      tukey_df$comparison <- rownames(tukey_df)
      return(tukey_df)
    }, error = function(e) {
      return(data.frame(
        diff = NA, lwr = NA, upr = NA, `p adj` = NA,
        comparison = "Error", error = as.character(e)
      ))
    })
  }) %>%
  mutate(
    variable = "Algae_density",
    analysis_by = "Temperature",
    test_type = "Tukey_HSD"
  )

# Tukey for Algae by Species
tukey_algae_by_species <- x1 %>%
  filter(species %in% passed_species_algae) %>%
  group_by(species) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`Cells/cm2` ~ Temp_Cat, data = .x)
      tukey_result <- TukeyHSD(anova_result)
      tukey_df <- as.data.frame(tukey_result$Temp_Cat)
      tukey_df$comparison <- rownames(tukey_df)
      return(tukey_df)
    }, error = function(e) {
      return(data.frame(
        diff = NA, lwr = NA, upr = NA, `p adj` = NA,
        comparison = "Error", error = as.character(e)
      ))
    })
  }) %>%
  mutate(
    variable = "Algae_density",
    analysis_by = "Species",
    test_type = "Tukey_HSD",
    Temp_Cat = NA
  )

# Combine all Tukey results
all_tukey_results_algae <- bind_rows(
  tukey_algae_by_temp,
  tukey_algae_by_species
) %>%
  mutate(
    # Format for publication
    diff = round(diff, 2),
    lwr = round(lwr, 2),
    upr = round(upr, 2),
    `p adj` = case_when(
      `p adj` < 0.001 ~ "< 0.001",
      `p adj` < 0.01 ~ as.character(round(`p adj`, 4)),
      TRUE ~ as.character(round(`p adj`, 4))
    ),
    significance = case_when(
      as.numeric(gsub("< ", "", `p adj`)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", `p adj`)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", `p adj`)) <= 0.05 ~ "*",
      TRUE ~ ""
    ),
    conf_interval = paste0("[", lwr, ", ", upr, "]")
  ) %>%
  select(-lwr, -upr) %>%
  rename(
    "Mean_Difference" = diff,
    "95%_CI" = conf_interval,
    "p_adj" = `p adj`,
    "Comparison" = comparison
  )

# EXTRACT KRUSKAL-WALLIS RESULTS FOR NON-PARAMETRIC ANALYSES
# Kruskal-Wallis for Algae by Temperature
kruskal_algae_by_temp <- x1 %>%
  filter(Temp_Cat %in% failed_temp_algae) %>%
  group_by(Temp_Cat) %>%
  summarise(
    statistic = kruskal.test(`Cells/cm2` ~ species, data = cur_data())$statistic,
    p.value = kruskal.test(`Cells/cm2` ~ species, data = cur_data())$p.value,
    parameter = kruskal.test(`Cells/cm2` ~ species, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Algae_density",
    analysis_by = "Temperature",
    test_type = "Kruskal-Wallis"
  )

# Kruskal-Wallis for Algae by Species - with data validation
kruskal_algae_by_species <- x1 %>%
  filter(species %in% failed_species_algae) %>%
  group_by(species) %>%
  summarise(
    n_temps = n_distinct(Temp_Cat),
    n_obs = n(),
    statistic = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                      kruskal.test(`Cells/cm2` ~ Temp_Cat, data = cur_data())$statistic, 
                      NA),
    p.value = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                    kruskal.test(`Cells/cm2` ~ Temp_Cat, data = cur_data())$p.value, 
                    NA),
    parameter = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                      kruskal.test(`Cells/cm2` ~ Temp_Cat, data = cur_data())$parameter, 
                      NA),
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Algae_density",
    analysis_by = "Species",
    test_type = "Kruskal-Wallis",
    Temp_Cat = NA,
    error = ifelse(is.na(statistic), "Insufficient groups or data", NA)
  ) %>%
  select(-n_temps, -n_obs)

# Combine all Kruskal-Wallis results
all_kruskal_results_algae <- bind_rows(
  kruskal_algae_by_temp,
  kruskal_algae_by_species
) %>%
  mutate(
    statistic = round(statistic, 3),
    p.value = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ as.character(round(p.value, 3)),
      TRUE ~ as.character(round(p.value, 3))
    ),
    significance = case_when(
      as.numeric(gsub("< ", "", p.value)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", p.value)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", p.value)) <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rename("Chi_squared" = statistic, "df" = parameter)

# 4. EXTRACT DUNN'S TEST RESULTS FOR NON-PARAMETRIC POST-HOC
# Dunn's test for Algae by Temperature
dunn_algae_by_temp <- data.frame()

for(temp in failed_temp_algae) {
  temp_data <- x1 %>% filter(Temp_Cat == temp)
  
  tryCatch({
    dunn_result <- dunnTest(`Cells/cm2` ~ species, data = temp_data, method = "bonferroni")
    
    # Extract results and add metadata
    temp_results <- dunn_result$res %>%
      mutate(
        Temp_Cat = temp,
        variable = "Algae_density",
        analysis_by = "Temperature",
        test_type = "Dunn_Test"
      )
    
    dunn_algae_by_temp <- bind_rows(dunn_algae_by_temp, temp_results)
    
  }, error = function(e) {
    cat("Error for Algae temperature", temp, ":", e$message, "\n")
    
    error_row <- data.frame(
      Comparison = "Error",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      Temp_Cat = temp,
      variable = "Algae_density",
      analysis_by = "Temperature",
      test_type = "Dunn_Test",
      error = as.character(e$message)
    )
    
    dunn_algae_by_temp <- bind_rows(dunn_algae_by_temp, error_row)
  })
}

# Dunn's test for Algae by Species - with data validation
dunn_algae_by_species <- data.frame()

for(sp in failed_species_algae) {
  species_data <- x1 %>% filter(species == sp)
  
  # Check if there are enough groups and observations
  if(n_distinct(species_data$Temp_Cat) > 1 && nrow(species_data) > 2) {
    tryCatch({
      dunn_result <- dunnTest(`Cells/cm2` ~ Temp_Cat, data = species_data, method = "bonferroni")
      
      species_results <- dunn_result$res %>%
        mutate(
          species = sp,
          variable = "Algae_density",
          analysis_by = "Species",
          test_type = "Dunn_Test",
          Temp_Cat = NA
        )
      
      dunn_algae_by_species <- bind_rows(dunn_algae_by_species, species_results)
      
    }, error = function(e) {
      cat("Error for Algae species", sp, ":", e$message, "\n")
      
      error_row <- data.frame(
        Comparison = "Error",
        Z = NA,
        P.unadj = NA,
        P.adj = NA,
        species = sp,
        variable = "Algae_density",
        analysis_by = "Species",
        test_type = "Dunn_Test",
        Temp_Cat = NA,
        error = as.character(e$message)
      )
      
      dunn_algae_by_species <- bind_rows(dunn_algae_by_species, error_row)
    })
  } else {
    # Not enough data for testing
    error_row <- data.frame(
      Comparison = "Insufficient data",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      species = sp,
      variable = "Algae_density",
      analysis_by = "Species",
      test_type = "Dunn_Test",
      Temp_Cat = NA,
      error = paste("Only", n_distinct(species_data$Temp_Cat), "temperature groups with", nrow(species_data), "observations")
    )
    
    dunn_algae_by_species <- bind_rows(dunn_algae_by_species, error_row)
  }
}

# Combine all Dunn's test results
all_dunn_results_algae <- bind_rows(
  dunn_algae_by_temp,
  dunn_algae_by_species
) %>%
  mutate(
    Z = round(as.numeric(Z), 3),
    P.unadj = case_when(
      is.na(P.unadj) ~ "Error",
      P.unadj < 0.001 ~ "< 0.001",
      P.unadj < 0.01 ~ as.character(round(P.unadj, 4)),
      TRUE ~ as.character(round(P.unadj, 4))
    ),
    P.adj = case_when(
      is.na(P.adj) ~ "Error",
      P.adj < 0.001 ~ "< 0.001",
      P.adj < 0.01 ~ as.character(round(P.adj, 4)),
      TRUE ~ as.character(round(P.adj, 4))
    ),
    significance = case_when(
      P.adj == "Error" ~ "",
      as.numeric(gsub("< ", "", P.adj)) <= 0.001 ~ "***",
      as.numeric(gsub("< ", "", P.adj)) <= 0.01 ~ "**",
      as.numeric(gsub("< ", "", P.adj)) <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rename(
    "Z_statistic" = Z,
    "p_unadj" = P.unadj,
    "p_adj" = P.adj
  )

# Print the final combined results
print(all_dunn_results_algae)

# 5. CREATE SUMMARY STATISTICS
summary_stats_algae <- x1 %>%
  group_by(species, Temp_Cat) %>%
  summarise(
    N = n(),
    Mean = round(mean(`Cells/cm2`, na.rm = TRUE), 2),
    SD = round(sd(`Cells/cm2`, na.rm = TRUE), 2),
    SE = round(sd(`Cells/cm2`, na.rm = TRUE) / sqrt(n()), 2),
    Median = round(median(`Cells/cm2`, na.rm = TRUE), 2),
    IQR = round(IQR(`Cells/cm2`, na.rm = TRUE), 2),
    Min = round(min(`Cells/cm2`, na.rm = TRUE), 2),
    Max = round(max(`Cells/cm2`, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(Variable = "Algae_density")

# Create methods summary
statistical_methods_summary_algae <- data.frame(
  Variable = "Algae_density",
  Units = "Cells per cm²",
  Normality_Test = "Shapiro-Wilk",
  Homogeneity_Test = "Levene's test",
  Parametric_Test = "ANOVA + Tukey HSD",
  NonParametric_Test = "Kruskal-Wallis + Dunn's test",
  Significance_Level = "α = 0.05",
  Multiple_Comparisons = "Bonferroni correction",
  Outlier_Removal = "1.5 × IQR method"
)

# Create final Excel workbook
algae_density_stats_tables <- list(
  "Summary_Statistics" = summary_stats_algae,
  "ANOVA_Results" = all_anova_results_algae,
  "Tukey_HSD_Results" = all_tukey_results_algae,
  "Kruskal_Wallis_Results" = all_kruskal_results_algae,
  "Dunn_PostHoc_Results" = all_dunn_results_algae,
  "Statistical_Methods_Used" = statistical_methods_summary_algae,
  "Shapiro_Wilk_Results" = shapiro_results_algae %>%
    rename("Shapiro_p_value" = p_value) %>%
    mutate(
      Normal_Distribution = ifelse(Shapiro_p_value > 0.05, "Yes", "No"),
      Variable = "Algae_density"
    ),
  "Levene_Test_Results" = bind_rows(
    levene_results_by_temp_algae %>% mutate(Variable = "Algae_density", Analysis = "by_Temperature"),
    levene_results_by_species_algae %>% mutate(Variable = "Algae_density", Analysis = "by_Species")
  ),
  "Raw_Data_Algae_Density" = x1 %>% 
    select(colony_id, species, Temp_Cat, `Cells/cm2`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "Temperature" = Temp_Cat,
      "Algae_Cells_per_cm2" = `Cells/cm2`
    )
)

# Create output directory and save
write_xlsx(algae_density_stats_tables, "Algae_Density_Statistical_Analysis.xlsx")
