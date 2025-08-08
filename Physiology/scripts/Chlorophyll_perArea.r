
############# "Chlorophyll a and c concentration measurements" #############

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
Data <- read_delim(here("Chlorophyll_raw.csv"), delim = ";")
Data_algae <- read_delim(here("Algae_Data.csv"), delim = ",")

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

## Add a new column "Tot Chla ug" to Combined_Data - normalize by homog volume
Combined_Data <- Combined_Data %>%
  mutate(`Tot chla ug` = `Chla_ug/ml` * `PBS ml`)

## Add a new column "Tot Chlc ug" to Combined_Data - normalize by homog volume
Combined_Data <- Combined_Data %>%
  mutate(`Tot chlc ug` = `Chlc_ug/ml` * `PBS ml`)

## Add a new column "Chla per area" to Combined_Data - normalize by skeleton surface area
Combined_Data <- Combined_Data %>%
  mutate(`chla per area_ug/cm2` = `Tot chla ug` / `surface.area.cm2`)

## Add a new column "Chlc per area" to Combined_Data - normalize by skeleton surface area
Combined_Data <- Combined_Data %>%
  mutate(`chlc per area_ug/cm2` = `Tot chlc ug` / `surface.area.cm2`)

# Plot results by temp - Chla
Combined_Data %>% 
  ggplot(aes(x = species, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
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

# Plot results by temp - Chlc
Combined_Data %>% 
  ggplot(aes(x = species, y = `chlc per area_ug/cm2`, color = species, fill =species)) + 
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



# Plot results by species - Chla

Combined_Data$`Temp Cat`= as.factor(Combined_Data$`Temp Cat`)

Combined_Data %>% 
  ggplot(aes(x = `Temp Cat`, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
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

# Plot results by species - Chlc
Combined_Data %>% 
  ggplot(aes(x = `Temp Cat`, y = `chlc per area_ug/cm2`, color = species, fill =species)) + 
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

Quants_chla <- Combined_Data %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chla per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chla per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`chla per area_ug/cm2`))

Quants_chlc <- Combined_Data %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chlc per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chlc per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`chlc per area_ug/cm2`))


#Calculate Quantile upper and lower ranges
Quants_chla$upper_chla <-  Quants_chla$quant75+1.5*Quants_chla$IQRbyGroup # Upper Range
Quants_chla$lower_chla <- Quants_chla$quant25-1.5*Quants_chla$IQRbyGroup # Lower Range

Quants_chlc$upper_chlc <-  Quants_chlc$quant75+1.5*Quants_chlc$IQRbyGroup # Upper Range
Quants_chlc$lower_chlc <- Quants_chlc$quant25-1.5*Quants_chlc$IQRbyGroup # Lower Range


# Join quantile ranges back to the original dataset
Combined_Data_with_bounds <- Combined_Data %>%
  left_join(Quants_chla, by = c("species", "Temp Cat"))

Combined_Data_with_bounds <- Combined_Data_with_bounds %>%
  left_join(Quants_chlc, by = c("species", "Temp Cat"))

# Filter out outliers based on the quantile ranges
x1 <- Combined_Data_with_bounds %>%
  filter(`chla per area_ug/cm2` < upper_chla & `chla per area_ug/cm2` > lower_chla)

x2 <- Combined_Data_with_bounds %>%
  filter(`chlc per area_ug/cm2` < upper_chlc & `chlc per area_ug/cm2` > lower_chlc)

# Plot results by temp - chla
x1 %>% 
  ggplot(aes(x = species, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
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

# Plot results by temp - chlc
x2 %>% 
  ggplot(aes(x = species, y = `chlc per area_ug/cm2`, color = species, fill =species)) + 
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

# Plot results by species - chla

x1$`Temp Cat`= as.factor(x1$`Temp Cat`)

x1 %>% 
  ggplot(aes(x = `Temp Cat`, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
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

# Plot results by species - chlc

x2$`Temp Cat`= as.factor(x2$`Temp Cat`)

x2 %>% 
  ggplot(aes(x = `Temp Cat`, y = `chlc per area_ug/cm2`, color = species, fill =species)) + 
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


#### Remove outliers round 2 (outliers after 1st round are re-calculated based on the samples left)

# Create a new table x3 with selected columns from x1
x3 <- x1 %>% select(1, 4, 6, 10)

x4 <- x2 %>% select(1, 4, 6, 11)


Quants_chla <- x3 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chla per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chla per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`chla per area_ug/cm2`))

Quants_chlc <- x4 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chlc per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chlc per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`chlc per area_ug/cm2`))


#Calculate Quantile upper and lower ranges
Quants_chla$upper_chla <-  Quants_chla$quant75+1.5*Quants_chla$IQRbyGroup # Upper Range
Quants_chla$lower_chla <- Quants_chla$quant25-1.5*Quants_chla$IQRbyGroup # Lower Range

Quants_chlc$upper_chlc <-  Quants_chlc$quant75+1.5*Quants_chlc$IQRbyGroup # Upper Range
Quants_chlc$lower_chlc <- Quants_chlc$quant25-1.5*Quants_chlc$IQRbyGroup # Lower Range


# Join quantile ranges back to the original dataset
Chla_with_bounds <- x3 %>%
  left_join(Quants_chla, by = c("species", "Temp Cat"))

Chlc_with_bounds <- x4 %>%
  left_join(Quants_chlc, by = c("species", "Temp Cat"))

# Filter out outliers based on the quantile ranges
x5 <- Chla_with_bounds %>%
  filter(`chla per area_ug/cm2` < upper_chla & `chla per area_ug/cm2` > lower_chla)

x6 <- Chlc_with_bounds %>%
  filter(`chlc per area_ug/cm2` < upper_chlc & `chlc per area_ug/cm2` > lower_chlc)

# Plot results - chla
x5 %>% 
  ggplot(aes(x = species, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
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


# Plot results by temp - chlc
plot_x6 <-  x6 %>% 
  ggplot(aes(x = species, y = `chlc per area_ug/cm2`, color = species, fill =species)) + 
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


## by species

x6$Temp_Cat= as.factor(x6$Temp_Cat)

plot_x6 <-  x6 %>% 
  ggplot(aes(x = Temp_Cat, y = `chlc per area_ug/cm2`, color = species, fill =species)) + 
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


plot_x6


#### Remove outliers round 3 - chla

# Create a new table x3 with selected columns from x1
x7 <- x5 %>% select(1, 2,3,4)


Quants_chla <- x7 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chla per area_ug/cm2`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chla per area_ug/cm2`, probs = q[2]),
            IQRbyGroup=IQR(`chla per area_ug/cm2`))

#Calculate Quantile upper and lower ranges
Quants_chla$upper_chla <-  Quants_chla$quant75+1.5*Quants_chla$IQRbyGroup # Upper Range
Quants_chla$lower_chla <- Quants_chla$quant25-1.5*Quants_chla$IQRbyGroup # Lower Range


# Join quantile ranges back to the original dataset
Chla_with_bounds <- x7 %>%
  left_join(Quants_chla, by = c("species", "Temp Cat"))


# Filter out outliers based on the quantile ranges
x8 <- Chla_with_bounds %>%
  filter(`chla per area_ug/cm2` < upper_chla & `chla per area_ug/cm2` > lower_chla)


# Plot results - chla
plot_x8 <- x8 %>% 
  ggplot(aes(x = species, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
  facet_wrap(~`Temp Cat`, scales = "free_y", nrow = 1) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.title.x = element_blank()) + 
  # scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
  #                             "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) +
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


## by species

x8$Temp_Cat= as.factor(x8$Temp_Cat)

plot_x8 <-  x8 %>% 
  ggplot(aes(x = Temp_Cat, y = `chla per area_ug/cm2`, color = species, fill =species)) + 
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


plot_x8


# Load the patchwork package
library(patchwork)

# Combine the plots into a single figure
combined_plot <- plot_x8 / plot_x6

# Display the combined plot
combined_plot


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
x8 <- x8 %>%
  rename(Temp_Cat = `Temp Cat`)

x6 <- x6 %>%
  rename(Temp_Cat = `Temp Cat`)

  # Ensure Temp.Cat is a factor
x8$Temp_Cat <- as.factor(x8$Temp_Cat)
x6$Temp_Cat <- as.factor(x6$Temp_Cat)
x8$species <- as.factor(x8$species)
x6$species <- as.factor(x6$species)

  # Perform Shapiro-Wilk test for normality for each species
  shapiro_results_chla <- x8 %>% #chla
    group_by(species, Temp_Cat) %>%
    summarise(
      p_value = shapiro.test(`chla per area_ug/cm2`)$p.value
    )
  
  # View the results
  print(shapiro_results_chla) #all normally distributed

  shapiro_results_chlc <- x6 %>% #chlc
    group_by(species, Temp_Cat) %>%
    summarise(
      p_value = shapiro.test(`chlc per area_ug/cm2`)$p.value
    )
  
  # View the results
  print(shapiro_results_chlc) #all normally distributed
  
    # Generate Q-Q plots for each species
  qqplots <- x8 %>%
    group_by(species, Temp_Cat) %>%
    nest() %>%
    mutate(
      qqplot = map(data, ~ ggqqplot(.x$`chla per area_ug/cm2`, 
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
  
  ###### Perform the Levene test for homogeneity of variances
  
  library(car)
  
  # Perform Levene's test for homogeneity of variances 
  levene_results_chla <- x8 %>%
    group_by(Temp_Cat) %>%  # Group by Temp_Cat
    summarise(
      levene_p_value = leveneTest(`chla per area_ug/cm2` ~ species, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  
    # View the results
  print(levene_results_chla) #12, 26.8 and 30 do not pass the test

  
  levene_results_chlc <- x6 %>%
    group_by(Temp_Cat) %>%  # Group by Temp_Cat
    summarise(
      levene_p_value = leveneTest(`chlc per area_ug/cm2` ~ species, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  # Set the number format to display decimals instead of scientific notation
  options(scipen = 999, digits = 4)
    # View the results
  print(levene_results_chlc) #12, 26.8 and 35 do not pass the test
  
  
  levene_results_chla <- x8 %>%
    group_by(species) %>%  # Group by species
    summarise(
      levene_p_value = leveneTest(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  
  # View the results
  print(levene_results_chla)
  # A tibble: 3 × 2
  # species            levene_p_value
  # <fct>                       <dbl>
  # 1 Montipora capitata         0.0130
  # 2 Pocillopora acuta          0.381 
  # 3 Porites compressa          0.0799
  
  
  levene_results_chlc <- x6 %>%
    group_by(species) %>%  # Group by species
    summarise(
      levene_p_value = leveneTest(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
 
  # View the results
  print(levene_results_chlc) #12, 26.8 and 35 do not pass the test
  # A tibble: 3 × 2
  # species            levene_p_value
  # <fct>                       <dbl>
  #   1 Montipora capitata      0.151    
  # 2 Pocillopora acuta       0.540    
  # 3 Porites compressa       0.0000895
  
  
  # Split Temp_Cat into those that passed and failed the Levene test
  passed_temp_cat_chla <- levene_results_chla %>%
    filter(levene_p_value > 0.05) %>%
    pull(Temp_Cat)
  
  failed_temp_cat_chla <- levene_results_chla %>%
    filter(levene_p_value <= 0.05) %>%
    pull(Temp_Cat)
  
  passed_temp_cat_chlc <- levene_results_chlc %>%
    filter(levene_p_value > 0.05) %>%
    pull(Temp_Cat)
  
  failed_temp_cat_chlc <- levene_results_chlc %>%
    filter(levene_p_value <= 0.05) %>%
    pull(Temp_Cat)
  
  # Split species into those that passed and failed the Levene test
  passed_temp_cat_chla <- levene_results_chla %>%
    filter(levene_p_value > 0.05) %>%
    pull(species)
  
  failed_temp_cat_chla <- levene_results_chla %>%
    filter(levene_p_value <= 0.05) %>%
    pull(species)
  
  passed_temp_cat_chlc <- levene_results_chlc %>%
    filter(levene_p_value > 0.05) %>%
    pull(species)
  
  failed_temp_cat_chlc <- levene_results_chlc %>%
    filter(levene_p_value <= 0.05) %>%
    pull(species)
  
  
  ####### ANOVA (for normally-distributed data and with homogeneous variance)
  
  # Perform ANOVA for Temp_Cat that passed the Levene test
  anova_results_chla <- x8 %>%
    filter(Temp_Cat %in% passed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      anova_p_value = summary(aov(`chla per area_ug/cm2` ~ species, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chla)
  # # A tibble: 3 × 2
  # Temp_Cat anova_p_value
  # <fct>            <dbl>
  #   1 18         0.000000437
  # 2 25         0.000000177
  # 3 35         0.00000281 
  
  
  anova_results_chlc <- x6 %>%
    filter(Temp_Cat %in% passed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      anova_p_value = summary(aov(`chlc per area_ug/cm2` ~ species, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chlc)
  # # A tibble: 3 × 2
  # Temp_Cat anova_p_value
  # <fct>            <dbl>
  # 1 18       0.00000810   
  # 2 25       0.00000000435
  # 3 30       0.00000000799
  
  
  ### Perform ANOVA for species that passed the Levene test
  anova_results_chla <- x8 %>%
    filter(species %in% passed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      anova_p_value = summary(aov(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chla)
  # A tibble: 2 × 2
  # species           anova_p_value
  # <fct>                     <dbl>
  #   1 Pocillopora acuta         0.508
  # 2 Porites compressa         0.223
  
  
  anova_results_chlc <- x6 %>%
    filter(species %in% passed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      anova_p_value = summary(aov(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chlc)
  # # A tibble: 2 × 2
  # species            anova_p_value
  # <fct>                      <dbl>
  #   1 Montipora capitata       0.00382
  # 2 Pocillopora acuta        0.150  
  
  
  
  # Perform Tukey's HSD test for Temp_Cat that passed the Levene test
  tukey_results_chla <- x8 %>%
    filter(Temp_Cat %in% passed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chla per area_ug/cm2` ~ species, data = cur_data()))),
      .groups = "drop"
    ) %>%
    mutate(
      tukey_table = map(tukey, ~ as.data.frame(.x$species) %>%
                          mutate(Comparison = rownames(.x$species)))
    ) %>%
    select(Temp_Cat, tukey_table) %>%
    unnest(tukey_table)
  
  # View Tukey's HSD results
  print(tukey_results_chla)
  
  
  tukey_results_chlc <- x6 %>%
    filter(Temp_Cat %in% passed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chlc per area_ug/cm2` ~ species, data = cur_data()))),
      .groups = "drop"
    ) %>%
    mutate(
      tukey_table = map(tukey, ~ as.data.frame(.x$species) %>%
                          mutate(Comparison = rownames(.x$species)))
    ) %>%
    select(Temp_Cat, tukey_table) %>%
    unnest(tukey_table)
  
  # View Tukey's HSD results
  print(tukey_results_chlc)
  
  
  
  # Perform Tukey's HSD test for species that passed the Levene test
  tukey_results_chla <- x8 %>%
    filter(species%in% passed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data()))),
      .groups = "drop"
    ) %>%
    mutate(
      tukey_table = map(tukey, ~ as.data.frame(.x$Temp_Cat ) %>%
                          mutate(Comparison = rownames(.x$Temp_Cat )))
    ) %>%
    select(species, tukey_table) %>%
    unnest(tukey_table)
  
  # View Tukey's HSD results
  print(tukey_results_chla)
  
  
  tukey_results_chlc <- x6 %>%
    filter(species %in% passed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data()))),
      .groups = "drop"
    ) %>%
    mutate(
      tukey_table = map(tukey, ~ as.data.frame(.x$Temp_Cat ) %>%
                          mutate(Comparison = rownames(.x$Temp_Cat )))
    ) %>%
    select(species, tukey_table) %>%
    unnest(tukey_table)
  
  # View Tukey's HSD results
  print(tukey_results_chlc)
  
  
  
  ########## Perform Kruskal-Wallis test for Temp_Cat that failed the Levene test
  kruskal_results_chla <- x8 %>%
    filter(Temp_Cat %in% failed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chla per area_ug/cm2` ~ species, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chla)
  # # A tibble: 3 × 2
  # Temp_Cat kruskal_p_value
  # <fct>              <dbl>
  #   1 12              0.00656 
  # 2 26.8            0.000477
  # 3 30              0.000685
  
  kruskal_results_chlc <- x6 %>%
    filter(Temp_Cat %in% failed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chlc per area_ug/cm2` ~ species, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chlc)
  # # A tibble: 3 × 2
  # Temp_Cat kruskal_p_value
  # <fct>              <dbl>
  #   1 12              0.00878 
  # 2 26.8            0.000642
  # 3 35              0.000458
  
  
  ###by species
  kruskal_results_chla <- x8 %>%
    filter(species %in% failed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chla)
  # # A tibble: 1 × 2
  # species            kruskal_p_value
  # <fct>                        <dbl>
  #   1 Montipora capitata         0.00271
  
  kruskal_results_chlc <- x6 %>%
    filter(species %in% failed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chlc)
  # # A tibble: 1 × 2
  # species           kruskal_p_value
  # <fct>                       <dbl>
  #   1 Porites compressa           0.207
  
  
  
  
  # Load the FSA package
  library(FSA)
  
  
  # Perform Dunn's test for Temp_Cat that failed the Levene test
  dunn_results_chla <- x8 %>%
    filter(Temp_Cat %in% failed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      dunn_test = list(dunnTest(`chla per area_ug/cm2` ~ species, data = cur_data(), method = "bonferroni")),
      .groups = "drop"
    ) %>%
    mutate(
      dunn_table = map(dunn_test, ~ as.data.frame(.x$res) %>%
                         mutate(comparison = rownames(.x$res)))
    ) %>%
    select(Temp_Cat, dunn_table) %>%
    unnest(dunn_table)
  
  # View Dunn's test results
  print(dunn_results_chla)
  
 
  dunn_results_chlc <- x6 %>%
    filter(Temp_Cat %in% failed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      dunn_test = list(dunnTest(`chlc per area_ug/cm2` ~ species, data = cur_data(), method = "bonferroni")),
      .groups = "drop"
    ) %>%
    mutate(
      dunn_table = map(dunn_test, ~ as.data.frame(.x$res) %>%
                         mutate(comparison = rownames(.x$res)))
    ) %>%
    select(Temp_Cat, dunn_table) %>%
    unnest(dunn_table)
  
  # View Dunn's test results
  print(dunn_results_chlc)
  
  
  ### by species
  dunn_results_chla <- x8 %>%
    filter(species%in% failed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      dunn_test = list(dunnTest(`chla per area_ug/cm2` ~ Temp_Cat , data = cur_data(), method = "bonferroni")),
      .groups = "drop"
    ) %>%
    mutate(
      dunn_table = map(dunn_test, ~ as.data.frame(.x$res) %>%
                         mutate(comparison = rownames(.x$res)))
    ) %>%
    select(species, dunn_table) %>%
    unnest(dunn_table)
  
  # View Dunn's test results
  print(dunn_results_chla)
  
  
  dunn_results_chlc <- x6 %>%
    filter(species %in% failed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      dunn_test = list(dunnTest(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data(), method = "bonferroni")),
      .groups = "drop"
    ) %>%
    mutate(
      dunn_table = map(dunn_test, ~ as.data.frame(.x$res) %>%
                         mutate(comparison = rownames(.x$res)))
    ) %>%
    select(species, dunn_table) %>%
    unnest(dunn_table)
  
  # View Dunn's test results
  print(dunn_results_chlc)
  
  
  

  # Combine Dunn and Tukey results into a single data frame
  significant_results_chla <- bind_rows(
    dunn_results_chla %>%
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
    tukey_results_chla %>%
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
  
  #by species
  significant_results_chla <- bind_rows(
    dunn_results_chla %>%
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
      select(species, group1, group2, label),
    tukey_results_chla %>%
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
      select(species, group1, group2, label)
  )
  
  
  #by species
  significant_results_chlc <- bind_rows(
    dunn_results_chlc %>%
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
      select(species, group1, group2, label),
    tukey_results_chlc %>%
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
      select(species, group1, group2, label)
  )
  
  

  # Add y.position and bar_y columns to significant_results_chl
  significant_results_chla <- significant_results_chla %>%
    left_join(
      x8 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group1" = "species")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x8 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group2" = "species")
    ) %>%
    rename(max_y_group2 = max_y) %>%  # Rename for clarity
    mutate(
      base_y = pmax(max_y_group1, max_y_group2) + 1  # Use the higher max_y of the two groups and add padding
    ) %>%
    group_by(Temp_Cat) %>%
    mutate(
      bar_y = base_y + (row_number() - 1) * 1.5,  # Increment bar_y for each comparison within the same Temp_Cat
      y.position = bar_y + 1  # Add extra padding for asterisks
    ) %>%
    ungroup()
  
  # Ensure group1 and group2 are factors
  significant_results_chla$group1 <- as.factor(significant_results_chla$group1)
  significant_results_chla$group2 <- as.factor(significant_results_chla$group2)
  
  
  #chlc
  significant_results_chlc <- significant_results_chlc %>%
    left_join(
      x6 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group1" = "species")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x6 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group2" = "species")
    ) %>%
    rename(max_y_group2 = max_y) %>%  # Rename for clarity
    mutate(
      base_y = pmax(max_y_group1, max_y_group2) + 0.8  # distance from bars and data points
    ) %>%
    group_by(Temp_Cat) %>%
    mutate(
      bar_y = base_y + (row_number() - 1) * 0.4,  # Increment bar_y for each comparison within the same Temp_Cat
      y.position = bar_y + 0.5  # Add extra padding for asterisks
    ) %>%
    ungroup()
  
  # Ensure group1 and group2 are factors
  significant_results_chlc$group1 <- as.factor(significant_results_chlc$group1)
  significant_results_chlc$group2 <- as.factor(significant_results_chlc$group2)
  
  
  
  #### by species
  significant_results_chla <- significant_results_chla %>%
    left_join(
      x8 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group1" = "Temp_Cat")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x8 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group2" = "Temp_Cat")
    ) %>%
    rename(max_y_group2 = max_y) %>%  # Rename for clarity
    mutate(
      base_y = pmax(max_y_group1, max_y_group2) + 1  # Use the higher max_y of the two groups and add padding
    ) %>%
    group_by(species) %>%
    mutate(
      bar_y = base_y + (row_number() - 1) * 1.5,  # Increment bar_y for each comparison within the same species
      y.position = bar_y + 1  # Add extra padding for asterisks
    ) %>%
    ungroup()
  
  # Ensure group1 and group2 are factors
  significant_results_chla$group1 <- as.factor(significant_results_chla$group1)
  significant_results_chla$group2 <- as.factor(significant_results_chla$group2)
  
  
  #chlc
  significant_results_chlc <- significant_results_chlc %>%
    left_join(
      x6 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group1" = "Temp_Cat")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x6 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc per area_ug/cm2`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group2" = "Temp_Cat")
    ) %>%
    rename(max_y_group2 = max_y) %>%  # Rename for clarity
    mutate(
      base_y = pmax(max_y_group1, max_y_group2) + 1  # Use the higher max_y of the two groups and add padding
    ) %>%
    group_by(species) %>%
    mutate(
      bar_y = base_y + (row_number() - 1) * 1.5,  # Increment bar_y for each comparison within the same species
      y.position = bar_y + 1  # Add extra padding for asterisks
    ) %>%
    ungroup()
  
  # Ensure group1 and group2 are factors
  significant_results_chlc$group1 <- as.factor(significant_results_chlc$group1)
  significant_results_chlc$group2 <- as.factor(significant_results_chlc$group2)
  
  
  
  
  
  
  # Ensure the ggpubr package is loaded
  library(ggpubr)
  
  # Prepare the annotations for statistical significance
  annotations <- significant_results_chla %>%
    mutate(
      group1 = as.character(group1),  # Ensure group1 is a character
      group2 = as.character(group2)   # Ensure group2 is a character
    )
  

  annotations <- significant_results_chlc %>%
    mutate(
      group1 = as.character(group1),  # Ensure group1 is a character
      group2 = as.character(group2)   # Ensure group2 is a character
    )
  
  
  plot_x8_stat <- x8 %>%
    ggplot(aes(x = species, y = `chla per area_ug/cm2`, color = species, fill = species)) +
    facet_wrap(~Temp_Cat, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
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
      significant_results_chla,
      label = "label",  # Use the label column for significance levels
      xmin = "group1",  # Start of the comparison
      xmax = "group2",  # End of the comparison
      y.position = significant_results_chla$y.position - 1,  # Decrease space between bars and data points
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )
  
  # Display the updated plot
  plot_x8_stat
  
  
  
  
  plot_x6_stat <- x6 %>%
    ggplot(aes(x = species, y = `chlc per area_ug/cm2`, color = species, fill = species)) +
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
      significant_results_chlc,
      label = "label",  # Use the label column for significance levels
      xmin = "group1",  # Start of the comparison
      xmax = "group2",  # End of the comparison
      y.position = significant_results_chlc$y.position - 1,  # Decrease space between bars and data points
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )
  # scale_y_continuous(
  #   limits = c(0, max(x6$`chlc per area_ug/cm2`, na.rm = TRUE) + 0.7),  # Dynamic upper limit
  #   expand = expansion(mult = c(0, 0.05))  # Remove extra space below and slightly above
  # )
  # 
  
  # Display the updated plot
  plot_x6_stat
  
  
  # Load the patchwork package
  library(patchwork)
  
  # Combine the plots into a single figure
  combined_plot_stat <- plot_x8_stat / plot_x6_stat
  
  # Display the combined plot
  combined_plot_stat
  
  
  
  ####### by species
  
 plot_x8_stat <- x8 %>%
    ggplot(aes(x = Temp_Cat, y = `chla per area_ug/cm2`, color = species, fill = species)) +
    facet_wrap(~species, scales = "free_y", nrow = 3) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
   theme_classic() +
   theme(
     legend.position = "top",
     axis.text.x = element_blank(),
     axis.title.x = element_blank(),
     strip.text = element_blank()  # Remove facet titles
   )  +
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
      significant_results_chla,
      label = "label",  # Use the label column for significance levels
      xmin = "group1",  # Start of the comparison
      xmax = "group2",  # End of the comparison
      y.position = significant_results_chla$y.position - 1,  # Decrease space between bars and data points
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )
  
  # Display the updated plot
  plot_x8_stat
  
  
  
  
plot_x6_stat <- x6 %>%
    ggplot(aes(x = Temp_Cat, y = `chlc per area_ug/cm2`, color = species, fill = species)) +
    facet_wrap(~species, scales = "free_y", nrow = 3) +
    scale_color_manual(values = c("green", "cyan", "orange")) +
    scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank()  # Remove facet titles
  ) +
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
      significant_results_chlc,
      label = "label",  # Use the label column for significance levels
      xmin = "group1",  # Start of the comparison
      xmax = "group2",  # End of the comparison
      y.position = significant_results_chlc$y.position - 1,  # Decrease space between bars and data points
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )
  # scale_y_continuous(
  #   limits = c(0, max(x6$`chlc per area_ug/cm2`, na.rm = TRUE) + 0.7),  # Dynamic upper limit
  #   expand = expansion(mult = c(0, 0.05))  # Remove extra space below and slightly above
  # )
  # 
  
  # Display the updated plot
  plot_x6_stat

  
  # Load the patchwork package
  library(patchwork)
  
  # Combine the plots into a single figure
  combined_plot_stat <- plot_x8_stat / plot_x6_stat
  
  # Display the combined plot
  combined_plot_stat
  
  
  
######### Save statistical tests resutls for supplementary 
  
library(writexl)
library(broom)
library(dplyr)
library(purrr)
library(car) 

# Check what you actually have
cat("passed_temp_cat_chla contains:", paste(passed_temp_cat_chla, collapse = ", "), "\n")
cat("failed_temp_cat_chla contains:", paste(failed_temp_cat_chla, collapse = ", "), "\n")

# For SPECIES analysis (comparing temperatures within each species):
passed_species_chla <- passed_temp_cat_chla  
failed_species_chla <- failed_temp_cat_chla  

passed_species_chlc <- passed_temp_cat_chlc   
failed_species_chlc <- failed_temp_cat_chlc  

# For TEMPERATURE analysis, need to extract the correct values
# re-run the Levene's test by temperature to get the correct categories:

# Levene's test by Temperature Category (comparing species within each temp)
levene_results_by_temp_chla <- x8 %>%
  group_by(Temp_Cat) %>%
  summarise(
    levene_p_value = leveneTest(`chla per area_ug/cm2` ~ species, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

levene_results_by_temp_chlc <- x6 %>%
  group_by(Temp_Cat) %>%
  summarise(
    levene_p_value = leveneTest(`chlc per area_ug/cm2` ~ species, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

# get the correct temperature categories
passed_temp_chla <- levene_results_by_temp_chla %>%
  filter(levene_p_value > 0.05) %>%
  pull(Temp_Cat)

failed_temp_chla <- levene_results_by_temp_chla %>%
  filter(levene_p_value <= 0.05) %>%
  pull(Temp_Cat)

passed_temp_chlc <- levene_results_by_temp_chlc %>%
  filter(levene_p_value > 0.05) %>%
  pull(Temp_Cat)

failed_temp_chlc <- levene_results_by_temp_chlc %>%
  filter(levene_p_value <= 0.05) %>%
  pull(Temp_Cat)

cat("Correct temperature categories:\n")
cat("Passed temperatures for Chla:", paste(passed_temp_chla, collapse = ", "), "\n")
cat("Failed temperatures for Chla:", paste(failed_temp_chla, collapse = ", "), "\n")

# EXTRACT ANOVA RESULTS
# 1. ANOVA by Temperature (comparing species within each temperature)
anova_results_chla_by_temp <- x8 %>%
  filter(Temp_Cat %in% passed_temp_chla) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$species)) < 2) {
      return(data.frame(
        term = "species", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chla per area_ug/cm2` ~ species, data = .x)
      tidy_result <- tidy(anova_result)
      return(tidy_result)
    }, error = function(e) {
      return(data.frame(
        term = "species", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = as.character(e)
      ))
    })
  }) %>%
  filter(term != "Residuals") %>%
  mutate(
    variable = "Chlorophyll_a",
    analysis_by = "Temperature",
    test_type = "ANOVA"
  )

# 2. ANOVA by Species (comparing temperatures within each species)
anova_results_chla_by_species <- x8 %>%
  filter(species %in% passed_species_chla) %>%
  group_by(species) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$Temp_Cat)) < 2) {
      return(data.frame(
        term = "Temp_Cat", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chla per area_ug/cm2` ~ Temp_Cat, data = .x)
      tidy_result <- tidy(anova_result)
      return(tidy_result)
    }, error = function(e) {
      return(data.frame(
        term = "Temp_Cat", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = as.character(e)
      ))
    })
  }) %>%
  filter(term != "Residuals") %>%
  mutate(
    variable = "Chlorophyll_a",
    analysis_by = "Species",
    test_type = "ANOVA",
    Temp_Cat = NA  # Add this column for consistency
  )

# Similar for Chlorophyll c
anova_results_chlc_by_temp <- x6 %>%
  filter(Temp_Cat %in% passed_temp_chlc) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$species)) < 2) {
      return(data.frame(
        term = "species", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chlc per area_ug/cm2` ~ species, data = .x)
      tidy_result <- tidy(anova_result)
      return(tidy_result)
    }, error = function(e) {
      return(data.frame(
        term = "species", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = as.character(e)
      ))
    })
  }) %>%
  filter(term != "Residuals") %>%
  mutate(
    variable = "Chlorophyll_c",
    analysis_by = "Temperature",
    test_type = "ANOVA"
  )

anova_results_chlc_by_species <- x6 %>%
  filter(species %in% passed_species_chlc) %>%
  group_by(species) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$Temp_Cat)) < 2) {
      return(data.frame(
        term = "Temp_Cat", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chlc per area_ug/cm2` ~ Temp_Cat, data = .x)
      tidy_result <- tidy(anova_result)
      return(tidy_result)
    }, error = function(e) {
      return(data.frame(
        term = "Temp_Cat", 
        df = NA, 
        sumsq = NA, 
        meansq = NA, 
        statistic = NA, 
        p.value = NA,
        error = as.character(e)
      ))
    })
  }) %>%
  filter(term != "Residuals") %>%
  mutate(
    variable = "Chlorophyll_c",
    analysis_by = "Species", 
    test_type = "ANOVA",
    Temp_Cat = NA  # Add this column for consistency
  )

# Print results to check
print(anova_results_chla_by_temp)
print(anova_results_chlc_by_temp)
print(anova_results_chla_by_species)
print(anova_results_chlc_by_species)



# 1. COMBINE ALL ANOVA RESULTS
all_anova_results <- bind_rows(
  anova_results_chla_by_temp,
  anova_results_chlc_by_temp,
  anova_results_chla_by_species,
  anova_results_chlc_by_species
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

# 2. EXTRACT TUKEY HSD RESULTS FOR PARAMETRIC ANALYSES
# Tukey for Chla by Temperature
tukey_chla_by_temp <- x8 %>%
  filter(Temp_Cat %in% passed_temp_chla) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chla per area_ug/cm2` ~ species, data = .x)
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
    variable = "Chlorophyll_a",
    analysis_by = "Temperature",
    test_type = "Tukey_HSD"
  )

# Tukey for Chlc by Temperature
tukey_chlc_by_temp <- x6 %>%
  filter(Temp_Cat %in% passed_temp_chlc) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chlc per area_ug/cm2` ~ species, data = .x)
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
    variable = "Chlorophyll_c",
    analysis_by = "Temperature",
    test_type = "Tukey_HSD"
  )

# Tukey for Chla by Species
tukey_chla_by_species <- x8 %>%
  filter(species %in% passed_species_chla) %>%
  group_by(species) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chla per area_ug/cm2` ~ Temp_Cat, data = .x)
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
    variable = "Chlorophyll_a",
    analysis_by = "Species",
    test_type = "Tukey_HSD",
    Temp_Cat = NA
  )

# Tukey for Chlc by Species
tukey_chlc_by_species <- x6 %>%
  filter(species %in% passed_species_chlc) %>%
  group_by(species) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chlc per area_ug/cm2` ~ Temp_Cat, data = .x)
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
    variable = "Chlorophyll_c",
    analysis_by = "Species",
    test_type = "Tukey_HSD",
    Temp_Cat = NA
  )

# Combine all Tukey results
all_tukey_results <- bind_rows(
  tukey_chla_by_temp,
  tukey_chlc_by_temp,
  tukey_chla_by_species,
  tukey_chlc_by_species
) %>%
  mutate(
    # Format for publication
    diff = round(diff, 4),
    lwr = round(lwr, 4),
    upr = round(upr, 4),
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

# 3. EXTRACT KRUSKAL-WALLIS RESULTS FOR NON-PARAMETRIC ANALYSES
# Kruskal-Wallis for Chla by Temperature
kruskal_chla_by_temp <- x8 %>%
  filter(Temp_Cat %in% failed_temp_chla) %>%
  group_by(Temp_Cat) %>%
  summarise(
    statistic = kruskal.test(`chla per area_ug/cm2` ~ species, data = cur_data())$statistic,
    p.value = kruskal.test(`chla per area_ug/cm2` ~ species, data = cur_data())$p.value,
    parameter = kruskal.test(`chla per area_ug/cm2` ~ species, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_a",
    analysis_by = "Temperature",
    test_type = "Kruskal-Wallis"
  )

# Kruskal-Wallis for Chlc by Temperature
kruskal_chlc_by_temp <- x6 %>%
  filter(Temp_Cat %in% failed_temp_chlc) %>%
  group_by(Temp_Cat) %>%
  summarise(
    statistic = kruskal.test(`chlc per area_ug/cm2` ~ species, data = cur_data())$statistic,
    p.value = kruskal.test(`chlc per area_ug/cm2` ~ species, data = cur_data())$p.value,
    parameter = kruskal.test(`chlc per area_ug/cm2` ~ species, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_c",
    analysis_by = "Temperature",
    test_type = "Kruskal-Wallis"
  )

# Kruskal-Wallis for Chla by Species
kruskal_chla_by_species <- x8 %>%
  filter(species %in% failed_species_chla) %>%
  group_by(species) %>%
  summarise(
    statistic = kruskal.test(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data())$statistic,
    p.value = kruskal.test(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data())$p.value,
    parameter = kruskal.test(`chla per area_ug/cm2` ~ Temp_Cat, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_a",
    analysis_by = "Species",
    test_type = "Kruskal-Wallis",
    Temp_Cat = NA
  )

# Kruskal-Wallis for Chlc by Species
kruskal_chlc_by_species <- x6 %>%
  filter(species %in% failed_species_chlc) %>%
  group_by(species) %>%
  summarise(
    statistic = kruskal.test(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data())$statistic,
    p.value = kruskal.test(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data())$p.value,
    parameter = kruskal.test(`chlc per area_ug/cm2` ~ Temp_Cat, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_c",
    analysis_by = "Species",
    test_type = "Kruskal-Wallis",
    Temp_Cat = NA
  )

# Combine all Kruskal-Wallis results
all_kruskal_results <- bind_rows(
  kruskal_chla_by_temp,
  kruskal_chlc_by_temp,
  kruskal_chla_by_species,
  kruskal_chlc_by_species
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
# Dunn's test for Chla by Temperature
dunn_chla_by_temp <- data.frame()

for(temp in failed_temp_chla) {
  temp_data <- x8 %>% filter(Temp_Cat == temp)
  
  tryCatch({
    dunn_result <- dunnTest(`chla per area_ug/cm2` ~ species, data = temp_data, method = "bonferroni")
    
    # Extract results and add metadata
    temp_results <- dunn_result$res %>%
      mutate(
        Temp_Cat = temp,
        variable = "Chlorophyll_a",
        analysis_by = "Temperature",
        test_type = "Dunn_Test"
      )
    
    dunn_chla_by_temp <- bind_rows(dunn_chla_by_temp, temp_results)
    
  }, error = function(e) {
    cat("Error for Chla temperature", temp, ":", e$message, "\n")
    
    error_row <- data.frame(
      Comparison = "Error",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      Temp_Cat = temp,
      variable = "Chlorophyll_a",
      analysis_by = "Temperature",
      test_type = "Dunn_Test",
      error = as.character(e$message)
    )
    
    dunn_chla_by_temp <- bind_rows(dunn_chla_by_temp, error_row)
  })
}

# Dunn's test for Chlc by Temperature
dunn_chlc_by_temp <- data.frame()

for(temp in failed_temp_chlc) {
  temp_data <- x6 %>% filter(Temp_Cat == temp)
  
  tryCatch({
    dunn_result <- dunnTest(`chlc per area_ug/cm2` ~ species, data = temp_data, method = "bonferroni")
    
    temp_results <- dunn_result$res %>%
      mutate(
        Temp_Cat = temp,
        variable = "Chlorophyll_c",
        analysis_by = "Temperature",
        test_type = "Dunn_Test"
      )
    
    dunn_chlc_by_temp <- bind_rows(dunn_chlc_by_temp, temp_results)
    
  }, error = function(e) {
    cat("Error for Chlc temperature", temp, ":", e$message, "\n")
    
    error_row <- data.frame(
      Comparison = "Error",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      Temp_Cat = temp,
      variable = "Chlorophyll_c",
      analysis_by = "Temperature",
      test_type = "Dunn_Test",
      error = as.character(e$message)
    )
    
    dunn_chlc_by_temp <- bind_rows(dunn_chlc_by_temp, error_row)
  })
}

# Dunn's test for Chla by Species
dunn_chla_by_species <- data.frame()

for(sp in failed_species_chla) {
  species_data <- x8 %>% filter(species == sp)
  
  tryCatch({
    dunn_result <- dunnTest(`chla per area_ug/cm2` ~ Temp_Cat, data = species_data, method = "bonferroni")
    
    species_results <- dunn_result$res %>%
      mutate(
        species = sp,
        variable = "Chlorophyll_a",
        analysis_by = "Species",
        test_type = "Dunn_Test",
        Temp_Cat = NA
      )
    
    dunn_chla_by_species <- bind_rows(dunn_chla_by_species, species_results)
    
  }, error = function(e) {
    cat("Error for Chla species", sp, ":", e$message, "\n")
    
    error_row <- data.frame(
      Comparison = "Error",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      species = sp,
      variable = "Chlorophyll_a",
      analysis_by = "Species",
      test_type = "Dunn_Test",
      Temp_Cat = NA,
      error = as.character(e$message)
    )
    
    dunn_chla_by_species <- bind_rows(dunn_chla_by_species, error_row)
  })
}

# Dunn's test for Chlc by Species
dunn_chlc_by_species <- data.frame()

for(sp in failed_species_chlc) {
  species_data <- x6 %>% filter(species == sp)
  
  tryCatch({
    dunn_result <- dunnTest(`chlc per area_ug/cm2` ~ Temp_Cat, data = species_data, method = "bonferroni")
    
    species_results <- dunn_result$res %>%
      mutate(
        species = sp,
        variable = "Chlorophyll_c",
        analysis_by = "Species",
        test_type = "Dunn_Test",
        Temp_Cat = NA
      )
    
    dunn_chlc_by_species <- bind_rows(dunn_chlc_by_species, species_results)
    
  }, error = function(e) {
    cat("Error for Chlc species", sp, ":", e$message, "\n")
    
    error_row <- data.frame(
      Comparison = "Error",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      species = sp,
      variable = "Chlorophyll_c",
      analysis_by = "Species",
      test_type = "Dunn_Test",
      Temp_Cat = NA,
      error = as.character(e$message)
    )
    
    dunn_chlc_by_species <- bind_rows(dunn_chlc_by_species, error_row)
  })
}

# Check the results
print(dunn_chla_by_temp)
print(dunn_chlc_by_temp)
print(dunn_chla_by_species)
print(dunn_chlc_by_species)

# Combine all Dunn's test results
all_dunn_results <- bind_rows(
  dunn_chla_by_temp,
  dunn_chlc_by_temp,
  dunn_chla_by_species,
  dunn_chlc_by_species
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
print(all_dunn_results)

# 5. CREATE SUMMARY STATISTICS
chlorophyll_stats_tables <- list(
  "Summary_Statistics" = combined_summary_stats,
  "ANOVA_Results" = all_anova_results,
  "Tukey_HSD_Results" = all_tukey_results,
  "Kruskal_Wallis_Results" = all_kruskal_results,
  "Dunn_PostHoc_Results" = all_dunn_results,  
  "Statistical_Methods_Used" = statistical_methods_summary,
  "Levene_Test_Results" = bind_rows(
    levene_results_by_temp_chla %>% mutate(Variable = "Chlorophyll_a", Analysis = "by_Temperature"),
    levene_results_by_temp_chlc %>% mutate(Variable = "Chlorophyll_c", Analysis = "by_Temperature")
  ),
  "Raw_Data_Chla" = x8 %>% 
    select(colony_id, species, Temp_Cat, `chla per area_ug/cm2`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "Temperature" = Temp_Cat,
      "Chlorophyll_a_ug_cm2" = `chla per area_ug/cm2`
    ),
  "Raw_Data_Chlc" = x6 %>% 
    select(colony_id, species, Temp_Cat, `chlc per area_ug/cm2`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "Temperature" = Temp_Cat,
      "Chlorophyll_c_ug_cm2" = `chlc per area_ug/cm2`
    )
)

# Save to Excel
write_xlsx(chlorophyll_stats_tables, "Chlorophyll_Concentration_Statistical_Analysis.xlsx")
