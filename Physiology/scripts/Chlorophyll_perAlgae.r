
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
  mutate(`chla (ug)/cell` = `Tot chla ug` / `surface.area.cm2`)

## Add a new column "Chlc per area" to Combined_Data - normalize by skeleton surface area
Combined_Data <- Combined_Data %>%
  mutate(`chlc (ug)/cell` = `Tot chlc ug` / `surface.area.cm2`)

Combined_Data <- Combined_Data %>%
  left_join(Data_algae, by = "colony_id") %>%
  distinct(colony_id, .keep_all = TRUE)

## Add a new column "Chla per algae cell" to Combined_Data - normalize algae cell
Combined_Data <- Combined_Data %>%
  mutate(`chla (ug)/cell` = `chla (ug)/cell` / `Cells/cm2`)
## Add a new column "Chlc per algae cell" to Combined_Data - normalize algae cell
Combined_Data <- Combined_Data %>%
  mutate(`chlc (ug)/cell` = `chlc (ug)/cell` / `Cells/cm2`)


# Plot results by temp - Chla
Combined_Data %>% 
  ggplot(aes(x = species, y = `chla (ug)/cell`, color = species, fill =species)) + 
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
  ggplot(aes(x = species, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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
  ggplot(aes(x = `Temp Cat`, y = `chla (ug)/cell`, color = species, fill =species)) + 
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
  ggplot(aes(x = `Temp Cat`, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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
  summarise(quant25 = quantile(`chla (ug)/cell`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chla (ug)/cell`, probs = q[2]),
            IQRbyGroup=IQR(`chla (ug)/cell`))

Quants_chlc <- Combined_Data %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chlc (ug)/cell`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chlc (ug)/cell`, probs = q[2]),
            IQRbyGroup=IQR(`chlc (ug)/cell`))


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
  filter(`chla (ug)/cell` < upper_chla & `chla (ug)/cell` > lower_chla)

x2 <- Combined_Data_with_bounds %>%
  filter(`chlc (ug)/cell` < upper_chlc & `chlc (ug)/cell` > lower_chlc)

# Plot results by temp - chla
x1 %>% 
  ggplot(aes(x = species, y = `chla (ug)/cell`, color = species, fill =species)) + 
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
  ggplot(aes(x = species, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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
  ggplot(aes(x = `Temp Cat`, y = `chla (ug)/cell`, color = species, fill =species)) + 
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
  ggplot(aes(x = `Temp Cat`, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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
x3 <- x1 %>% select(1, 4, 6, 13)

x4 <- x2 %>% select(1, 4, 6, 14)


Quants_chla <- x3 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chla (ug)/cell`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chla (ug)/cell`, probs = q[2]),
            IQRbyGroup=IQR(`chla (ug)/cell`))

Quants_chlc <- x4 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chlc (ug)/cell`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chlc (ug)/cell`, probs = q[2]),
            IQRbyGroup=IQR(`chlc (ug)/cell`))


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
  filter(`chla (ug)/cell` < upper_chla & `chla (ug)/cell` > lower_chla)

x10 <- Chlc_with_bounds %>%
  filter(`chlc (ug)/cell` < upper_chlc & `chlc (ug)/cell` > lower_chlc)

# Plot results - chla
x5 %>% 
  ggplot(aes(x = species, y = `chla (ug)/cell`, color = species, fill =species)) + 
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
plot_x10 <-  x10 %>% 
  ggplot(aes(x = species, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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

x10$`Temp Cat`= as.factor(x10$`Temp Cat`)

plot_x10 <-  x10 %>% 
  ggplot(aes(x = `Temp Cat`, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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


plot_x10


#### Remove outliers round 3

# Create a new table x3 with selected columns from x1
x7 <- x5 %>% select(1, 2,3,4)
x9 <- x10 %>% select(1, 2,3,4)

Quants_chla <- x7 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chla (ug)/cell`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chla (ug)/cell`, probs = q[2]),
            IQRbyGroup=IQR(`chla (ug)/cell`))

Quants_chlc <- x9 %>%
  group_by(species, `Temp Cat`) %>%
  summarise(quant25 = quantile(`chlc (ug)/cell`, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(`chlc (ug)/cell`, probs = q[2]),
            IQRbyGroup=IQR(`chlc (ug)/cell`))

#Calculate Quantile upper and lower ranges
Quants_chla$upper_chla <-  Quants_chla$quant75+1.5*Quants_chla$IQRbyGroup # Upper Range
Quants_chla$lower_chla <- Quants_chla$quant25-1.5*Quants_chla$IQRbyGroup # Lower Range

Quants_chlc$upper_chlc <-  Quants_chlc$quant75+1.5*Quants_chlc$IQRbyGroup # Upper Range
Quants_chlc$lower_chlc <- Quants_chlc$quant25-1.5*Quants_chlc$IQRbyGroup # Lower Range


# Join quantile ranges back to the original dataset
Chla_with_bounds <- x7 %>%
  left_join(Quants_chla, by = c("species", "Temp Cat"))

Chlc_with_bounds <- x9 %>%
  left_join(Quants_chlc, by = c("species", "Temp Cat"))

# Filter out outliers based on the quantile ranges
x9 <- Chla_with_bounds %>%
  filter(`chla (ug)/cell` < upper_chla & `chla (ug)/cell` > lower_chla)

x10 <- Chlc_with_bounds %>%
  filter(`chlc (ug)/cell` < upper_chlc & `chlc (ug)/cell` > lower_chlc)


# Plot results - chla
plot_x9 <- x9 %>% 
  ggplot(aes(x = species, y = `chla (ug)/cell`, color = species, fill =species)) + 
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

plot_x9


#chlc
plot_x10 <- x10 %>% 
  ggplot(aes(x = species, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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

plot_x10




## by species

x9$`Temp Cat`= as.factor(x9$`Temp Cat`)

plot_x9 <-  x9 %>% 
  ggplot(aes(x = `Temp Cat`, y = `chla (ug)/cell`, color = species, fill =species)) + 
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


plot_x9


#chlc
x10$`Temp Cat`= as.factor(x10$`Temp Cat`)

plot_x10 <-  x10 %>% 
  ggplot(aes(x = `Temp Cat`, y = `chlc (ug)/cell`, color = species, fill =species)) + 
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


plot_x10



# # Load the patchwork package
# library(patchwork)
# 
# # Combine the plots into a single figure
# combined_plot <- plot_x9 / plot_x10
# 
# # Display the combined plot
# combined_plot


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
x9 <- x9 %>%
  rename(Temp_Cat = `Temp Cat`)

x10 <- x10 %>%
  rename(Temp_Cat = `Temp Cat`)

  # Ensure Temp.Cat is a factor
x9$Temp_Cat <- as.factor(x9$Temp_Cat)
x10$Temp_Cat <- as.factor(x10$Temp_Cat)
x9$species <- as.factor(x9$species)
x10$species <- as.factor(x10$species)

  # Perform Shapiro-Wilk test for normality for each species
  shapiro_results_chla <- x9 %>% #chla
    group_by(species, Temp_Cat) %>%
    summarise(
      p_value = shapiro.test(`chla (ug)/cell`)$p.value
    )
  
  # View the results
  print(shapiro_results_chla) #all normally distributed

  shapiro_results_chlc <- x10 %>% #chlc
    group_by(species, Temp_Cat) %>%
    summarise(
      p_value = shapiro.test(`chlc (ug)/cell`)$p.value
    )
  
  # View the results
  print(shapiro_results_chlc) #all normally distributed
  
    # Generate Q-Q plots for each species
  qqplots <- x9 %>%
    group_by(species, Temp_Cat) %>%
    nest() %>%
    mutate(
      qqplot = map(data, ~ ggqqplot(.x$`chla (ug)/cell`, 
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
  levene_results_chla <- x9 %>%
    group_by(Temp_Cat) %>%  # Group by Temp_Cat
    summarise(
      levene_p_value = leveneTest(`chla (ug)/cell` ~ species, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  
    # View the results
  print(levene_results_chla) #18, 26.8 pass the test

  
  levene_results_chlc <- x10 %>%
    group_by(Temp_Cat) %>%  # Group by Temp_Cat
    summarise(
      levene_p_value = leveneTest(`chlc (ug)/cell` ~ species, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  # Set the number format to display decimals instead of scientific notation
  options(scipen = 999, digits = 4)
    # View the results
  print(levene_results_chlc) #26.8, 30, 18 pass the test
  
  
  
  levene_results_chla <- x9 %>%
    group_by(species) %>%  # Group by species
    summarise(
      levene_p_value = leveneTest(`chla (ug)/cell` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
  
  # View the results
  print(levene_results_chla)
  # # A tibble: 3 × 2
  # species            levene_p_value
  # <fct>                       <dbl>
  #   1 Montipora capitata        0.0405 
  # 2 Pocillopora acuta         0.00279
  # 3 Porites compressa         0.0510 
  
  
  levene_results_chlc <- x10 %>%
    group_by(species) %>%  # Group by species
    summarise(
      levene_p_value = leveneTest(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
      .groups = "drop"
    )
 
  # View the results
  print(levene_results_chlc) #12, 26.8 and 35 do not pass the test
  # # A tibble: 3 × 2
  # species            levene_p_value
  # <fct>                       <dbl>
  #   1 Montipora capitata       0.763   
  # 2 Pocillopora acuta        0.000313
  # 3 Porites compressa        0.0999  
  
  
  # Split Temp_Cat into those that passed and failed the Levene test - temp_cat
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
  
  # Split species into those that passed and failed the Levene test - species
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
  anova_results_chla <- x9 %>%
    filter(Temp_Cat %in% passed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      anova_p_value = summary(aov(`chla (ug)/cell` ~ species, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chla)
  # A tibble: 2 × 2
  # Temp_Cat anova_p_value
  # <fct>            <dbl>
  #   1 18             0.00309
  # 2 26.8           0.00833
  
  
  anova_results_chlc <- x10 %>%
    filter(Temp_Cat %in% passed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      anova_p_value = summary(aov(`chlc (ug)/cell` ~ species, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chlc)
  # # A tibble: 3 × 2
  # Temp_Cat anova_p_value
  # <fct>            <dbl>
  #   1 18              0.0143
  # 2 26.8            0.0147
  # 3 30              0.164 
  
  
  
  ### Perform ANOVA for species that passed the Levene test
  anova_results_chla <- x9 %>%
    filter(species %in% passed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      anova_p_value = summary(aov(`chla (ug)/cell` ~ Temp_Cat, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chla)
  # # A tibble: 1 × 2
  # species           anova_p_value
  # <fct>                     <dbl>
  #   1 Porites compressa        0.0221
  
  
  anova_results_chlc <- x10 %>%
    filter(species %in% passed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      anova_p_value = summary(aov(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data()))[[1]][["Pr(>F)"]][1],
      .groups = "drop"
    )
  
  # View ANOVA results
  print(anova_results_chlc)
  # # A tibble: 2 × 2
  # species            anova_p_value
  # <fct>                      <dbl>
  #   1 Montipora capitata         0.541
  # 2 Porites compressa          0.193
  
  
  
  # Perform Tukey's HSD test for Temp_Cat that passed the Levene test
  tukey_results_chla <- x9 %>%
    filter(Temp_Cat %in% passed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chla (ug)/cell` ~ species, data = cur_data()))),
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
  # A tibble: 6 × 6
  # Temp_Cat         diff         lwr          upr `p adj` Comparison                          
  # <fct>           <dbl>       <dbl>        <dbl>   <dbl> <chr>                               
  #   1 18       -0.00000439  -0.00000752 -0.00000127  0.00760 Pocillopora acuta-Montipora capitata
  # 2 18       -0.00000484  -0.00000787 -0.00000182  0.00316 Porites compressa-Montipora capitata
  # 3 18       -0.000000449 -0.00000304  0.00000214  0.887   Porites compressa-Pocillopora acuta 
  # 4 26.8     -0.00000711  -0.0000129  -0.00000131  0.0153  Pocillopora acuta-Montipora capitata
  # 5 26.8     -0.00000639  -0.0000118  -0.000000991 0.0192  Porites compressa-Montipora capitata
  # 6 26.8      0.000000722 -0.00000491  0.00000636  0.943   Porites compressa-Pocillopora acuta 
  
  
  tukey_results_chlc <- x10 %>%
    filter(Temp_Cat %in% passed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chlc (ug)/cell` ~ species, data = cur_data()))),
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
  # A tibble: 9 × 6
  # Temp_Cat     diff         lwr          upr `p adj` Comparison                          
  # <fct>       <dbl>       <dbl>        <dbl>   <dbl> <chr>                               
  #   1 18       -2.10e-6 -0.00000445  0.000000258  0.0836 Pocillopora acuta-Montipora capitata
  # 2 18       -2.97e-6 -0.00000533 -0.000000617  0.0139 Porites compressa-Montipora capitata
  # 3 18       -8.75e-7 -0.00000334  0.00000159   0.626  Porites compressa-Pocillopora acuta 
  # 4 26.8     -1.66e-6 -0.00000321 -0.000000108  0.0350 Pocillopora acuta-Montipora capitata
  # 5 26.8     -1.66e-6 -0.00000310 -0.000000218  0.0228 Porites compressa-Montipora capitata
  # 6 26.8     -1.78e-9 -0.00000151  0.00000150   1.00   Porites compressa-Pocillopora acuta 
  # 7 30       -1.99e-7 -0.00000160  0.00000120   0.931  Pocillopora acuta-Montipora capitata
  # 8 30       -9.98e-7 -0.00000236  0.000000360  0.176  Porites compressa-Montipora capitata
  # 9 30       -7.99e-7 -0.00000216  0.000000559  0.316  Porites compressa-Pocillopora acuta 
  
  
  
  # Perform Tukey's HSD test for species that passed the Levene test
  tukey_results_chla <- x9 %>%
    filter(species%in% passed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chla (ug)/cell` ~ Temp_Cat, data = cur_data()))),
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
  
  
  tukey_results_chlc <- x10 %>%
    filter(species %in% passed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      tukey = list(TukeyHSD(aov(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data()))),
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
  kruskal_results_chla <- x9 %>%
    filter(Temp_Cat %in% failed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chla (ug)/cell` ~ species, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chla)
  # Temp_Cat kruskal_p_value
  # <fct>              <dbl>
  #   1 12               0.0129 
  # 2 25               0.00516
  # 3 30               0.0361 
  # 4 35               0.00313
  
  kruskal_results_chlc <- x10 %>%
    filter(Temp_Cat %in% failed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chlc (ug)/cell` ~ species, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chlc)
  # Temp_Cat kruskal_p_value
  # <fct>              <dbl>
  #   1 12               0.0362 
  # 2 25               0.00835
  # 3 35               0.00539
  
  
  
  ###by species
  kruskal_results_chla <- x9 %>%
    filter(species %in% failed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chla (ug)/cell` ~ Temp_Cat, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chla)
  # # A tibble: 2 × 2
  # species            kruskal_p_value
  # <fct>                        <dbl>
  #   1 Montipora capitata           0.221
  # 2 Pocillopora acuta            0.611
  
  
  kruskal_results_chlc <- x10 %>%
    filter(species %in% failed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      kruskal_p_value = kruskal.test(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data())$p.value,
      .groups = "drop"
    )
  
  # View Kruskal-Wallis results
  print(kruskal_results_chlc)
  # # A tibble: 1 × 2
  # species           kruskal_p_value
  # <fct>                       <dbl>
  #   1 Pocillopora acuta           0.945
  
  
  
  
  # Load the FSA package
  library(FSA)
  
  
#### Perform Dunn's test for Temp_Cat that failed the Levene test
  dunn_results_chla <- x9 %>%
    filter(Temp_Cat %in% failed_temp_cat_chla) %>%
    group_by(Temp_Cat) %>%
    summarise(
      dunn_test = list(dunnTest(`chla (ug)/cell` ~ species, data = cur_data(), method = "bonferroni")),
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
  # # A tibble: 12 × 6
  # Temp_Cat Comparison                                 Z  P.unadj   P.adj comparison
  # <fct>    <chr>                                  <dbl>    <dbl>   <dbl> <chr>     
  #   1 12       Montipora capitata - Pocillopora acuta  1.58 0.114    0.343   1         
  # 2 12       Montipora capitata - Porites compressa  2.93 0.00334  0.0100  2         
  # 3 12       Pocillopora acuta - Porites compressa   1.38 0.168    0.503   3         
  # 4 25       Montipora capitata - Pocillopora acuta  1.24 0.216    0.648   1         
  # 5 25       Montipora capitata - Porites compressa  3.22 0.00129  0.00388 2         
  # 6 25       Pocillopora acuta - Porites compressa   1.98 0.0477   0.143   3         
  # 7 30       Montipora capitata - Pocillopora acuta -1.13 0.258    0.775   1         
  # 8 30       Montipora capitata - Porites compressa  1.45 0.146    0.438   2         
  # 9 30       Pocillopora acuta - Porites compressa   2.58 0.00995  0.0299  3         
  # 10 35       Montipora capitata - Pocillopora acuta  1.80 0.0714   0.214   1         
  # 11 35       Montipora capitata - Porites compressa  3.39 0.000689 0.00207 2         
  # 12 35       Pocillopora acuta - Porites compressa   1.59 0.112    0.335   3
  
  
 
  dunn_results_chlc <- x10 %>%
    filter(Temp_Cat %in% failed_temp_cat_chlc) %>%
    group_by(Temp_Cat) %>%
    summarise(
      dunn_test = list(dunnTest(`chlc (ug)/cell` ~ species, data = cur_data(), method = "bonferroni")),
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
  # # A tibble: 9 × 6
  # Temp_Cat Comparison                                 Z P.unadj   P.adj comparison
  # <fct>    <chr>                                  <dbl>   <dbl>   <dbl> <chr>     
  #   1 12       Montipora capitata - Pocillopora acuta 1.84  0.0663  0.199   1         
  # 2 12       Montipora capitata - Porites compressa 2.47  0.0135  0.0406  2         
  # 3 12       Pocillopora acuta - Porites compressa  0.549 0.583   1       3         
  # 4 25       Montipora capitata - Pocillopora acuta 0.565 0.572   1       1         
  # 5 25       Montipora capitata - Porites compressa 2.98  0.00287 0.00860 2         
  # 6 25       Pocillopora acuta - Porites compressa  2.05  0.0404  0.121   3         
  # 7 35       Montipora capitata - Pocillopora acuta 1.34  0.179   0.537   1         
  # 8 35       Montipora capitata - Porites compressa 3.22  0.00129 0.00388 2         
  # 9 35       Pocillopora acuta - Porites compressa  1.87  0.0610  0.183   3      
  
  
  ### by species
  dunn_results_chla <- x9 %>%
    filter(species%in% failed_temp_cat_chla) %>%
    group_by(species) %>%
    summarise(
      dunn_test = list(dunnTest(`chla (ug)/cell` ~ Temp_Cat , data = cur_data(), method = "bonferroni")),
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
  
  
  dunn_results_chlc <- x10 %>%
    filter(species %in% failed_temp_cat_chlc) %>%
    group_by(species) %>%
    summarise(
      dunn_test = list(dunnTest(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data(), method = "bonferroni")),
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
  
  
  

### Combine Dunn and Tukey results into a single data frame
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
      select(Temp_Cat, group1, group2, label),
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
  
  

  
  ###################### Add y.position column to significant_results_chla
  significant_results_chla <- significant_results_chla %>%
    left_join(
      x9 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group1" = "species")
    ) %>%
    rename(max_y_group1 = max_y) %>%
    left_join(
      x9 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group2" = "species")
    ) %>%
    rename(max_y_group2 = max_y) %>%
    mutate(
      y.position = pmax(max_y_group1, max_y_group2) + 0.000001  # Add a small padding above the highest data point
    )
  
  # Update y.position for significant_results_chla
  significant_results_chla <- significant_results_chla %>%
    group_by(Temp_Cat) %>%
    mutate(
      y.position = y.position + row_number() * 0.00001  # Increment height for each comparison
    ) %>%
    ungroup()
  
  # Update y.position for significant_results_chla
  significant_results_chla <- significant_results_chla %>%
    mutate(
      y.position = pmax(max_y_group1, max_y_group2) + 0.00000009  # Reduce padding above data points
    ) %>%
    group_by(Temp_Cat) %>%
    mutate(
      y.position = y.position + row_number() * 0.0000009  # Adjust increment for comparisons
    ) %>%
    ungroup()
  
  
  # Ensure group1 and group2 are factors
  significant_results_chla$group1 <- as.factor(significant_results_chla$group1)
  significant_results_chla$group2 <- as.factor(significant_results_chla$group2)
  
  
  #chlc
  significant_results_chlc <- significant_results_chlc %>%
    left_join(
      x10 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group1" = "species")
    ) %>%
    rename(max_y_group1 = max_y) %>%
    left_join(
      x10 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("Temp_Cat", "group2" = "species")
    ) %>%
    rename(max_y_group2 = max_y) %>%
    mutate(
      y.position = pmax(max_y_group1, max_y_group2) + 0.000001   # Add a small padding above the highest data point
    )
  
  # Update y.position for significant_results_chlc
  significant_results_chlc <- significant_results_chlc %>%
    group_by(Temp_Cat) %>%
    mutate(
      y.position = y.position + row_number() * 0.00001  # Increment height for each comparison
    ) %>%
    ungroup()
  
  # Update y.position for significant_results_chlc
  significant_results_chlc <- significant_results_chlc %>%
    mutate(
      y.position = pmax(max_y_group1, max_y_group2) + 0.00000005  # Reduce padding above data points
    ) %>%
    group_by(Temp_Cat) %>%
    mutate(
      y.position = y.position + row_number() * 0.0000005 # Adjust increment for comparisons
    ) %>%
    ungroup()
  

  # Ensure group1 and group2 are factors
  significant_results_chlc$group1 <- as.factor(significant_results_chlc$group1)
  significant_results_chlc$group2 <- as.factor(significant_results_chlc$group2)
  
  
  
  #### by species
  significant_results_chla <- significant_results_chla %>%
    left_join(
      x9 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group1" = "Temp_Cat")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x9 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chla (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group2" = "Temp_Cat")
    ) %>%
    rename(max_y_group2 = max_y) %>%
    mutate(
      y.position = pmax(max_y_group1, max_y_group2) + 0.000001  # Add a small padding above the highest data point
    )
  
  # Ensure group1 and group2 are factors
  significant_results_chla$group1 <- as.factor(significant_results_chla$group1)
  significant_results_chla$group2 <- as.factor(significant_results_chla$group2)
  
  
  #chlc
  significant_results_chlc <- significant_results_chlc %>%
    left_join(
      x10 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group1" = "Temp_Cat")
    ) %>%
    rename(max_y_group1 = max_y) %>%  # Rename for clarity
    left_join(
      x10 %>%
        group_by(Temp_Cat, species) %>%
        summarise(max_y = max(`chlc (ug)/cell`, na.rm = TRUE), .groups = "drop"),  # Get max y-value for each species
      by = c("species", "group2" = "Temp_Cat")
    ) %>%
    rename(max_y_group2 = max_y) %>%
    mutate(
      y.position = pmax(max_y_group1, max_y_group2) + 0.000001  # Add a small padding above the highest data point
    )
  
  # Ensure group1 and group2 are factors
  significant_results_chlc$group1 <- as.factor(significant_results_chlc$group1)
  significant_results_chlc$group2 <- as.factor(significant_results_chlc$group2)
  
  
  
  
  
  
  # Ensure the ggpubr package is loaded
  library(ggpubr)
  
  # # Prepare the annotations for statistical significance
  # annotations <- significant_results_chla %>%
  #   mutate(
  #     group1 = as.character(group1),  # Ensure group1 is a character
  #     group2 = as.character(group2)   # Ensure group2 is a character
  #   )
  # 
  # 
  # annotations <- significant_results_chlc %>%
  #   mutate(
  #     group1 = as.character(group1),  # Ensure group1 is a character
  #     group2 = as.character(group2)   # Ensure group2 is a character
  #   )
  
  
  plot_x9_stat <- x9 %>%
    ggplot(aes(x = species, y = `chla (ug)/cell`, color = species, fill = species)) +
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
      y.position = "y.position",  # Use dynamically calculated y.position
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )
  
    plot_x9_stat 
  
  
  
 plot_x10_stat <- x10 %>%
    ggplot(aes(x = species, y = `chlc (ug)/cell`, color = species, fill = species)) +
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
      y.position = "y.position",  # Use dynamically calculated y.position
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )

  
  plot_x10_stat
  
  
  
  # Load the patchwork package
  library(patchwork)
  
  # Combine the plots into a single figure
  combined_plot_stat <- plot_x9_stat / plot_x10_stat
  
  # Display the combined plot
  combined_plot_stat
  
  
  
  ####### by species
  
plot_x9_stat <- x9 %>%
    ggplot(aes(x = Temp_Cat, y = `chla (ug)/cell`, color = species, fill = species)) +
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
      y.position = "y.position",  # Use dynamically calculated y.position  # Decrease space between bars and data points
      inherit.aes = FALSE,
      hide.ns = TRUE  # Hide non-significant comparisons
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Ensure y-axis starts at 0
      expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
    )
  
  # Display the updated plot
  plot_x9_stat
  
  
  
plot_x10_stat <- x10 %>%
    ggplot(aes(x = Temp_Cat, y = `chlc (ug)/cell`, color = species, fill = species)) +
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
    ) 
# +
#     # Add statistical annotations
#     stat_pvalue_manual(
#       significant_results_chlc,
#       label = "label",  # Use the label column for significance levels
#       xmin = "group1",  # Start of the comparison
#       xmax = "group2",  # End of the comparison
#       y.position = "y.position",  # Use dynamically calculated y.position   # Decrease space between bars and data points
#       inherit.aes = FALSE,
#       hide.ns = TRUE  # Hide non-significant comparisons
#     ) +
#     scale_y_continuous(
#       limits = c(0, NA),  # Ensure y-axis starts at 0
#       expand = expansion(mult = c(0, 0.1))  # Remove extra space below and slightly above
#     )

  # Display the updated plot
  plot_x10_stat
  
  
  # Load the patchwork package
  library(patchwork)
  
  # Combine the plots into a single figure
  combined_plot_stat <- plot_x9_stat / plot_x10_stat
  
  # Display the combined plot
  combined_plot_stat
  
  
  
  
######### Save statistical tests resutls for supplementary 
  
library(writexl)
library(broom)
library(dplyr)
library(purrr)
library(car) 
library(FSA)


# Separation of TEMPERATURE vs SPECIES analysis

# 1. TEMPERATURE ANALYSIS: Run Levene's test by temperature (comparing species within each temp)
levene_results_by_temp_chla <- x9 %>%
  group_by(Temp_Cat) %>%
  summarise(
    levene_p_value = leveneTest(`chla (ug)/cell` ~ species, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

levene_results_by_temp_chlc <- x10 %>%
  group_by(Temp_Cat) %>%
  summarise(
    levene_p_value = leveneTest(`chlc (ug)/cell` ~ species, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

# Extract temperature categories that passed/failed Levene's test
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

# 2. SPECIES ANALYSIS: Run Levene's test by species (comparing temperatures within each species)
levene_results_by_species_chla <- x9 %>%
  group_by(species) %>%
  summarise(
    levene_p_value = leveneTest(`chla (ug)/cell` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

levene_results_by_species_chlc <- x10 %>%
  group_by(species) %>%
  summarise(
    levene_p_value = leveneTest(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  )

# Extract species that passed/failed Levene's test
passed_species_chla <- levene_results_by_species_chla %>%
  filter(levene_p_value > 0.05) %>%
  pull(species)

failed_species_chla <- levene_results_by_species_chla %>%
  filter(levene_p_value <= 0.05) %>%
  pull(species)

passed_species_chlc <- levene_results_by_species_chlc %>%
  filter(levene_p_value > 0.05) %>%
  pull(species)

failed_species_chlc <- levene_results_by_species_chlc %>%
  filter(levene_p_value <= 0.05) %>%
  pull(species)



# 1. ANOVA by Temperature (comparing species within each temperature)
anova_results_chla_by_temp <- x9 %>%
  filter(Temp_Cat %in% passed_temp_chla) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$species)) < 2) {
      return(data.frame(
        term = "species", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chla (ug)/cell` ~ species, data = .x)
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
    variable = "Chlorophyll_a_per_algae",
    analysis_by = "Temperature",
    test_type = "ANOVA"
  )

# 2. ANOVA by Species (comparing temperatures within each species)
nova_results_chla_by_species <- x9 %>%
  filter(species %in% passed_species_chla) %>%
  group_by(species) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$Temp_Cat)) < 2) {
      return(data.frame(
        term = "Temp_Cat", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chla (ug)/cell` ~ Temp_Cat, data = .x)
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
    variable = "Chlorophyll_a_per_algae",
    analysis_by = "Species",
    test_type = "ANOVA",
    Temp_Cat = NA
  )



# Similar for Chlorophyll c per algae
anova_results_chlc_by_temp <- x10 %>%
  filter(Temp_Cat %in% passed_temp_chlc) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$species)) < 2) {
      return(data.frame(
        term = "species", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chlc (ug)/cell` ~ species, data = .x)
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
    variable = "Chlorophyll_c_per_algae",
    analysis_by = "Temperature",
    test_type = "ANOVA"
  )

anova_results_chlc_by_species <- x10 %>%
  filter(species %in% passed_species_chlc) %>%
  group_by(species) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || length(unique(.x$Temp_Cat)) < 2) {
      return(data.frame(
        term = "Temp_Cat", df = NA, sumsq = NA, meansq = NA, 
        statistic = NA, p.value = NA, error = "Insufficient data"
      ))
    }
    
    tryCatch({
      anova_result <- aov(`chlc (ug)/cell` ~ Temp_Cat, data = .x)
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
    variable = "Chlorophyll_c_per_algae",
    analysis_by = "Species", 
    test_type = "ANOVA",
    Temp_Cat = NA
  )



# COMBINE ALL ANOVA RESULTS
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


# EXTRACT TUKEY HSD RESULTS FOR PARAMETRIC ANALYSES
# Tukey for Chla by Temperature
tukey_chla_by_temp <- x9 %>%
  filter(Temp_Cat %in% passed_temp_chla) %>%
  group_by(Temp_Cat) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chla (ug)/cell` ~ species, data = .x)
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
    variable = "Chlorophyll_a_per_algae",
    analysis_by = "Temperature",
    test_type = "Tukey_HSD"
  )

# Tukey for Chlc by Temperature
tukey_chla_by_species <- x9 %>%
  filter(species %in% passed_species_chla) %>%
  group_by(species) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chla (ug)/cell` ~ Temp_Cat, data = .x)
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
    variable = "Chlorophyll_a_per_algae",  
    analysis_by = "Species",
    test_type = "Tukey_HSD",
    Temp_Cat = NA
  )

# Tukey for Chlc by Species
tukey_chlc_by_species <- x10 %>%
  filter(species %in% passed_species_chlc) %>%
  group_by(species) %>%
  group_modify(~ {
    tryCatch({
      anova_result <- aov(`chlc (ug)/cell` ~ Temp_Cat, data = .x)
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
    variable = "Chlorophyll_c_per_algae", 
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


# EXTRACT KRUSKAL-WALLIS RESULTS FOR NON-PARAMETRIC ANALYSES
# Kruskal-Wallis for Chla by Temperature
kruskal_chla_by_temp <- x9 %>%
  filter(Temp_Cat %in% failed_temp_chla) %>%
  group_by(Temp_Cat) %>%
  summarise(
    statistic = kruskal.test(`chla (ug)/cell` ~ species, data = cur_data())$statistic,
    p.value = kruskal.test(`chla (ug)/cell` ~ species, data = cur_data())$p.value,
    parameter = kruskal.test(`chla (ug)/cell` ~ species, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_a_per_algae", 
    analysis_by = "Temperature",
    test_type = "Kruskal-Wallis"
  )

# Kruskal-Wallis for Chlc by Temperature
kruskal_chlc_by_temp <- x10 %>%
  filter(Temp_Cat %in% failed_temp_chlc) %>%
  group_by(Temp_Cat) %>%
  summarise(
    statistic = kruskal.test(`chlc (ug)/cell` ~ species, data = cur_data())$statistic,
    p.value = kruskal.test(`chlc (ug)/cell` ~ species, data = cur_data())$p.value,
    parameter = kruskal.test(`chlc (ug)/cell` ~ species, data = cur_data())$parameter,
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_c_per_algae", 
    analysis_by = "Temperature",
    test_type = "Kruskal-Wallis"
  )

# Kruskal-Wallis for Chla by Species
kruskal_chla_by_species <- x9 %>%
  filter(species %in% failed_species_chla) %>%
  group_by(species) %>%
  summarise(
    n_temps = n_distinct(Temp_Cat),
    n_obs = n(),
    statistic = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                      kruskal.test(`chla (ug)/cell` ~ Temp_Cat, data = cur_data())$statistic, 
                      NA),
    p.value = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                    kruskal.test(`chla (ug)/cell` ~ Temp_Cat, data = cur_data())$p.value, 
                    NA),
    parameter = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                      kruskal.test(`chla (ug)/cell` ~ Temp_Cat, data = cur_data())$parameter, 
                      NA),
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_a_per_algae",
    analysis_by = "Species",
    test_type = "Kruskal-Wallis",
    Temp_Cat = NA,
    error = ifelse(is.na(statistic), "Insufficient groups or data", NA)
  ) %>%
  select(-n_temps, -n_obs)

# Kruskal-Wallis for Chlc by Species
kruskal_chlc_by_species <- x10 %>%
  filter(species %in% failed_species_chlc) %>%
  group_by(species) %>%
  summarise(
    n_temps = n_distinct(Temp_Cat),
    n_obs = n(),
    statistic = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                      kruskal.test(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data())$statistic, 
                      NA),
    p.value = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                    kruskal.test(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data())$p.value, 
                    NA),
    parameter = ifelse(n_distinct(Temp_Cat) > 1 & n() > 2, 
                      kruskal.test(`chlc (ug)/cell` ~ Temp_Cat, data = cur_data())$parameter, 
                      NA),
    .groups = "drop"
  ) %>%
  mutate(
    variable = "Chlorophyll_c_per_algae",
    analysis_by = "Species",
    test_type = "Kruskal-Wallis",
    Temp_Cat = NA,
    error = ifelse(is.na(statistic), "Insufficient groups or data", NA)
  ) %>%
  select(-n_temps, -n_obs)


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
  temp_data <- x9 %>% filter(Temp_Cat == temp)  
  
  tryCatch({
    dunn_result <- dunnTest(`chla (ug)/cell` ~ species, data = temp_data, method = "bonferroni")
    
    # Extract results and add metadata
    temp_results <- dunn_result$res %>%
      mutate(
        Temp_Cat = temp,
        variable = "Chlorophyll_a_per_algae",
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
      variable = "Chlorophyll_a_per_algae",
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
  temp_data <- x10 %>% filter(Temp_Cat == temp)
  
  tryCatch({
    dunn_result <- dunnTest(`chlc (ug)/cell` ~ species, data = temp_data, method = "bonferroni")
    
    temp_results <- dunn_result$res %>%
      mutate(
        Temp_Cat = temp,
        variable = "Chlorophyll_c_per_algae",
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
      variable = "Chlorophyll_c_per_algae",
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
  species_data <- x9 %>% filter(species == sp) 
  
  # Check if there are enough groups and observations
  if(n_distinct(species_data$Temp_Cat) > 1 && nrow(species_data) > 2) {
    tryCatch({
      dunn_result <- dunnTest(`chla (ug)/cell` ~ Temp_Cat, data = species_data, method = "bonferroni")
      
      species_results <- dunn_result$res %>%
        mutate(
          species = sp,
          variable = "Chlorophyll_a_per_algae",
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
        variable = "Chlorophyll_a_per_algae",
        analysis_by = "Species",
        test_type = "Dunn_Test",
        Temp_Cat = NA,
        error = as.character(e$message)
      )
      
      dunn_chla_by_species <- bind_rows(dunn_chla_by_species, error_row)
    })
  } else {
    # Not enough data for testing
    error_row <- data.frame(
      Comparison = "Insufficient data",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      species = sp,
      variable = "Chlorophyll_a_per_algae",
      analysis_by = "Species",
      test_type = "Dunn_Test",
      Temp_Cat = NA,
      error = paste("Only", n_distinct(species_data$Temp_Cat), "temperature groups with", nrow(species_data), "observations")
    )
    
    dunn_chla_by_species <- bind_rows(dunn_chla_by_species, error_row)
  }
}


# Dunn's test for Chlc by Species
dunn_chlc_by_species <- data.frame()

for(sp in failed_species_chlc) {
  species_data <- x10 %>% filter(species == sp) 
  
  # Check if there are enough groups and observations
  if(n_distinct(species_data$Temp_Cat) > 1 && nrow(species_data) > 2) {
    tryCatch({
      dunn_result <- dunnTest(`chlc (ug)/cell` ~ Temp_Cat, data = species_data, method = "bonferroni")
      
      species_results <- dunn_result$res %>%
        mutate(
          species = sp,
          variable = "Chlorophyll_c_per_algae",
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
        variable = "Chlorophyll_c_per_algae",
        analysis_by = "Species",
        test_type = "Dunn_Test",
        Temp_Cat = NA,
        error = as.character(e$message)
      )
      
      dunn_chlc_by_species <- bind_rows(dunn_chlc_by_species, error_row)
    })
  } else {
    # Not enough data for testing
    error_row <- data.frame(
      Comparison = "Insufficient data",
      Z = NA,
      P.unadj = NA,
      P.adj = NA,
      species = sp,
      variable = "Chlorophyll_c_per_algae",
      analysis_by = "Species",
      test_type = "Dunn_Test",
      Temp_Cat = NA,
      error = paste("Only", n_distinct(species_data$Temp_Cat), "temperature groups with", nrow(species_data), "observations")
    )
    
    dunn_chlc_by_species <- bind_rows(dunn_chlc_by_species, error_row)
  }
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
summary_stats_chla <- x9 %>%
  group_by(species, Temp_Cat) %>%
  summarise(
    N = n(),
    Mean = round(mean(`chla (ug)/cell`, na.rm = TRUE), 6),
    SD = round(sd(`chla (ug)/cell`, na.rm = TRUE), 6),
    SE = round(sd(`chla (ug)/cell`, na.rm = TRUE) / sqrt(n()), 6),
    Median = round(median(`chla (ug)/cell`, na.rm = TRUE), 6),
    IQR = round(IQR(`chla (ug)/cell`, na.rm = TRUE), 6),
    Min = round(min(`chla (ug)/cell`, na.rm = TRUE), 6),
    Max = round(max(`chla (ug)/cell`, na.rm = TRUE), 6),
    .groups = "drop"
  ) %>%
  mutate(Variable = "Chlorophyll_a_per_algae")

summary_stats_chlc <- x10 %>%
  group_by(species, Temp_Cat) %>%
  summarise(
    N = n(),
    Mean = round(mean(`chlc (ug)/cell`, na.rm = TRUE), 6),
    SD = round(sd(`chlc (ug)/cell`, na.rm = TRUE), 6),
    SE = round(sd(`chlc (ug)/cell`, na.rm = TRUE) / sqrt(n()), 6),
    Median = round(median(`chlc (ug)/cell`, na.rm = TRUE), 6),
    IQR = round(IQR(`chlc (ug)/cell`, na.rm = TRUE), 6),
    Min = round(min(`chlc (ug)/cell`, na.rm = TRUE), 6),
    Max = round(max(`chlc (ug)/cell`, na.rm = TRUE), 6),
    .groups = "drop"
  ) %>%
  mutate(Variable = "Chlorophyll_c_per_algae")

combined_summary_stats <- bind_rows(summary_stats_chla, summary_stats_chlc)

# Create methods summary
statistical_methods_summary <- data.frame(
  Variable = c("Chlorophyll_a_per_algae", "Chlorophyll_c_per_algae"),
  Units = c("μg chlorophyll a per algae cell", "μg chlorophyll c per algae cell"),
  Normality_Test = c("Shapiro-Wilk", "Shapiro-Wilk"),
  Homogeneity_Test = c("Levene's test", "Levene's test"),
  Parametric_Test = c("ANOVA + Tukey HSD", "ANOVA + Tukey HSD"),
  NonParametric_Test = c("Kruskal-Wallis + Dunn's test", "Kruskal-Wallis + Dunn's test"),
  Significance_Level = c("α = 0.05", "α = 0.05"),
  Multiple_Comparisons = c("Bonferroni correction", "Bonferroni correction")
)

# Create final Excel workbook
chlorophyll_per_algae_stats_tables <- list(
  "Summary_Statistics" = combined_summary_stats,
  "ANOVA_Results" = all_anova_results,
  "Tukey_HSD_Results" = all_tukey_results,
  "Kruskal_Wallis_Results" = all_kruskal_results,
  "Dunn_PostHoc_Results" = all_dunn_results,
  "Statistical_Methods_Used" = statistical_methods_summary,
  "Levene_Test_Results" = bind_rows(
    levene_results_by_temp_chla %>% mutate(Variable = "Chlorophyll_a_per_algae", Analysis = "by_Temperature"),
    levene_results_by_temp_chlc %>% mutate(Variable = "Chlorophyll_c_per_algae", Analysis = "by_Temperature"),
    levene_results_by_species_chla %>% mutate(Variable = "Chlorophyll_a_per_algae", Analysis = "by_Species"),
    levene_results_by_species_chlc %>% mutate(Variable = "Chlorophyll_c_per_algae", Analysis = "by_Species")
  ),
  "Raw_Data_Chla_per_algae" = x9 %>% 
    select(colony_id, species, Temp_Cat, `chla (ug)/cell`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "Temperature" = Temp_Cat,
      "Chlorophyll_a_ug_per_algae_cell" = `chla (ug)/cell`
    ),
  "Raw_Data_Chlc_per_algae" = x10 %>% 
    select(colony_id, species, Temp_Cat, `chlc (ug)/cell`) %>%
    rename(
      "Colony_ID" = colony_id,
      "Species" = species,
      "Temperature" = Temp_Cat,
      "Chlorophyll_c_ug_per_algae_cell" = `chlc (ug)/cell`
    )
)

# Create output directory and save
write_xlsx(chlorophyll_per_algae_stats_tables, "Chlorophyll_Per_Algae_Cell_Statistical_Analysis.xlsx")
