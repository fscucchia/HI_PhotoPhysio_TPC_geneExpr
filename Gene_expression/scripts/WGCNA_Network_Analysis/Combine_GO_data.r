
########### Combined GO enrichment for Mcap, Pacu and Pcom

library(dplyr)
library(ggplot2)
library(forcats)

# Load tables for each species
topGO_pcom <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Pcom/Top10_GO_p01Fisher_Pcom.csv")
topGO_mcap <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/Top10_GO_p01Fisher_Mcap.csv")
topGO_pacu <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Pacu/Top10_GO_p01Fisher_Pacu.csv")

slim_pcom <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Pcom/GOslim_All_p01Fisher_Pcom.csv")
slim_mcap <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Mcap/GOslim_All_p01Fisher_Mcap.csv")
slim_pacu <- read.csv("/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/Functional_Enrichment/Pacu/GOslim_All_p01Fisher_Pacu.csv")

# Add species column
topGO_pcom <- topGO_pcom %>% mutate(species = "Pcom")
topGO_mcap <- topGO_mcap %>% mutate(species = "Mcap")
topGO_pacu <- topGO_pacu %>% mutate(species = "Pacu")

slim_pcom <- slim_pcom %>% mutate(species = "Pcom")
slim_mcap <- slim_mcap %>% mutate(species = "Mcap")
slim_pacu <- slim_pacu %>% mutate(species = "Pacu")

# Combine all tables
topGO_all <- bind_rows(topGO_pcom, topGO_mcap, topGO_pacu)
slim_all <- bind_rows(slim_pcom, slim_mcap, slim_pacu)

# View combined table
head(topGO_all)
head(slim_all)

##Prepare data for plotting

# Make sure temp is a factor with desired order
topGO_all <- topGO_all %>% mutate(temp = factor(temperature, levels = c("control", "30", "35")))
slim_all  <- slim_all  %>% mutate(temp = factor(temperature, levels = c("control", "30", "35")))

## topGO heatmap

ggplot(topGO_all, aes(x = species, y = term, fill = direction)) +
  geom_tile(color = "white") +
  facet_grid(rows = vars(temp), switch = "y", scales = "free_y", space = "free") +
  scale_fill_manual(
    name = "Direction",
    values = c("up" = "#FDE725", "down" = "#440154")
  ) +
  labs(x = "Species", y = "GO Term", title = "GO Terms by Species and Temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  )


# Filter for shared terms within each temperature group
shared_terms <- topGO_all %>%
  group_by(term, temp) %>%
  filter(n_distinct(species) > 1) %>%
  ungroup()

# Plot with free y-axis for each facet
ggplot(shared_terms, aes(x = species, y = term, fill = direction)) +
  geom_tile(color = "white") +
  facet_grid(rows = vars(temp), switch = "y", scales = "free_y", space = "free") +
  scale_fill_manual(
    name = "Direction",
    values = c("up" = "#FDE725", "down" = "#440154")
  ) +
  labs(x = "Species", y = "GO Term", title = "Shared GO Terms by Species and Temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  )


### slimGO heatmap

ggplot(slim_all, aes(x = species, y = slim_term, fill = direction)) +
  geom_tile(color = "white") +
  facet_grid(rows = vars(temp), switch = "y", scales = "free_y", space = "free") +
  scale_fill_manual(
    name = "Direction",
    values = c("up" = "#FDE725", "down" = "#440154")
  ) +
  labs(x = "Species", y = "GO Term", title = "GO Terms by Species and Temperature") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  )


  ## Shared terms (present in ≥2 species per temperature)
  slim_shared <- slim_all %>%
  group_by(slim_term, temp) %>%
  filter(n_distinct(species) > 1) %>%
  ungroup()

ggplot(slim_shared, aes(x = species, y = slim_term, fill = direction)) +
  geom_tile(color = "white") +
  facet_grid(rows = vars(temp), switch = "y", scales = "free_y", space = "free") +
  scale_fill_manual(
    name = "Direction",
    values = c("up" = "#fd255fff", "down" = "#240154ff")
  ) +
  labs(x = "Species", y = "GO Term", title = "Shared GO Terms by Species and Temperature") +
theme_minimal() +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  strip.placement = "outside",
  strip.text.y.left = element_text(angle = 0),
  panel.grid = element_blank()
)

## Non-shared terms (present in only 1 species per temperature)
slim_nonshared <- slim_all %>%
  group_by(slim_term, temp) %>%
  filter(n_distinct(species) == 1) %>%
  ungroup()

ggplot(slim_nonshared, aes(x = species, y = slim_term, fill = direction)) +
  geom_tile(color = "white") +
  facet_grid(rows = vars(temp), switch = "y", scales = "free_y", space = "free") +
  scale_fill_manual(
    name = "Direction",
    values = c("up" = "#fd255fff", "down" = "#240154ff")
  ) +
  labs(x = "Species", y = "GO Term", title = "Unique GO Terms by Species and Temperature") +
theme_minimal() +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  strip.placement = "outside",
  strip.text.y.left = element_text(angle = 0),
  panel.grid = element_blank()
)