
################################# "Oxygen flux rate extractions - generate PI curves"

## install packages if you dont already have them in your library
if (!require("devtools")) install.packages("devtools")
if (!require("furrr")) install.packages("furrr")
if (!require("future")) install.packages("future")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("lubridate")) install.packages("lubridate")
if (!require("cowplot")) install.packages("cowplot")
if (!require("LoLinR")) install_github('colin-olito/LoLinR') 

## load libraries
library(devtools)
library(LoLinR)
library(tidyverse)
library(lubridate)
library(cowplot)
library(broom)
library(plotrix)
library(Hmisc)
library(rTPC)
library(nls.multstart)
library(here)
library(nlraa)
library(ggplot2)

## libraries for parallel processing
library(future)
library(furrr)


# Define data   

## specify data 
Data<-read_csv(here("output","pi_curve_extracted_rates_meanTemp_thinned20.csv"))


Data$PAR <- as.numeric(Data$Light_Value.x) 
Data$Pc <- as.numeric(Data$micromol.cm2.s) 
#Data$Pc <- as.numeric(Data$micromol.cm2.h) 

Data %>% 
ggplot(aes(x=PAR, y=Pc, color=Species))+ 
  geom_point() 


Data %>% 
  ggplot(aes(x=PAR, y=Pc, color=Species, group=colony_id))+ 
  geom_point()+
  geom_line()+
  facet_wrap(~ Species*Temp.Cat.x,  ncol = 6)+
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    "Montipora capitata" = "green",
    "Pocillopora acuta" = "cyan",
    "Porites compressa" = "orange"
  )) +
  labs(y = expression(paste(mu, "mol O"[2], " ", cm^{-2}, " ", h^{-1})))

Data %>% 
  ggplot(aes(x = PAR, y = Pc, color = Species, group = colony_id)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~ Temp.Cat.x, ncol = 6, labeller = label_value) +
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
  )) +
  labs(y = expression(paste(mu, "mol O"[2], " ", cm^{-2}, " ", h^{-1})))



### create PI curve plots per each temperature

# Create a directory to save the plots if it doesn't exist
output_dir <- "output/RAW_PIcurves_PlotPerTemp_lostTempsRemoved"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate and save one plot per Temp.Cat.x
unique_temp_cat <- unique(Data$Temp.Cat.x)

for (temp_cat in unique_temp_cat) {
  plot <- Data %>%
    filter(Temp.Cat.x == temp_cat) %>%
    ggplot(aes(x = PAR, y = Pc, color = Species, group = colony_id)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c(
      "Montipora capitata" = "green",
      "Pocillopora acuta" = "cyan",
      "Porites compressa" = "orange")) +
    labs(title = paste("Temp.Cat.x:", temp_cat)) +
    theme(legend.position = "right")
  
  ggsave(filename = paste0(output_dir, "/pi_curve_raw_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
}

ests <- Data%>%
  group_by(Species, Temp.Cat.x, colony_id)%>%
  summarise(Rd= min(Pc),
            Pmax = max(Pc))


dirty.tpc.plot.Rd <- ests %>% 
  ggplot(aes(x=Temp.Cat.x, y=-Rd, color=Species, group=Species))+ 
  geom_point()+
  #geom_line()+
  facet_wrap(~ Species)+
  theme(legend.position = "none")

dirty.tpc.plot.Rd  


dirty.tpc.plot.Pmax <- ests %>% 
  ggplot(aes(x=Temp.Cat.x, y=Pmax, color=Species, group=Species))+ 
  geom_point()+
  #geom_line()+
  facet_wrap(~ Species)+
  theme(legend.position = "none")

dirty.tpc.plot.Pmax 



# Define PI curve function as a nonlinear Least Squares regression of a quadratic fit, test nls fit 
#Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), data=., start=list(Am=(max(.$Pc)-min(.$Pc)),  AQY=0.001, Rd=-min(.$Pc)) 
#Aquatic Photosynthesis, Falkowski  
#PAR = irradiance from 400-700nm (also called I or E) 
#PC = oxygen flux rate 
#Pmax = max photosynthesis (also called Am)   
#alpha = quantum yeild (also called AQY)   
#I or E or PAR = irradiance  
#Rd = dark respiration 
#Ik (saturating irradiance) is the point at which photosynthesis reaches the max of initial slope = Am/AQY 
#Ic=(Am*Rd)/(AQY*(sqrt(Am^2-Rd^2))) 
#Equation for Ic derived from quadratic equation above. Ic = Par when Pc = 0 (x intercept). Ic = light compensation point; point at which photosynthesis is released from carbon limitation.  

##### Run nls model Using flexible initial values based on input data: 

############################################# Example
require(ggplot2)
library(nlraa)
set.seed(12345)

Data$PAR <- as.numeric(Data$Light_Value.x) 
Data$Pc <- as.numeric(Data$micromol.cm2.s) 

x <- Data$PAR
y <- Data$Pc
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSnrh(x, asym, phi, theta, rd), data = dat)

## Visualize observed and simulated
ggplot(data = dat, aes(x = x, y = y)) + 
  geom_point() + 
  geom_line(aes(y = fitted(fit)))
## Testing predict function
prd <- predict_nls(fit, interval = "confidence")
datA <- cbind(dat, prd)
## Plotting
ggplot(data = datA, aes(x = x, y = y)) + 
  geom_point() + 
  geom_line(aes(y = fitted(fit))) + 
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), 
  fill = "purple", alpha = 0.3)
#########################################

#Step 1: Define a Function for Colony-Specific Fitting
#Create a function to fit the nls model for a single colony and return the results.
#Am = asym, AQY = phi, Rd = rd

library(minpack.lm)
library(broom)

# Define the fitting function
fit_colony_unscaled <- function(data) {
  tryCatch({
    # Fit the model
    fit <- nlsLM(Pc ~ (Am * ((AQY * PAR) / sqrt(Am^2 + (AQY * PAR)^2)) - Rd),
                 data = data,
                 start = list(Am = 0.7, AQY = 0.001, Rd = 0.4), #based on the data
                 control = nls.lm.control(maxiter = 1000))
    
    # Extract model parameters and goodness-of-fit metrics
    params <- tidy(fit)
    metrics <- glance(fit)
    
    # Combine parameters and metrics into one table
    results <- bind_cols(params, metrics)
    results
  }, error = function(e) {
    # Return NA for failed fits
    tibble(term = NA, estimate = NA, std.error = NA, statistic = NA,
           p.value = NA, r.squared = NA, adj.r.squared = NA,
           sigma = NA, AIC = NA, BIC = NA)
  })
}

#Step 2: Apply the Function to Each Colony
#Use split and map to apply the function to each colony separately.

library(purrr)

# Split the data by colony_id
colony_data <- split(Data, Data$colony_id)

# Apply the fitting function to each colony
colony_results <- map(colony_data, fit_colony_unscaled)

#Step 3: Combine the Results
#Combine the results into a single table and include the colony_id.

# Add colony_id to each result
colony_results_with_id <- map2_dfr(colony_results, names(colony_results), ~ mutate(.x, colony_id = .y))

# View the combined results
colony_results_with_id

#Step 4: Save the Results
#Write the combined results to a CSV file.

write_csv(colony_results_with_id, "output/nls_fitted_params_per_colony_Falkowsky.csv")
#write_csv(colony_results_with_id, "output/nls_fitted_params_per_colony_minusFailed.csv")

# #Step 5: Visualize the Results (Optional)
# #Visualize the fitted parameters to identify patterns or diagnose issues.
# 
# library(ggplot2)
# 
# # Plot parameter estimates
# colony_results_with_id %>%
#   filter(!is.na(term)) %>%
#   ggplot(aes(x = colony_id, y = estimate, fill = term)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "Parameter Estimates by Colony", x = "Colony ID", y = "Estimate")

#Step 6: Extract the Fitted Parameters
#extract the parameter estimates (asym, phi, theta, and rd) for each colony and compute the derived metrics (Am, AQY, Ik, Ic).

# Extract and reshape parameter estimates
colony_params <- colony_results_with_id %>%
  filter(!is.na(term)) %>%         # Remove rows with failed fits
  select(colony_id, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) 

# Compute derived metrics (Ik and Ic)
colony_params <- colony_params %>%
  mutate(
    Ik = Am / AQY,  # Saturating irradiance
    Ic = (Am * Rd) / (AQY * sqrt(Am^2 - Rd^2))  # Compensation point
  )

write_csv(colony_params, "output/pi_curve_fitted_metrics_thinned20.csv")
write_csv(colony_params, "output/pi_curve_fitted_metrics_Falkowsky.csv")


# #### Identify colony IDs that failed the fitting
# failed_fits <- setdiff(unique(Data$colony_id), unique(colony_params$colony_id))
# 
# #### Output the failed fits in a datatable
# library(DT)
# datatable(data.frame(Failed_Colony_IDs = failed_fits))
# 
# ##### Write the failed fits to a CSV file
# write_csv(data.frame(Failed_Colony_IDs = failed_fits), "output/failed_fits_thinned20.csv")


#Step 7: Plot the fitted curve over the data points for each colony to groundtruth output

# Define a function to generate the fitted curve
generate_fitted_curve <- function(PAR, Am, AQY, theta, Rd) {
  (Am * ((AQY * PAR) / (sqrt(Am^2 + (AQY * PAR)^2))) - Rd)
}

# Plot the data points and fitted curve for each colony
plot_fitted_curve <- function(data, params) {
  ggplot(data, aes(x = PAR, y = Pc)) +
    geom_point() +
    geom_line(aes(y = generate_fitted_curve(PAR, params$Am, params$AQY, params$theta, params$Rd)), color = "blue") +
    labs(title = paste("Colony ID:", unique(data$colony_id)),
         x = "PAR",
         y = "Pc") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_line(color = "gray95"))
}

# Filter colony_data to include only colonies with successful fits
successful_colonies <- colony_params$colony_id
filtered_colony_data <- colony_data[successful_colonies]

# Apply the plotting function to each colony and save the plots
walk2(filtered_colony_data, split(colony_params, colony_params$colony_id), function(data, params) {
  plot <- plot_fitted_curve(data, params)
  ggsave(filename = paste0("C:/Users/Federica/OneDrive - University of Rhode Island/Postdoc_URI/TPC_Hawaii/R analysis/Hawaii/output/Fitted_IndividualPlots_nlsLM_updated/plot_colony_", unique(data$colony_id), ".png"), plot = plot)
})


####### Pacu-C1, Pacu-D6, Pacu-H6  look pretty weird!!!


# Step 8: Generate a plot with all curves for all colony_id values, colored by species
# Combine data and parameters for plotting
combined_data <- bind_rows(filtered_colony_data, .id = "colony_id")
combined_params <- left_join(combined_data, colony_params, by = "colony_id")

# Generate the fitted curves for all colonies
combined_params <- combined_params %>%
  #rowwise() %>%
  mutate(fitted_curve = generate_fitted_curve(PAR, Am, AQY, theta, Rd))
#mutate(fitted_curve = generate_fitted_curve(PAR, Pmax.gross, AQY, theta, Rd))

# Plot all curves colored by species
ggplot(combined_params, aes(x = PAR, y = fitted_curve, color = Species, group = colony_id)) +
  geom_line() +
  labs(title = "Fitted Curves for All Colonies",
       x = "PAR",
       y = "Fitted Pc") +
  theme_minimal() +
  scale_color_manual(values = c(
    "Montipora capitata" = "green",
    "Pocillopora acuta" = "cyan",
    "Porites compressa" = "orange")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"))

# Save the combined plot
ggsave(filename = "output/all_colonies_fitted_curves_missingTempsRemoved.png")



combined_params %>% 
  filter(colony_id != "PACT-C1") %>% 
  ggplot(aes(x = PAR, y = fitted_curve, color = Species, group = colony_id)) + 
  geom_point() +
  geom_line() +
  #geom_text(aes(label = colony_id), hjust = 1.5, vjust = 1.5, size = 3) +  # Add colony_id labels
  facet_wrap(~ Temp.Cat.x, ncol = 6, labeller = label_value) +
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
  )) +
  scale_y_continuous(breaks = function(limits) pretty(limits, n = 6)) +
  labs(y = expression(paste(mu, "mol O"[2], " ", cm^{-2}, " ", h^{-1})))

##### label colonies
combined_params %>% 
  group_by(colony_id) %>% 
  mutate(label_position = 492) %>%  # Choose the position for the label 
  ggplot(aes(x = PAR, y = fitted_curve, color = Species, group = colony_id)) + 
  geom_point() +
  geom_line() +
  geom_text(
    data = . %>% filter(PAR == label_position),  # Filter to label only once per curve
    aes(label = colony_id),
    hjust = -0.2, vjust = 0.5, size = 3
  ) +
  facet_wrap(~ Temp.Cat.x, ncol = 6, labeller = label_value) +
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
  )) +
  labs(y = expression(paste(mu, "mol O"[2], " ", cm^{-2}, " ", h^{-1})))


# Create a directory to save the plots if it doesn't exist
output_dir <- "output/PI_curvePlots_nlsLM_perTempCat_missingTempsRemoved"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate and save one plot per Temp.Cat with all three species together
unique_temp_cat <- unique(combined_params$Temp.Cat.x)

for (temp_cat in unique_temp_cat) {
  plot <- combined_params %>%
    filter(Temp.Cat.x == temp_cat) %>%
    ggplot(aes(x = PAR, y = Pc, color = Species, group = colony_id)) +
    geom_point() +
    geom_line(aes(y = fitted_curve)) +
    scale_color_manual(values=c("green", "cyan", "orange"))+
    #geom_text(aes(label = colony_id), hjust = 1.5, vjust = 1.5, size = 3) +
    labs(title = paste("Temperature Category:", temp_cat),
         x = expression(paste('PAR (', mu, "mol photons m"^-2, 's'^-1,")")),
         y = expression(paste('Photosynthetic rate (', mu, "mol cm"^-2, 's'^-1,")"))) +
    theme_classic() +
    theme(legend.position = "right")
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/pi_curve_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
}





####### PI curve fitted parameters for each individual 

#pars <- nls_data 
#pars <- nls_data 
pars <- colony_params
md <- read_csv("data/1_pi_curves/coral_metadata.csv") 

#### remove missing temperature samples
# Remove rows where the "Run" column has values from 2 to 10
md <- md %>% 
  filter(!(Run %in% 2:10))
df <- left_join(pars, md) 
df <- df %>% 
  pivot_longer(cols=Am:Ik, names_to="variable") 


#### PACT-D6 and Pacu-C1 have very weird PI curve plots, I will remove them from the dataset

library(dplyr)

df <- df %>% 
  filter(!grepl('PACT-D6', colony_id))
df <- df %>% 
  filter(!grepl('PACT-C1', colony_id))

# select one variable at the time
df_Am <- filter(df, variable == "Am")
df_AQY <- filter(df, variable == "AQY")
df_theta <- filter(df, variable == "theta")
df_Ik <- filter(df, variable == "Ik")
df_Rd <- filter(df, variable == "Rd")


###### Plot before outliers removal


# # Facet grid for each variable 
# 
# PI.params <- df %>% 
#   ggplot(aes(x = species, y = value, color =species)) + 
#   geom_boxplot(alpha = 1) + 
#   #scale_color_manual(values=c("green", "cyan", "orange"))+ 
#   facet_wrap(~variable*Temp.Cat, scales = "free_y", nrow=5) + 
#   geom_jitter(aes(color=species), width = 0.1)+ 
#   scale_color_manual(values=c("green", "cyan", "orange"))+
#   theme_classic()+ 
#   theme(legend.position="none")+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                               "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params 
# 
# 
# # Facet grid for each variable separately
# 
# #### theta
# PI.params_theta <- df_theta %>% 
#   ggplot(aes(x = species, y = value, color = species)) + 
#   geom_boxplot() + 
#   facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green", "cyan", "orange"))+
#   theme_classic()+ 
#   theme(legend.position="none")+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                             "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_theta 
# 
# 
# # each species separate across temperatures
# 
# df_filtered_Pacu <- df_theta %>% filter(species == "Pocillopora acuta")
# # Ensure Temp.Cat is a factor
# df_filtered_Pacu$Temp.Cat <- as.factor(df_filtered_Pacu$Temp.Cat)
# 
# PI.params_theta <- df_filtered_Pacu %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("cyan"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "theta")
#   #scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                            # "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# 
# 
# df_filtered_Mcap <- df_theta %>% filter(species == "Montipora capitata")
# # Ensure Temp.Cat is a factor
# df_filtered_Mcap$Temp.Cat <- as.factor(df_filtered_Mcap$Temp.Cat)
# 
# PI.params_theta <- df_filtered_Mcap %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "theta")
# PI.params_theta 
# 
# df_filtered_Pcomp <- df_theta %>% filter(species == "Porites compressa")
# # Ensure Temp.Cat is a factor
# df_filtered_Pcomp$Temp.Cat <- as.factor(df_filtered_Pcomp$Temp.Cat)
# 
# PI.params_theta <- df_filtered_Pcomp %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("orange"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "theta")
# PI.params_theta 
# 
# # Create a directory to save the plots if it doesn't exist
# output_dir <- "output/theta_plots"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }
# 
# # Generate and save one plot per Temp.Cat with all species together
# unique_temp_cat <- unique(df_theta$Temp.Cat)
# 
# for (temp_cat in unique_temp_cat) {
#   plot <- df_theta %>%
#     dplyr::filter(Temp.Cat == temp_cat) %>%
#     ggplot(aes(x = species, y = value, color = species)) + 
#     geom_boxplot() + 
#     facet_wrap(~variable*Temp.Cat, nrow=2) + 
#     geom_jitter(aes(color = species), width = 0.1)+ 
#     scale_color_manual(values=c("green", "cyan", "orange"))+
#     theme_classic()+ 
#     theme(legend.position="none")+ 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#     scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                               "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))
#   
#   # Save the plot
#   ggsave(filename = paste0(output_dir, "/pi_curve_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
# }
# 
# 
# #### Am
# PI.params_Am <- df_Am %>% 
#   ggplot(aes(x = species, y = value, color = species)) + 
#   geom_boxplot() + 
#   facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green", "cyan", "orange"))+
#   theme_classic()+ 
#   theme(legend.position="none")+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                             "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_Am
# 
# # each species separate across temperatures
# 
# df_filtered_Pacu <- df_Am %>% filter(species == "Pocillopora acuta")
# # Ensure Temp.Cat is a factor
# df_filtered_Pacu$Temp.Cat <- as.factor(df_filtered_Pacu$Temp.Cat)
# 
# PI.params_Am <- df_filtered_Pacu %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("cyan"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Am")
# #scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
# # "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_Am
# 
# df_filtered_Mcap <- df_Am %>% filter(species == "Montipora capitata")
# # Ensure Temp.Cat is a factor
# df_filtered_Mcap$Temp.Cat <- as.factor(df_filtered_Mcap$Temp.Cat)
# 
# PI.params_theta <- df_filtered_Mcap %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Am")
# PI.params_theta 
# 
# df_filtered_Pcomp <- df_Am %>% filter(species == "Porites compressa")
# # Ensure Temp.Cat is a factor
# df_filtered_Pcomp$Temp.Cat <- as.factor(df_filtered_Pcomp$Temp.Cat)
# 
# PI.params_theta <- df_filtered_Pcomp %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("orange"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Am")
# PI.params_theta 
# 
# # Create a directory to save the plots if it doesn't exist
# output_dir <- "output/Am_plots"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }
# 
# 
# for (temp_cat in unique_temp_cat) {
#   plot <- df_Am %>%
#     dplyr::filter(Temp.Cat == temp_cat) %>%
#     ggplot(aes(x = species, y = value, color = species)) + 
#     geom_boxplot() + 
#     facet_wrap(~variable*Temp.Cat, nrow=2) + 
#     geom_jitter(aes(color = species), width = 0.1)+ 
#     scale_color_manual(values=c("green", "cyan", "orange"))+
#     theme_classic()+ 
#     theme(legend.position="none")+ 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#     scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                               "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))
#   
#   # Save the plot
#   ggsave(filename = paste0(output_dir, "/pi_curve_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
# }
# 
# 
# ############ AQY
# PI.params_AQY <- df_AQY %>% 
#   ggplot(aes(x = species, y = value, color = species)) + 
#   geom_boxplot() + 
#   facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green", "cyan", "orange"))+
#   theme_classic()+ 
#   theme(legend.position="none")+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                             "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_AQY
# 
# # each species separate across temperatures
# 
# df_filtered_Pacu <- df_AQY %>% filter(species == "Pocillopora acuta")
# # Ensure Temp.Cat is a factor
# df_filtered_Pacu$Temp.Cat <- as.factor(df_filtered_Pacu$Temp.Cat)
# 
# PI.params_AQY <- df_filtered_Pacu %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("cyan"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "AQY")
# #scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
# # "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_AQY
# 
# df_filtered_Mcap <- df_AQY %>% filter(species == "Montipora capitata")
# # Ensure Temp.Cat is a factor
# df_filtered_Mcap$Temp.Cat <- as.factor(df_filtered_Mcap$Temp.Cat)
# 
# PI.params_AQY <- df_filtered_Mcap %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "AQY")
# PI.params_AQY
# 
# df_filtered_Pcomp <- df_AQY %>% filter(species == "Porites compressa")
# # Ensure Temp.Cat is a factor
# df_filtered_Pcomp$Temp.Cat <- as.factor(df_filtered_Pcomp$Temp.Cat)
# 
# PI.params_AQY <- df_filtered_Pcomp %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("orange"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "AQY")
# PI.params_AQY 
# 
# # Create a directory to save the plots if it doesn't exist
# output_dir <- "output/AQY_plots"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }
# 
# 
# for (temp_cat in unique_temp_cat) {
#   plot <- df_AQY %>%
#     dplyr::filter(Temp.Cat == temp_cat) %>%
#     ggplot(aes(x = species, y = value, color = species)) + 
#     geom_boxplot() + 
#     facet_wrap(~variable*Temp.Cat, nrow=2) + 
#     geom_jitter(aes(color = species), width = 0.1)+ 
#     scale_color_manual(values=c("green", "cyan", "orange"))+
#     theme_classic()+ 
#     theme(legend.position="none")+ 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#     scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                               "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))
#   
#   # Save the plot
#   ggsave(filename = paste0(output_dir, "/pi_curve_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
# }
# 
# 
# ################## Ik
# PI.params_Ik <- df_Ik %>% 
#   ggplot(aes(x = species, y = value, color = species)) + 
#   geom_boxplot() + 
#   facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green", "cyan", "orange"))+
#   theme_classic()+ 
#   theme(legend.position="none")+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                             "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_Ik
# 
# # each species separate across temperatures
# 
# df_filtered_Pacu <- df_Ik %>% filter(species == "Pocillopora acuta")
# # Ensure Temp.Cat is a factor
# df_filtered_Pacu$Temp.Cat <- as.factor(df_filtered_Pacu$Temp.Cat)
# 
# PI.params_Ik <- df_filtered_Pacu %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("cyan"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Ik")
# #scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
# # "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_Ik
# 
# df_filtered_Mcap <- df_Ik %>% filter(species == "Montipora capitata")
# # Ensure Temp.Cat is a factor
# df_filtered_Mcap$Temp.Cat <- as.factor(df_filtered_Mcap$Temp.Cat)
# 
# PI.params_Ik <- df_filtered_Mcap %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Ik")
# PI.params_Ik
# 
# df_filtered_Pcomp <- df_Ik %>% filter(species == "Porites compressa")
# # Ensure Temp.Cat is a factor
# df_filtered_Pcomp$Temp.Cat <- as.factor(df_filtered_Pcomp$Temp.Cat)
# 
# PI.params_Ik <- df_filtered_Pcomp %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("orange"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Ik")
# PI.params_Ik 
# 
# 
# # Create a directory to save the plots if it doesn't exist
# output_dir <- "output/Ik_plots"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }
# 
# for (temp_cat in unique_temp_cat) {
#   plot <- df_Ik %>%
#     dplyr::filter(Temp.Cat == temp_cat) %>%
#     ggplot(aes(x = species, y = value, color = species)) + 
#     geom_boxplot() + 
#     facet_wrap(~variable*Temp.Cat, nrow=2) + 
#     geom_jitter(aes(color = species), width = 0.1)+ 
#     scale_color_manual(values=c("green", "cyan", "orange"))+
#     theme_classic()+ 
#     theme(legend.position="none")+ 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#     scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                               "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))
#   
#   # Save the plot
#   ggsave(filename = paste0(output_dir, "/pi_curve_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
# }
# 
# 
# ############### Rd
# PI.params_Rd <- df_Rd %>% 
#   ggplot(aes(x = species, y = value, color = species)) + 
#   geom_boxplot() + 
#   facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green", "cyan", "orange"))+
#   theme_classic()+ 
#   theme(legend.position="none")+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                             "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_Rd
# 
# 
# # each species separate across temperatures
# 
# df_filtered_Pacu <- df_Rd %>% filter(species == "Pocillopora acuta")
# # Ensure Temp.Cat is a factor
# df_filtered_Pacu$Temp.Cat <- as.factor(df_filtered_Pacu$Temp.Cat)
# 
# PI.params_Rd <- df_filtered_Pacu %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("cyan"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Rd")
# #scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
# # "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
# PI.params_Rd
# 
# df_filtered_Mcap <- df_Rd%>% filter(species == "Montipora capitata")
# # Ensure Temp.Cat is a factor
# df_filtered_Mcap$Temp.Cat <- as.factor(df_filtered_Mcap$Temp.Cat)
# 
# PI.params_Rd <- df_filtered_Mcap %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("green"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Rd")
# PI.params_Rd
# 
# df_filtered_Pcomp <- df_Rd %>% filter(species == "Porites compressa")
# # Ensure Temp.Cat is a factor
# df_filtered_Pcomp$Temp.Cat <- as.factor(df_filtered_Pcomp$Temp.Cat)
# 
# PI.params_Rd <- df_filtered_Pcomp %>% 
#   ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
#   geom_boxplot() + 
#   #facet_wrap(~variable*Temp.Cat, nrow=2) + 
#   geom_jitter(aes(color = species), width = 0.1)+ 
#   scale_color_manual(values=c("orange"))+
#   theme_classic()+ 
#   #theme(legend.position="none")+ 
#   #theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
#   labs(y = "Rd")
# PI.params_Rd 
# 
# 
# # Create a directory to save the plots if it doesn't exist
# output_dir <- "output/Rd_plots"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }
# 
# 
# for (temp_cat in unique_temp_cat) {
#   plot <- df_Rd %>%
#     dplyr::filter(Temp.Cat == temp_cat) %>%
#     ggplot(aes(x = species, y = value, color = species)) + 
#     geom_boxplot() + 
#     facet_wrap(~variable*Temp.Cat, nrow=2) + 
#     geom_jitter(aes(color = species), width = 0.1)+ 
#     scale_color_manual(values=c("green", "cyan", "orange"))+
#     theme_classic()+ 
#     theme(legend.position="none")+ 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#     scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
#                               "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))
#   
#   # Save the plot
#   ggsave(filename = paste0(output_dir, "/pi_curve_plot_", temp_cat, ".pdf"), plot = plot, width = 10, height = 6)
# }




###### Load rate data and metadata 
 
#pars <- read_csv("output/pi_curve_pars_nls.csv") 
pars <- read_csv("output/pi_curve_fitted_metrics_Falkowsky.csv") 
md <- read_csv("data/1_pi_curves/coral_metadata.csv") 
df <- left_join(pars, md) 

#### PACT-D6 and Pacu-C1 have very weird PI curve plots, I will remove them from the dataset
df <- df %>% 
  filter(!grepl('PACT-D6', colony_id))
df <- df %>% 
  filter(!grepl('PACT-C1', colony_id))


df$Am.log10 <-log10(df$Am+1) 
df$AQY.log10 <-log10(df$AQY+1) 
df$Rd.log10 <-log10(df$Rd+1) 
df$Ik.log10 <-log10(df$Ik+1) 
df$Ic.log10 <-log10(df$Ic+1) 


# add a group name by pasting Temperature and Treatment 
df$group <- paste0(df$species,"_", df$Temp.Cat) 
#df$id_group <- paste0(df$colony_id,"_", df$Temp.Cat) 

# df2 <- df %>% 
#   pivot_longer(cols=Am:Ic, names_to="variable") 

df2 <- df %>% 
  pivot_longer(cols=Am.log10:Ic.log10, names_to="variable") 

#Visualize data  
r_plot<-df2 %>% 
  ggplot(., aes(x = Temp.Cat, y = value, colour=species)) + 
  #geom_point(aes(fill=species, group=species), pch = 21, size=1, alpha=0.5) + 
  geom_point(aes(fill=species)) + 
  xlab("Temperature") +  
  facet_wrap(vars(variable), scales = "free_y", nrow = 2, strip.position = "top") + 
  theme_classic() +  
  theme(legend.position="none", 
        axis.title=element_text(face="bold", size=16), 
        axis.text=element_text(size=12, color="black"),  
        legend.title=element_text(face="bold", size=14),  
        legend.text=element_text(size=12)); r_plot 


############ View and remove outliers 

#identify outliers by Temperature and Treatment groups 
library(ggstatsplot)
library(pals)

outlier.plot <- ggbetweenstats(df2, group, AQY, outlier.tagging = TRUE, 
                               package = "pals", palette = "polychrome") 

outlier.plot <- ggbetweenstats(df2, group, Ik, outlier.tagging = TRUE, 
                               package = "pals", palette = "polychrome") 

# outlier.plot <- ggbetweenstats(df, group, AQY.log10, outlier.tagging = TRUE, 
#                                plot.type = "box",
#                                outlier.coef = 2,
#                                outlier.label = title,
#                                map_signif_level = TRUE,
#                                results.subtitle = TRUE,
#                                package = "pals", palette = "polychrome"
#                                ) 

#ggsave("output/Pmax_outliersbygroup.pdf", outlier.plot, dpi=300, w=16, h=8, units="in") 

#set quantile values 
q <- c(0.25, 0.75) 
  
# calculate quantile values by Temperature and Treatment groups 
# Quants_AQY <- df2 %>%
#   group_by(species, Temp.Cat) %>%
#   summarise(quant25 = quantile(AQY.log10, probs = q[1]),   #make sure to use summarise not summarize!
#             quant75 = quantile(AQY.log10, probs = q[2]),
#             IQRbyGroup=IQR(AQY.log10))

Quants_AQY <- df2 %>%
  group_by(species, Temp.Cat) %>%
  summarise(quant25 = quantile(AQY, probs = q[1]),   #make sure to use summarise not summarize!
            quant75 = quantile(AQY, probs = q[2]),
            IQRbyGroup=IQR(AQY))

# Quants_Am <- df2 %>%
#   group_by(species, Temp.Cat) %>%
#   summarise(quant25 = quantile(Am.log10, probs = q[1]),   #make sure to use summarise not summarize!
#             quant75 = quantile(Am.log10, probs = q[2]),
#             IQRbyGroup=IQR(Am.log10))

Quants_Am <- df2 %>%
  group_by(species, Temp.Cat) %>%
  summarise(quant25 = quantile(Am, probs = q[1], na.rm = TRUE),   #make sure to use summarise not summarize!
            quant75 = quantile(Am, probs = q[2], na.rm = TRUE),
            IQRbyGroup=IQR(Am, na.rm = TRUE))

# Quants_Ik <- df2 %>%
#   group_by(species, Temp.Cat) %>%
#   summarise(quant25 = quantile(Ik.log10, probs = q[1]),   #make sure to use summarise not summarize!
#             quant75 = quantile(Ik.log10, probs = q[2]),
#             IQRbyGroup=IQR(Ik.log10))

Quants_Ik <- df2 %>%
  group_by(species, Temp.Cat) %>%
  summarise(quant25 = quantile(Ik, probs = q[1], na.rm = TRUE),   #make sure to use summarise not summarize!
            quant75 = quantile(Ik, probs = q[2], na.rm = TRUE),
            IQRbyGroup=IQR(Ik, na.rm = TRUE))

Quants_Rd <- df2 %>%
  group_by(species, Temp.Cat) %>%
  summarise(quant25 = quantile(Rd, probs = q[1], na.rm = TRUE),   #make sure to use summarise not summarize!
            quant75 = quantile(Rd, probs = q[2], na.rm = TRUE),
            IQRbyGroup=IQR(Rd, na.rm = TRUE))

# Quants_Rd <- df2 %>%
#   group_by(species) %>%
#   summarise(quant25 = quantile(Rd, probs = q[1], na.rm = TRUE),   #make sure to use summarise not summarize!
#             quant75 = quantile(Rd, probs = q[2], na.rm = TRUE),
#             IQRbyGroup=IQR(Rd, na.rm = TRUE))


# add a group name by pasting Temperature and Treatment
Quants_AQY$group <-paste0(Quants_AQY$species,"_", Quants_AQY$Temp.Cat)
Quants_Ik$group <-paste0(Quants_Ik$species,"_", Quants_Ik$Temp.Cat)
Quants_Am$group <-paste0(Quants_Am$species,"_", Quants_Am$Temp.Cat)
Quants_Rd$group <-paste0(Quants_Rd$species,"_", Quants_Rd$Temp.Cat)

#Calculate Quantile upper and lower ranges
Quants_AQY$upper <-  Quants_AQY$quant75+1.5*Quants_AQY$IQRbyGroup # Upper Range
Quants_AQY$lower <- Quants_AQY$quant25-1.5*Quants_AQY$IQRbyGroup # Lower Range

Quants_Am$upper <-  Quants_Am$quant75+1.5*Quants_Am$IQRbyGroup # Upper Range
Quants_Am$lower <- Quants_Am$quant25-1.5*Quants_Am$IQRbyGroup # Lower Range

Quants_Ik$upper <-  Quants_Ik$quant75+1.5*Quants_Ik$IQRbyGroup # Upper Range
Quants_Ik$lower <- Quants_Ik$quant25-1.5*Quants_Ik$IQRbyGroup # Lower Range

Quants_Rd$upper <-  Quants_Rd$quant75+1.5*Quants_Rd$IQRbyGroup # Upper Range
Quants_Rd$lower <- Quants_Rd$quant25-1.5*Quants_Rd$IQRbyGroup # Lower Range


### remove outliers from rates
# x1 <- df3 %>%
#   filter(AQY.log10 < upper) %>%
#   filter(AQY.log10 > lower) #13 colonies removed

# #join outlier cutoffs with rate data 
#df <- left_join(df, Quants_AQY, by="group") 
#df2 <- left_join(df2, Quants_AQY, by="group") 

x1 <- df3 %>%
  filter(AQY< upper) %>%
  filter(AQY > lower) 

x2 <- df4 %>%
  filter(Am< upper) %>%
  filter(Am> lower) 

x3 <- df5 %>%
  filter(Ik< upper) %>%
  filter(Ik > lower)

x4 <- df6 %>%
  filter(Rd< upper) %>%
  filter(Rd > lower)

# Create a table with the values that were filtered out
filtered_out_AQY <- anti_join(df2, x1, by = "AQY")
filtered_out_Am <- anti_join(df2, x2, by = "Am")
filtered_out_Ik <- anti_join(df2, x3, by = "Ik")
filtered_out_Rd <- anti_join(df2, x4, by = "Rd")


# Define a function to identify outliers
# identify_outliers <- function(data, value_col) {
#   data %>%
#     group_by(group) %>%
#     mutate(
#       Q1 = quantile(!!sym(value_col), 0.25, na.rm = TRUE),
#       Q3 = quantile(!!sym(value_col), 0.75, na.rm = TRUE),
#       IQR = Q3 - Q1,
#       lower_bound = Q1 - 1.5 * IQR,
#       upper_bound = Q3 + 1.5 * IQR,
#       is_outlier = (!!sym(value_col) < lower_bound) | (!!sym(value_col) > upper_bound)
#     )
# }
# 
# # Apply the function to identify outliers in df2$Am
# df2_outliers <- identify_outliers(df2, "Am")
# 
# # Separate outliers and non-outliers
# outliers <- df2_outliers %>% filter(is_outlier)
# non_outliers <- df2_outliers %>% filter(!is_outlier)


r_plot<- x1 %>%
  #dplyr::filter(variable == "AQY") %>%
  ggplot(., aes(x = Temp.Cat.x, y = AQY, color = species.x)) + 
  geom_point(aes(color=species.x, group=species.x), pch = 21, size=2) +
  xlab("Temperature") +
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  #theme(legend.position="none")+ 
  theme(
    axis.title=element_text(face="bold", size=16),
    axis.text=element_text(size=12, color="black"),
    legend.title=element_text(face="bold", size=14),
    legend.text=element_text(size=12)); r_plot


r_plot<- x2 %>%
  ggplot(., aes(x = Temp.Cat.x, y = Am, color = species.x)) + 
  geom_point(aes(color=species.x, group=species.x), pch = 21, size=2) +
  xlab("Temperature") +
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  #theme(legend.position="none")+ 
  theme(
    axis.title=element_text(face="bold", size=16),
    axis.text=element_text(size=12, color="black"),
    legend.title=element_text(face="bold", size=14),
    legend.text=element_text(size=12)); r_plot


r_plot<- x4 %>%
  #dplyr::filter(variable == "AQY") %>%
  ggplot(., aes(x = Temp.Cat.x, y =Rd, color = species.x)) + 
  geom_point(aes(color=species.x, group=species.x), pch = 21, size=2) +
  xlab("Temperature") +
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  #theme(legend.position="none")+ 
  theme(
    axis.title=element_text(face="bold", size=16),
    axis.text=element_text(size=12, color="black"),
    legend.title=element_text(face="bold", size=14),
    legend.text=element_text(size=12)); r_plot



r_plot<- x3 %>%
  #dplyr::filter(variable == "AQY") %>%
  ggplot(., aes(x = Temp.Cat.x, y =Ik, color = species.x)) + 
  geom_point(aes(color=species.x, group=species.x), pch = 21, size=2) +
  xlab("Temperature") +
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  #theme(legend.position="none")+ 
  theme(
    axis.title=element_text(face="bold", size=16),
    axis.text=element_text(size=12, color="black"),
    legend.title=element_text(face="bold", size=14),
    legend.text=element_text(size=12)); r_plot



#set outliers as NAs

# Replace all values in the AQY column with NA
# filtered_out_AQY$AQY <- NA
# filtered_out_Am$Am <- NA
# filtered_out_Ik$Ik <- NA
# filtered_out_Rd$Rd <- NA
# 
# 
# ######## re-make boxplots after outlier removal
# 
# # Perform the join and ensure values in filtered_out_AQY replace those in df2
# df2_updated <- df2 %>%
#   rows_update(filtered_out_Am, by = c(names(filtered_out_Am)[1], names(filtered_out_Am)[14]))
# 
# df2_updated <- df2_updated %>%
#   rows_update(filtered_out_AQY, by = c(names(filtered_out_AQY)[1], names(filtered_out_AQY)[14]))
# 
# df2_updated <- df2_updated %>%
#   rows_update(filtered_out_Ik, by = c(names(filtered_out_Ik)[1], names(filtered_out_Ik)[14]))
# 
# df2_updated <- df2_updated %>%
#   rows_update(filtered_out_Rd, by = c(names(filtered_out_Rd)[1], names(filtered_out_Rd)[14]))


# Set AQY values in df2 to NA for matching rows in filtered_out_AQY
df2_f <- df2 %>%
  mutate(AQY = if_else(row_number() %in% which(AQY %in% filtered_out_AQY$AQY), NA_real_, AQY))

# # Set Am values in df2 to NA for matching rows in filtered_out_Am based on the first and 14th columns
# df2_f  <- df2_f  %>%
#   mutate(Am = if_else(
#     paste0(.[[1]], .[[14]]) %in% paste0(filtered_out_Am[[1]], filtered_out_Am[[14]]),
#     NA_real_,
#     Am
#   ))
# 
# df2_f  <- df2_f  %>%
#   mutate(Ik = if_else(
#     paste0(.[[1]], .[[14]]) %in% paste0(filtered_out_Ik [[1]], filtered_out_Ik [[14]]),
#     NA_real_,
#     Ik 
#   ))
# 
# df2_f  <- df2_f  %>%
#   mutate(Rd= if_else(
#     paste0(.[[1]], .[[14]]) %in% paste0(filtered_out_Rd [[1]], filtered_out_Rd [[14]]),
#     NA_real_,
#     Rd 
#   ))

### set to NAs values per each parameters that are outliers
# Define the columns and their corresponding filtered data frames
columns_to_update <- c("AQY", "Am", "Ik", "Rd")
filtered_data <- list(filtered_out_AQY, filtered_out_Am, filtered_out_Ik, filtered_out_Rd)

# Loop through each column and update the values in df2_f
for (i in seq_along(columns_to_update)) {
  column <- columns_to_update[i]
  filtered <- filtered_data[[i]]
  
  df3_f <- df2_f %>%
    mutate(!!sym(column) := if_else(
      paste0(.[[1]], .[[14]]) %in% paste0(filtered[[1]], filtered[[14]]),
      NA_real_,
      .[[column]]
    ))
}


# filtered_out_AQY <- filtered_out_AQY %>% drop_na(value)
# merged_all <- full_join(filtered_out_AQY, filtered_out_Am, by = names(filtered_out_AQY)[1:14], suffix = c(".AQY", ".Am")) %>%
#   full_join(filtered_out_Ik, by = names(filtered_out_AQY)[1:14], suffix = c("", ".Ik")) %>%
#   full_join(filtered_out_Rd, by = names(filtered_out_AQY)[1:14], suffix = c("", ".Rd"))
# 
# ##### 13 samples have been filtered out as outliers
# 
# ## get new table with these outliers removed
# df2_outlier_filtered <- anti_join(df2, merged_all, by = names(df2)[1:14])


# Revert pivot_longer to a wide format dataframe
# df2_outlier_filtered <- df2_outlier_filtered %>% 
#   pivot_wider(names_from = variable, values_from = value)

df3_f <- df3_f %>% 
  pivot_wider(names_from = variable, values_from = value)

# df2_outlier_filtered <- df2_outlier_filtered %>% 
#   pivot_longer(cols=Am:Ic, names_to="variable") 

df3_f <- df3_f %>% 
  pivot_longer(cols=Am:Ic, names_to="variable") 

df_AQY_f <- filter(df3_f, variable == "AQY")
df_Am_f <- filter(df3_f, variable == "Am")
df_Ik_f <- filter(df3_f, variable == "Ik")
df_Rd_f <- filter(df3_f, variable == "Rd")


# Boxplots in ggplot2 automatically calculate outliers based on the interquartile range (IQR) of the data 
# being plotted. Even after filtering, the boxplot may identify new outliers based on the remaining data. 


#### AQY
PI.params_AQY_f <- df_AQY_f %>% 
  ggplot(aes(x = species, y = value, color = species)) + 
  #geom_boxplot(outlier.shape = NA) +  # Suppress outliers in the boxplot
  geom_boxplot() +
  #geom_text(aes(label = colony_id), hjust = 0.5, vjust = 1, size = 2) +  # Add colony_id labels
  facet_wrap(~variable*Temp.Cat, nrow=2, scales = "free_y") + 
  #geom_jitter(aes(color = species), width = 0.1)+ 
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
                            "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
PI.params_AQY_f
## Pcom-B6 turns out as a re-calculated outlier
#Pcom-D6
#Mcap-B6

#### Am
PI.params_Am_f <- df_Am_f %>% 
  ggplot(aes(x = species, y = value, color = species)) + 
  #geom_boxplot(outlier.shape = NA) +  # Suppress outliers in the boxplot
  geom_boxplot() +
  facet_wrap(~Temp.Cat, nrow=2, scales = "free_y") + 
  #geom_text(aes(label = colony_id), hjust = 0.5, vjust = 1, size = 2) +  # Add colony_id labels
  #geom_jitter(aes(color = species), width = 0.1)+ 
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
                            "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
PI.params_Am_f
## Pcom-B6 turns out as a re-calculated outlier


#### Ik
PI.params_Ik_f <- df_Ik_f %>% 
  ggplot(aes(x = species, y = value, color = species)) + 
  #geom_boxplot(outlier.shape = NA) +  # Suppress outliers in the boxplot
  geom_boxplot() +
  facet_wrap(~Temp.Cat, nrow=2, scales = "free_y") + 
  #geom_text(aes(label = colony_id), hjust = 0.5, vjust = 1, size = 2) +  # Add colony_id labels
  #geom_jitter(aes(color = species), width = 0.1)+ 
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
                            "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
PI.params_Ik_f
#Pcom-B8

#### Rd
PI.params_Rd_f <- df_Rd_f %>% 
  ggplot(aes(x = species, y = value, color = species)) + 
  geom_boxplot() + 
  facet_wrap(~Temp.Cat, nrow=2, scales = "free_y") + 
  #geom_text(aes(label = colony_id), hjust = 0.5, vjust = 1, size = 2) +  # Add colony_id labels
  #geom_jitter(aes(color = species), width = 0.1)+ 
  scale_color_manual(values=c("green", "cyan", "orange"))+
  theme_classic()+ 
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_x_discrete(labels=c("Montipora capitata" = "Mcap", 
                            "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom")) 
PI.params_Rd_f
## Pacu-H6, Mcap-B7, Pcom-G1 turns out as a re-calculated outliers
#I won't remove Pcom-G1 as it will then re-calculate 2 more outliers (too few data points left then)


##### I'm setting to NAs the re-calculated NAs, run before "pivot_wider"

df3_f <- df3_f %>%
  mutate(
    AQY = if_else(colony_id == "PCOM-B6", NA_real_, AQY),
    Am = if_else(colony_id == "PCOM-B6", NA_real_, Am),
    AQY = if_else(colony_id == "PCOM-D6", NA_real_, AQY),
    AQY = if_else(colony_id == "MCAP-B6", NA_real_, AQY),
    Ik = if_else(colony_id == "PCOM-B8", NA_real_, Ik),
  )

df3_f <- df3_f %>%
  mutate(
    Rd = if_else(colony_id == "PACT-H6", NA_real_, Rd),
    Rd = if_else(colony_id == "MCAP-B7", NA_real_, Rd)
      )


##################### Facet grid for each variable 
# Filter the dataframe to include only the desired variables
df2_outlier_filtered_2 <- df3_f %>%
  filter(variable %in% c("Am", "AQY", "Ik", "Rd"))

# Plot the filtered variables
PI.params <- df2_outlier_filtered_2 %>%
  ggplot(aes(x = species, y = value, color = species, fill = species)) + 
  geom_boxplot(alpha = 0.1) + 
  facet_wrap(~variable*Temp.Cat, scales = "free_y", nrow = 4) + 
  geom_jitter(aes(color = species), width = 0.1) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                              "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))

# Display the plot
PI.params


#https://z3tt.github.io/beyond-bar-and-box-plots/
g <- df2_outlier_filtered_2 %>%
  ggplot(aes(x = species, y = value, color = species, fill = species)) + 
  facet_wrap(~variable*Temp.Cat, scales = "free_y", nrow = 4) + 
  scale_color_manual(values = c("green", "cyan", "orange")) +
  scale_fill_manual(values = c("green", "cyan", "orange")) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_x_discrete(labels = c("Montipora capitata" = "Mcap", 
                              "Pocillopora acuta" = "Pacu", "Porites compressa" = "Pcom"))

g + 
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
    size = 2.5, stroke = .9, shape = 1, color = "white"
  )



## Plot each variable per each species separate across temperatures


# Define the species and colors
species_list <- c("Montipora capitata", "Pocillopora acuta", "Porites compressa")
colors <- c("Montipora capitata" = "green", "Pocillopora acuta" = "cyan","Porites compressa" = "orange")

# Create a directory to save the plots if it doesn't exist
output_dir <- "output/PIparams_perSpecies_minusOutliers_newPIfit"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each variable
variables <- c("Am", "AQY", "Ik", "Rd")

for (variable in variables) {
  # Filter the data for the current variable
  df_filtered <- df2_outlier_filtered_2 %>% filter(variable == variable)
  df_filtered$Temp.Cat <- as.factor(df_filtered$Temp.Cat)
  
  # Generate the plot
  plot <- df_filtered %>% 
    ggplot(aes(x = Temp.Cat, y = value, color = species)) + 
    geom_boxplot() + 
    facet_wrap(~species*variable, nrow = 3, scales = "free_y") + 
    geom_jitter(aes(color = species), width = 0.1) + 
    scale_color_manual(values = colors) +
    theme_classic() + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(y = "Value", title = paste("Variable:", variable))
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/plot_", variable, ".pdf"), plot = plot, width = 10, height = 12)
}


###### change graph type

# Define the species and colors
species_list <- c("Montipora capitata", "Pocillopora acuta", "Porites compressa")
colors <- c("Montipora capitata" = "green", "Pocillopora acuta" = "cyan", "Porites compressa" = "orange")

# Create a directory to save the plots if it doesn't exist
output_dir <- "output/PIparams_perSpecies_minusOutliers_newPIfit"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each variable
variables <- c("Am", "AQY", "Ik", "Rd")

for (variable in variables) {
  # Filter the data for the current variable
  df_filtered <- df2_outlier_filtered_2 %>% filter(variable == variable)
  df_filtered$Temp.Cat <- as.factor(df_filtered$Temp.Cat)
  
  # Generate the plot
  plot <- df_filtered %>% 
    ggplot(aes(x = Temp.Cat, y = value, color = species, fill = species)) + 
    facet_wrap(~species*variable, nrow = 3, scales = "free_y") + 
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_classic() + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    labs(y = "Value", title = paste("Variable:", variable)) +
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
      size = 2.5, stroke = .9, shape = 1, color = "white"
    )
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/plot2_", variable, ".pdf"), plot = plot, width = 10, height = 12)
}




###### **Statistical Analysis**     

## Load packages
library(dplyr)
library(ggpubr)
library(onewaytests)
library(purrr)
library(tidyr)


#select relevant columns
df2_outlier_filtered_selected <- df2_outlier_filtered_2 %>% 
  select(1:5, 13:14)%>%
  distinct()

# Revert pivot_longer to a wide format dataframe
df_wide <- df2_outlier_filtered_selected %>% 
  pivot_wider(names_from = variable, values_from = value) ### 130 samples here

library(writexl)
# Save df_wide as an Excel file
write_xlsx(df_wide, "output/df_wide.xlsx")


## Assess the normality of the data per each variable
#Shapiro-Wilk’s method is widely recommended for normality test and it provides better power than K-S. It is based on the correlation between the data and the corresponding normal scores.

#Assess the normality of the data using the Shapiro-Wilk test for each combination of species and 
#variable across Temp.Cat, and for each combination of variable and Temp.Cat across species, using 
#the group_by and nest functions from the dplyr package along with the map function from the purrr package.

# Ensure Temp.Cat is a factor
df_wide$Temp.Cat <- as.factor(df_wide$Temp.Cat)

# Function to perform Shapiro-Wilk test
perform_shapiro_test <- function(data) {
  shapiro_test <- shapiro.test(data$value)
  tidy(shapiro_test)
}

# a) Perform Shapiro-Wilk test per each species, combining all Temp.Cat values for each variable
shapiro_species_variable <- df_wide %>%
  pivot_longer(cols = Am:Ik, names_to = "variable", values_to = "value") %>%
  group_by(species, variable) %>%
  filter(n() >= 3) %>%  # Filter out groups with fewer than 3 observations
  nest() %>%
  mutate(shapiro_test = map(data, perform_shapiro_test)) %>%
  unnest(shapiro_test)

# b) Perform Shapiro-Wilk test per each Temp.Cat, combining all species for each variable
perform_shapiro_test_species <- function(data) {
  shapiro_test <- shapiro.test(data$value)
  tidy(shapiro_test)
}

shapiro_variable_temp_cat <- df_wide %>%
  pivot_longer(cols = Am:Ik, names_to = "variable", values_to = "value") %>%
  group_by(Temp.Cat, variable) %>%
  filter(n() >= 3) %>%
  nest() %>%
  mutate(shapiro_test = map(data, perform_shapiro_test_species)) %>%
  unnest(shapiro_test)

# View the results
shapiro_species_variable  #only AQY normally distributed across species
shapiro_variable_temp_cat

# Save the Shapiro-Wilk test results to CSV files
write_csv(shapiro_species_variable, "output/shapiro_species_variable.csv")
write_csv(shapiro_variable_temp_cat, "output/shapiro_variable_temp_cat.csv")



### Function to generate Q-Q plot
generate_qqplot <- function(data, title) {
  ggqqplot(data$value, title = title)
}

# a) Generate Q-Q plots per each species, combining all Temp.Cat values for each variable
qqplots_species_variable <- df_wide %>%
  pivot_longer(cols = Am:Ik, names_to = "variable", values_to = "value") %>%
  group_by(species, variable) %>%
  filter(n() >= 3) %>%  # Filter out groups with fewer than 3 observations
  nest() %>%
  mutate(qqplot = map(data, ~ generate_qqplot(.x, paste("Q-Q Plot for", unique(.x$species), "-", unique(.x$variable)))))

# b) Generate Q-Q plots per each Temp.Cat, combining all species for each variable
qqplots_variable_temp_cat <- df_wide %>%
  pivot_longer(cols = Am:Ik, names_to = "variable", values_to = "value") %>%
  group_by(Temp.Cat, variable) %>%
  filter(n() >= 3) %>%  # Filter out groups with fewer than 3 observations
  nest() %>%
  mutate(qqplot = map(data, ~ generate_qqplot(.x, paste("Q-Q Plot for Temp.Cat", unique(.x$Temp.Cat), "-", unique(.x$variable)))))

# Save the Q-Q plots to files
output_dir_species_variable <- "output/QQplots/qqplots_species_variable"
output_dir_variable_temp_cat <- "output/QQplots/qqplots_variable_temp_cat"

if (!dir.exists(output_dir_species_variable)) {
  dir.create(output_dir_species_variable, recursive = TRUE)
}

if (!dir.exists(output_dir_variable_temp_cat)) {
  dir.create(output_dir_variable_temp_cat, recursive = TRUE)
}

# Save Q-Q plots for species and variable
for (i in seq_len(nrow(qqplots_species_variable))) {
  plot <- qqplots_species_variable$qqplot[[i]]
  species <- qqplots_species_variable$species[i]
  variable <- qqplots_species_variable$variable[i]
  filename <- paste0(output_dir_species_variable, "/qqplot_", species, "_", variable, ".png")
  ggsave(filename = filename, plot = plot, width = 8, height = 6)
}

# Save Q-Q plots for Temp.Cat and variable
for (i in seq_len(nrow(qqplots_variable_temp_cat))) {
  plot <- qqplots_variable_temp_cat$qqplot[[i]]
  temp_cat <- qqplots_variable_temp_cat$Temp.Cat[i]
  variable <- qqplots_variable_temp_cat$variable[i]
  filename <- paste0(output_dir_variable_temp_cat, "/qqplot_", temp_cat, "_", variable, ".png")
  ggsave(filename = filename, plot = plot, width = 8, height = 6)
}


## perform the Levene test for homogeneity of variances

library(car)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)

# Ensure Temp.Cat is a factor
df_wide$Temp.Cat <- as.factor(df_wide$Temp.Cat)

# Function to perform Levene test
perform_levene_test <- function(data) {
  levene_test <- leveneTest(value ~ Temp.Cat, data = data)
  tidy(levene_test)
}

# a) Perform Levene test per each species, combining all Temp.Cat values for each variable
levene_species_variable <- df_wide %>%
  pivot_longer(cols = Am:Ik, names_to = "variable", values_to = "value") %>%
  group_by(species, variable) %>%
  filter(n() >= 3) %>%  # Filter out groups with fewer than 3 observations
  nest() %>%
  mutate(levene_test = map(data, perform_levene_test)) %>%
  unnest(levene_test)

# b) Perform Levene test per each Temp.Cat, combining all species for each variable
perform_levene_test_species <- function(data) {
  levene_test <- leveneTest(value ~ species, data = data)
  tidy(levene_test)
}

levene_variable_temp_cat <- df_wide %>%
  pivot_longer(cols = Am:Ik, names_to = "variable", values_to = "value") %>%
  group_by(Temp.Cat, variable) %>%
  filter(n() >= 3) %>%
  nest() %>%
  mutate(levene_test = map(data, perform_levene_test_species)) %>%
  unnest(levene_test)

# View the results
levene_species_variable
levene_variable_temp_cat

# Save the Levene test results to CSV files
write_csv(levene_species_variable, "output/levene_species_variable.csv")
write_csv(levene_variable_temp_cat, "output/levene_variable_temp_cat.csv")

# Merge the two data frames based on Temp.Cat and variable
merged_results <- inner_join(shapiro_variable_temp_cat, levene_variable_temp_cat, by = c("Temp.Cat", "variable"))

# Add the assumptions column
merged_results <- merged_results %>%
  mutate(assumptions = if_else(p.value.x > 0.05 & p.value.y > 0.05, "met", "non met"))

# View the merged results
print(merged_results)



####### ANOVA (for normally-distributed data and with homogeneous variance)

#To run ANOVA and Tukey HSD tests only for the combinations of Temp.Cat and variable where the assumptions column has the value "met", 
#you can filter the merged_results table and then loop through the filtered combinations

library(dplyr)
library(broom)

# Filter the merged_results table for rows where assumptions are "met"
valid_combinations <- merged_results %>%
  filter(assumptions == "met") %>%
  select(Temp.Cat, variable)

# Initialize an empty list to store results
anova_tukey_results <- list()

# Loop through each valid combination
for (i in seq_len(nrow(valid_combinations))) {
  temp_cat <- valid_combinations$Temp.Cat[i]
  variable <- valid_combinations$variable[i]
  
  # Filter the data for the current Temp.Cat and variable
  filt_data <- df_wide %>%
    filter(Temp.Cat == temp_cat) %>%
    select(species, !!sym(variable)) %>%
    rename(value = !!sym(variable))  # Rename the selected column to "value" for consistency
  
  # Check if there are enough non-NA values for ANOVA
  if (nrow(filt_data) > 2 && sum(!is.na(filt_data$value)) > 2) {
    # Run ANOVA
    aov_result <- aov(value ~ species, data = filt_data)
    
    # Run Tukey HSD test
    tukey_result <- TukeyHSD(aov_result)
    
    # Store the results
    anova_tukey_results[[paste(temp_cat, variable, sep = "_")]] <- list(
      Temp.Cat = temp_cat,
      variable = variable,
      anova = tidy(aov_result),
      tukey = tidy(tukey_result)
    )
  }
}

# Combine the results into a single data frame for easier viewing
anova_tukey_results_df <- bind_rows(
  lapply(anova_tukey_results, function(x) {
    tibble(
      Temp.Cat = x$Temp.Cat,
      variable = x$variable,
      anova = list(x$anova),
      tukey = list(x$tukey)
    )
  })
)

# View the results
print(anova_tukey_results_df)



#### check one compared group to make sure the loop works

# Filter the data for Temp.Cat == "12" and select the "Am" column
filtered_data_12Am <- df_wide %>%
  filter(Temp.Cat == "12") %>%
  select(species, Am)

# Perform ANOVA
anova_result_12Am  <- aov(Am ~ species, data = filtered_data_12Am)

###Display the summary of the ANOVA
summary(anova_result_12Am)
# Df    Sum Sq   Mean Sq F value   Pr(>F)    
# species      2 1.801e-08 9.004e-09   13.17 0.000416 ***
#   Residuals   16 1.094e-08 6.840e-10                     

# Perform Tukey HSD test
tukey_result_12Am  <- TukeyHSD(anova_result_12Am)

# View the Tukey HSD results
print(tukey_result_12Am)
# $species
# diff           lwr          upr     p adj
# Pocillopora acuta-Montipora capitata -2.446975e-05 -6.200207e-05 1.306257e-05 0.2421162
# Porites compressa-Montipora capitata  5.127790e-05  1.374558e-05 8.881022e-05 0.0074732
# Porites compressa-Pocillopora acuta   7.574765e-05  3.679852e-05 1.146968e-04 0.0003514

# Optionally, plot the Tukey HSD results
plot(tukey_result_12Am)



####### Kruskal test for non normally-distributed data and/or with no homogeneous variance

# Filter the merged_results table for rows where assumptions are "non met"
non_parametric_combinations <- merged_results %>%
  filter(assumptions == "non met") %>%
  select(Temp.Cat, variable)

# Initialize an empty list to store results
non_parametric_results <- list()

# Loop through each combination
for (i in seq_len(nrow(non_parametric_combinations))) {
  temp_cat <- non_parametric_combinations$Temp.Cat[i]
  variable <- non_parametric_combinations$variable[i]
  
  # Filter the data for the current Temp.Cat and variable
  filt_data <- df_wide %>%
    filter(Temp.Cat == temp_cat) %>%
    select(species, !!sym(variable)) %>%
    rename(value = !!sym(variable))  # Rename the selected column to "value" for consistency
  
  # Check if there are enough non-NA values for the test
  if (nrow(filt_data) > 2 && sum(!is.na(filt_data$value)) > 2) {
    # Run Kruskal-Wallis test
    kruskal_result <- kruskal.test(value ~ species, data = filt_data)
    
    # Run pairwise Dunn's test
    dunn_result <- dunnTest(
      filt_data$value, filt_data$species,
      method = "bonferroni"
    )
    
    # Store the results
    non_parametric_results[[paste(temp_cat, variable, sep = "_")]] <- list(
      Temp.Cat = temp_cat,
      variable = variable,
      kruskal = tidy(kruskal_result),
      dunn = dunn_result
    )
  }
}

# Combine the results into a single data frame for easier viewing
non_parametric_results_df <- bind_rows(
  lapply(non_parametric_results, function(x) {
    # Extract the Dunn's test results as a tidy data frame
    dunn_df <- as.data.frame(as.table(x$dunn$p.value)) %>%
      rename(group1 = Var1, group2 = Var2, p_value = Freq) %>%  # Rename columns
      mutate(W_value = NA)  # Add a placeholder for W values (not available in dunnTest)
    
    tibble(
      Temp.Cat = x$Temp.Cat,
      variable = x$variable,
      kruskal = list(x$kruskal),
      dunn = list(dunn_df)  # Store the tidy Dunn's results
    )
  })
)

# View the results
print(non_parametric_results_df)




### re-make plot with statistics

library(ggsignif)
library(ggpubr)
library(rstatix)
library(tidyr)

# Prepare the Tukey and Dunn's results for annotation
# Extract significant comparisons from Tukey results
tukey_annotations <- anova_tukey_results_df %>%
  unnest(tukey) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(test = "Tukey") %>%
  separate(contrast, into = c("group1", "group2"), sep = "-") %>%  # Split contrast into group1 and group2
  select(Temp.Cat, variable, group1, group2, adj.p.value, test) %>%
  rename(p.adj = adj.p.value)

# Extract significant comparisons from Dunn's results
dunn_annotations <- non_parametric_results_df %>%
  unnest(dunn) %>%
  filter(p_value < 0.05) %>%
  mutate(test = "Dunn's") %>%
  select(Temp.Cat, variable, group1, group2, p_value, test)%>%
  rename(p.adj = p_value)

# Combine Tukey and Dunn's annotations
annotations <- bind_rows(tukey_annotations, dunn_annotations)

# Add significance levels to the annotations
annotations <- annotations %>%
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Compute y.position values for each combination of variable and Temp.Cat
y_positions <- df2_outlier_filtered_2 %>%
  group_by(variable, Temp.Cat) %>%
  summarise(
    max_value = max(value, na.rm = TRUE),  # Maximum value for the group
    .groups = "drop"
  ) %>%
  mutate(
    y_base = max_value  # Add a base offset to the maximum value
  )

# Ensure Temp.Cat is a factor in both data frames
annotations <- annotations %>%
  mutate(Temp.Cat = as.factor(Temp.Cat))

y_positions <- y_positions %>%
  mutate(Temp.Cat = as.factor(Temp.Cat))

# Add y.position to annotations with dynamic spacing
annotations <- annotations %>%
  left_join(y_positions, by = c("variable", "Temp.Cat")) %>%
  group_by(variable, Temp.Cat) %>%
  mutate(
    y.position = y_base + (row_number() - 1) * (y_base * 0.2)  # Dynamic offset based on y_base
  ) %>%
  ungroup() %>%
  select(-max_value, -y_base)  # Remove intermediate columns


# Filter out rows with missing or non-finite values in the data
df2_outlier_filtered_3 <- df2_outlier_filtered_2 %>%
  filter(!is.na(value) & is.finite(value))

# Add annotations to the plot
g <- df2_outlier_filtered_3 %>%
  ggplot(aes(x = species, y = value, color = species, fill = species)) + 
  facet_wrap(~variable*Temp.Cat, scales = "free_y", nrow = 4) + 
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

# Add statistical annotations to the plot
g + stat_pvalue_manual(
  annotations,
  label = "p.adj.signif",  # Use asterisks for significance
  hide.ns = TRUE
) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))


