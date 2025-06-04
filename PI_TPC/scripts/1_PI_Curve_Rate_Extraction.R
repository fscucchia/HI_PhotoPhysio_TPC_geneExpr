
############# "Oxygen flux rate extractions " #############

## install packages if you dont already have them in your library
if (!require("devtools")) install.packages("devtools")
if (!require("furrr")) install.packages("furrr")
if (!require("future")) install.packages("future")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("lubridate")) install.packages("lubridate")
if (!require("cowplot")) install.packages("cowplot")
if (!require("LoLinR")) install_github('colin-olito/LoLinR') 
install.packages("rTPC")
install.packages("nls.multstart")

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

## libraries for parallel processing
library(future)
library(furrr)

## Import metadata
path.p <- "data/1_pi_curves/" #the location of all your respirometry files 

# List data files
file.names <- list.files(path = path.p, pattern = "csv$")  # list all csv file names in the folder
file.names <- file.names[!grepl("metadata", file.names)]   # omit metadata from files to be read in as data

# Load PI curve sample metadata (i.e., which corals were in which runs)
sample.info <- read_csv(file = "data/1_pi_curves/1_pi_curves_sample_metadata.csv")

# Load PI curve run metadata (i.e., light levels and interval times for each run)
run.info <- read_csv(file = "data/1_pi_curves/1_pi_curves_run_metadata.csv")

# Join all coral and run metadata
metadata <- full_join(sample.info, run.info) %>%
  mutate(Date = as_date(as.character(Date), format = "%Y%m%d", tz = "Atlantic"))

# Select only certain columns
metadata <- metadata %>%
  select(Species, colony_id, Run, Temp.Cat, Chamber.Vol.L, Date, Start.time, Stop.time, Light_Value, Light_Level)


# Read in all data files
#identify the sample id from the file name and read in data
df <- tibble(file.name = file.names) %>%
  mutate(colony_id = gsub("_.*", "", file.name),                              # Get colony_id from filename
          info = map(colony_id, ~filter(metadata, colony_id == .)),           # Get associated sample info
         data0 = map(file.name, ~read_csv(file.path(path.p, .), skip=1, col_types = cols(.default = "d", Time = "t"))))   # Get associated O2 data

# Select only Time, Value, and Temp columns from O2 data
df1 <- df %>%
  mutate(data0 = map(data0, ~select(., Time, Value, Temp)))%>%
  mutate(data0 = map(data0, ~(.x %>% filter(complete.cases(.))))) #remove NAs 


## Use the time breaks in the sample info to link O2 data with light levels
#Use start time of each light step from the metadata to separate data by light stop
df2 <- df1 %>%
  mutate(intervals = map2(data0, info, function(.x, .y) {
    split(.x, f = cut(as.numeric(.x$Time), breaks = as.numeric(c(.y$Start.time, last(.y$Stop.time))),
                      labels = as.character(.y$Light_Value)))})) %>%
  mutate(data = map(intervals, ~ unnest(tibble(.), .id = "Light_Value")))

#df2 <- df1 %>%
#  mutate(intervals = map2(data0, info, function(.x, .y) {
#    split(.x, f = cut(as.numeric(.x$Time), breaks = as.numeric(c(.y$Start.time, last(.y$Stop.time))),
#                      labels = as.character(.y$Light_Value)))})) %>%
#  mutate(data = map(intervals, ~ unnest(tibble(.), cols = any_of("data"), .id = "Light_Value")))


### Thin data
# Set thinning parameter
thin_par <- 20

### will try 10 seconds instead of 20 to see if fitting improves for failed colonies
#thin_par <- 10


# Thin data for all samples
df3 <- df2 %>%
  mutate(thin_data = map(data, ~ slice(., seq(1, nrow(.), thin_par))))

# Create plots for full dataset and thinned data
# df3 <- df3 %>%
#   mutate(data_plot = map2(data, colony_id, ~ ggplot(.x, aes(x = Time, y = Value)) + 
#                             facet_wrap(~ as.numeric(Light_Value), scales = "free") +
#                             geom_point() +
#                             labs(title = .y)),
#     thin_data_plot = map2(thin_data, colony_id, ~ ggplot(.x, aes(x = Time, y = Value)) + 
#                             facet_wrap(~ as.numeric(Light_Value), scales = "free") +
#                             geom_point() +
#                             labs(title = .y)))
# 
# # Example of plots
# cowplot::plot_grid(df3$data_plot[[1]], df3$thin_data_plot[[1]], nrow = 2,
#                    labels = c("Example plot: all data", "Example plot: thinned data"))
# 
# 
# #### The full or thinned data plot for any sample can be accessed like this:
# # (I checked the colonies that failed the NLS fitting (21 of them) and they indeed look weird in these plots)
# 
# 
# df3 %>%
#   filter(colony_id == "Sample-Name-Here") %>%
#   pull(thin_data_plot)
# 
# 
# df3 %>%
#   filter(colony_id == "PACT-E7") %>%
#   pull(thin_data_plot)
# 
# df3 %>%
#   filter(colony_id == "MCAP-C3") %>%
#   pull(thin_data_plot)
# 
# df3 %>%
#   filter(colony_id == "PCOM-H3") %>%
#   pull(thin_data_plot)


# Fit regressions to each interval for each sample
# Define function for fitting LoLinR regressions to be applied to all intervals for all samples
fit_reg <- function(df3) {
  rankLocReg(xall = as.numeric(df3$Time), yall = df3$Value, 
             alpha = 0.2, method = "pc", verbose = FALSE)
}

# fit_reg <- function(df3) {
#   rankLocReg(xall = as.numeric(df3$Time), yall = df3$Value, 
#              alpha = 0.2, method = "pc", verbose = TRUE)
# }


# Setup for parallel processing
future::plan(multisession)

# Map LoLinR function onto all intervals of each sample's thinned dataset
df3 <- df3 %>%
  mutate(regs = furrr::future_map(thin_data, function(.) {       # future_map executes function in parallel
    group_by(., Light_Value) %>%
    do(rankLcRg = fit_reg(.))
  }))

##### EXTRACT MEAN TEMPERATURE OF EACH LIGHT LEVEL

 df4 <- df3 %>%
   mutate(regs = furrr::future_map(thin_data, function(.) {       # future_map executes function in parallel
     group_by(., Light_Value) %>%
       do(rankLcRg = fit_reg(.),
          MeanTemp = summarise(., mean(Temp)),
          SDTemp = summarise(., sd(Temp)),
          SEMTemp = summarise(., sd(Temp) / sqrt(n())))
   }))

thin_data_table <- df4$thin_data

# Extract mean temperature, sd, and sem into a new data frame
temp_stats <- df4 %>%
  unnest(cols = c(regs)) %>%
  select(Light_Value, MeanTemp, SDTemp, SEMTemp) %>%
  unnest(cols = c(MeanTemp, SDTemp, SEMTemp))

# Write the data frame to a new CSV file
write.csv(temp_stats, "temperature_stats_perLightLevel_thinned20.csv", row.names = FALSE)


## Now 'regs' contains the fitted local regressions for each interval of each sample's thinned dataset
#### The diagnostics for any regression can be plotted by specifying a colony_id and the number of the light
#curve interval: Define function to pull out and plot regression diagnostics

plot_rankLcRg <- function(colony_id_input, interval_number) {
  df4 %>%
    filter(colony_id == colony_id_input) %>%
    pluck("regs", 1, "rankLcRg", interval_number) %>%
    plot()
}


pdf("output/PCOMP-D1_test.pdf")
plot_rankLcRg("PCOM-D1", 4)
# plot_rankLcRg("Past-D1", 2)
# plot_rankLcRg("Past-D1", 3)
# plot_rankLcRg("Past-D1", 4)
# plot_rankLcRg("Past-D1", 5)
# plot_rankLcRg("Past-D1", 6)
# plot_rankLcRg("Past-D1", 7)
# plot_rankLcRg("Past-D1", 8)
# plot_rankLcRg("Past-D1", 9)
# plot_rankLcRg("Past-D1", 10)
dev.off()


### Extract slope of best regression for each interval for each sample

#extract slope as rate
df.out <- df4 %>% 
  unnest(regs) %>%
  mutate(micromol.L.s = map_dbl(rankLcRg, ~ pluck(., "allRegs", "b1", 1)))

df.out2 <- df.out %>% 
  unnest(MeanTemp)

df.out2 <- df.out %>% 
  unnest(MeanTemp) %>% 
  unnest(SDTemp) %>%
  unnest(SEMTemp)

#verify sample numbers
unique(df.out2[[2]])
length(unique(df.out2[[2]])) ## 280

#select only the essential columns
xx <- select(df.out2,colony_id, Light_Value, micromol.L.s, "mean(Temp)")

#add a grouping id for each colony at each light level
xx$grouping.id <- paste0(xx$colony_id, "-",xx$Light_Value)
nrow(xx)
nrow(distinct(xx))

#select only the essential columns
mx <- select(metadata, Species, colony_id, Run, Chamber.Vol.L, Temp.Cat)
nrow(mx)
nrow(distinct(df))

#join rates with metadata
pr <- left_join(xx, mx, by="colony_id")
length(unique(pr$colony_id))
nrow(distinct(pr))
pr <-distinct(pr)
  
# Write raw data to output file
write.csv(pr, "/output/pi_curve_extracted_rates_Raw_thinned20.csv")


# Adjust rates by chamber volume, subtract blank, and normalize to surface area

# Rename the column Chamber.Vol.L to Chamber.Vol.mL
pr <- pr %>% 
  rename(Chamber.Vol.mL = Chamber.Vol.L) ### column name was wrong, these are mL not L

# Create a new column Chamber.Vol.L with the values divided by 1000
pr <- pr %>% 
  mutate(Chamber.Vol.L = Chamber.Vol.mL / 1000)

# Correct for chamber volume 
pr <- pr %>% mutate(micromol.s = micromol.L.s * Chamber.Vol.L)
length(unique(pr$colony_id))

# Correct for blank rates
# Get blank values -- average for each run and light value in case multiple blanks
blanks <- pr %>%
  filter(grepl("BK", colony_id)) %>%
  group_by(Temp.Cat, Light_Value) %>%
  summarise(micromol.s.blank=mean(micromol.s))%>%
  mutate(blank_id=paste0(Temp.Cat,"-",Light_Value))

### generate a key for the blank id
pr <- pr %>%
  mutate(blank_id=paste0(Temp.Cat,"-",Light_Value))
length(unique(pr$colony_id))

#plot blank values
blanks %>% ggplot(aes(x=as.numeric(Light_Value), y=micromol.s.blank,colour = as.factor(Temp.Cat)))+
  geom_point()+
  geom_line()

#examine the effects of light and temp on the blank rates
anova(lm(micromol.s.blank~as.factor(Light_Value)*Temp.Cat, data=blanks))

#join the data and the mean of the blanks per temperature for each specific light level
pr  <- left_join(pr ,blanks,by = "blank_id")
length(unique(pr$colony_id))

#subtract temp and light specific blank values from samples
pr <- pr %>%
  mutate(micromol.s.adj = micromol.s - micromol.s.blank) %>%
  # After correcting for blank values, remove blanks from data
  filter(!grepl("BK", colony_id))
length(unique(pr$colony_id))  ######## 240

# Import surface area data
sa <- read.csv("/output/1_surface_area_final.csv",sep = ",")

# Join surface area with rest of data
pr <- left_join(pr, select(sa, colony_id, surface.area.cm2))
length(unique(pr$colony_id))

# pr2  <- dplyr::inner_join(
#    pr ,
#    dplyr::select(sa, -any_of(names(pr)), colony_id),
#    by = "colony_id"
# )
# length(unique(pr2$colony_id))

# Normalize rates by surface area
pr <- pr %>%
  mutate(micromol.cm2.s = micromol.s.adj / surface.area.cm2,
         micromol.cm2.h = micromol.cm2.s * 3600)


length(unique(pr$colony_id))


#### remove missing temperature samples

# Remove rows where the "Run" column has values from 2 to 10
pr <- pr %>% 
  filter(!(Run %in% 2:10))

pr <- pr %>% 
  drop_na()

# Write extracted rates to output file

# Select variables to write to file
pr.out <- pr %>% 
  select(Species, colony_id, Light_Value.x, Temp.Cat.x, "mean(Temp)", Run, micromol.cm2.s, micromol.cm2.h)
length(unique(pr.out$colony_id))    ##### 134 samples (without lost temps)!

# Write to output file
#write.csv(pr.out, "output/pi_curve_extracted_rates_meanTemp_thinned20.csv")
write.csv(pr.out, "output/pi_curve_extracted_rates_missingTempRemoved_thinned20.csv")

#add column of means per each temp cat
colnames(pr.out)[5] ="meanTempByLight"

pr.out_mean <- pr.out %>% 
  group_by(Temp.Cat.x) %>% 
  mutate(Temp.Cat.mean = mean(meanTempByLight))

#with sd and sem
pr.out_mean <- pr.out %>% 
  group_by(Temp.Cat.x) %>% 
  mutate(
    Temp.Cat.mean = mean(meanTempByLight),
    Temp.Cat.sd = sd(meanTempByLight),
    Temp.Cat.sem = sd(meanTempByLight) / sqrt(n())
  )


# Write to output file
write.csv(pr.out_mean, "output/pi_curve_extracted_rates_meanTemp_missingTempRemoved_thinned20.csv")


# Extract mean temperature, sd, and sem per each Temp.Cat.x into a new data frame
temp_stats2 <- pr.out_mean %>%
  group_by(Temp.Cat.x) %>%
  summarise(
    MeanTemp = mean(meanTempByLight),
    SDTemp = sd(meanTempByLight),
    SEMTemp = sd(meanTempByLight) / sqrt(n())
  )

# Write the data frame to a new CSV file
write.csv(temp_stats2, "output/temperature_stats_per_TempCat_missingTempRemoved_thinned20.csv", row.names = FALSE)





