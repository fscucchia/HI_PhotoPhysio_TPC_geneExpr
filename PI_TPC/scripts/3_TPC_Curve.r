
######################### **TPC fitting**   
# Based on Padifeld et al **rTPC and nls.multstart: A new pipeline to fit thermal performance curves in r**   
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13585   
#https://padpadpadpad.github.io/rTPC/

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("devtools")
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library(tidyverse)
library(cowplot)
library(broom)
library(rTPC)
library(nls.multstart)
library(here)
library(nlraa)
library(car)
library(ggrepel)
library(MuMIn)
library(readxl)
#library(future)
#library(furrr)

# Define data
# load PI parameters output from the PI curve analysis
df_wide <- read_excel(here("output", "df_wide.xlsx"))

# choose model 
get_model_names() 

#### Mcap CURVE FIT - Pmax
Mcap.df <- df_wide %>%  
  filter(species=="Montipora capitata") 

Mcap.df$Temp.Cat <- as.numeric(as.character(Mcap.df$Temp.Cat))

# Calculate mean and sd values for Am, grouped by Temp.Cat

# Mcap.df_mean <- Mcap.df %>%
#   group_by(Temp.Cat) %>%
#   summarise(across(7:11, mean, na.rm = TRUE)) %>%
#   mutate(colony_id = "MCAP") %>%
#   select(Temp.Cat, Am, colony_id)

Mcap.df_mean <- Mcap.df %>%
  dplyr::select(Temp.Cat, Am, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Am,na.rm = TRUE), ave_rate = mean(Am,na.rm = TRUE),
            groups = 'drop')

# View the resulting dataframe
print(Mcap.df_mean)

Mcap.df_mean$Temp.Cat <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Mcap.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#                                             # include weights here!
#          pawar = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = Temp.Cat, r_tref, e, eh, topt, tref = 18),
#                                           data = .x,
#                                           iter = c(4,4,4,4),
#                                           start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018') - 10,
#                                           start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018') + 10,
#                                           lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018'),
#                                           upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018'),
#                                           supp_errors = 'Y',
#                                           convergence_count = FALSE,
#                                           modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Mcap.df_mean$Temp.Cat), max(Mcap.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   dplyr::select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Am',
#        title = 'Am across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   dplyr::select(-fit) %>%
#   unnest(info) %>%
#   dplyr::select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model

# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Am',
#        title = 'Am across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


############ Pawar 2018, modified version of Sharpe Schoolfield 1981 model 
#Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory.
#Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

#### we want to run the model taking into account means and standard deviations of the PI parameters
#https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
#https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


Mcap.df_mean <- Mcap.df %>%
  dplyr::select(Temp.Cat, Am, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Am),ave_rate = mean(Am),
            groups = 'drop')

Mcap.df_mean$temp <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #temp needs to be numeric

Mcap.df_mean <- Mcap.df_mean %>%
  mutate(Temp.Cat = temp) %>%
  select(temp, everything(), -Temp.Cat)

# d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                                                             data = .x,
#                                                                             iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                              lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 10,
                                             lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
                                             upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE,
                                             # include weights here!
                                             modelweights = 1/sd)))


# check they work
d_fit_Mcap$weighted[[1]]

# get predictions using augment
newdata_Mcap <- tibble(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))
d_preds_Mcap <- d_fit_Mcap %>%
  mutate(., preds = map(weighted, augment, newdata = newdata_Mcap)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Mcap) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
  geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Am')+
       #title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
  #ylim(c(-0.25, 3.5))



### Residual re-sampling weighted non-linear regression

#So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
#number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
#is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
#as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM
# fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                data = Mcap.df_mean,
#                                start = coef(d_fit_Mcap$weighted[[1]]),
#                                lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                weights = 1/sd)

fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
                                    data = Mcap.df_mean,
                                    start = coef(d_fit_Mcap$weighted[[1]]),
                                    lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'pawar_2018'),
                                    upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'pawar_2018'),
                                    weights = 1/sd)

#Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
#are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# perform residual bootstrap
boot2_Mcap <- Boot(fit_nlsLM_Mcap, method = 'residual')

# predict over new data
# boot2_preds_Mcap <- boot2_Mcap$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref,e,eh,th, tref = 18))

boot2_preds_Mcap <- boot2_Mcap$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 18))


# calculate bootstrapped confidence intervals
boot2_conf_preds_Mcap <- boot2_preds_Mcap %>%
  group_by(temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs - 95% CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Mcap, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
  geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Am',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
  #ylim(c(-3, 3.5))

p1
  


## Calculating confidence intervals of estimated and calculated parameters
#As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
#parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
#confidence intervals (using confint-MASS) and asymptotic confidence intervals 
#(using nlstools::confint2()).

library(nlstools)

# get parameters of fitted model
param_Mcap <- broom::tidy(fit_nlsLM_Mcap) %>%
  select(param = term, estimate)


# calculate confidence intervals of models
ci1_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'profile')
#> Waiting for profiling to be done...
#> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# profiling method fails
ci2_Mcap <- mutate(ci1_Mcap, method = 'profile',
              conf_lower = NA,
              conf_upper = NA)


# CIs from residual resampling
ci4_Mcap <- confint(boot2_Mcap, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


ci_Am_Mcap <- bind_rows(ci1_Mcap, ci2_Mcap, ci4_Mcap) %>%
  full_join(., param_Mcap, by = "param")
#> Joining with `by = join_by(param)`

ggplot(ci_Am_Mcap, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       )


# We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

extra_params <- calc_params(fit_nlsLM_Mcap) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_Am_case_Mcap <- left_join(ci_extra_params, extra_params)

ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


ci_extra_params_Am_residual_Mcap <- left_join(ci_extra_params, extra_params)


ggplot(ci_extra_params_Am_residual_Mcap, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
       )





#### Pacu CURVE FIT - Pmax
Pacu.df <- df_wide %>%  
  filter(species=="Pocillopora acuta") 

Pacu.df$Temp.Cat_numeric <- as.numeric(as.character(Pacu.df$Temp.Cat))

# Calculate mean and sd values for Am, grouped by Temp.Cat

Pacu.df_mean <- Pacu.df %>%
  select(Temp.Cat, Am, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Am, na.rm = TRUE),ave_rate = mean(Am, na.rm = TRUE),
            groups = 'drop')


# View the resulting dataframe
print(Pacu.df_mean)

Pacu.df_mean$Temp.Cat <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits_Pacu <- nest(Pacu.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          pawar = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = Temp.Cat, r_tref, e, eh, topt, tref = 18),
#                                           data = .x,
#                                           iter = c(4,4,4,4),
#                                           start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018') - 10,
#                                           start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018') + 10,
#                                           lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018'),
#                                           upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018'),
#                                           supp_errors = 'Y',
#                                           convergence_count = FALSE,
#                                           modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pacu.df_mean$Temp.Cat), max(Pacu.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Am',
#        title = 'Am across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model

# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Am',
#        title = 'Am across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


############ Sharpe Schoolfield 1981 model 
#Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

#### we want to run the model taking into account means and standard deviations of the PI parameters
#https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
#https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


Pacu.df_mean <- Pacu.df %>%
  select(Temp.Cat, Am, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Am),ave_rate = mean(Am),
            groups = 'drop')


Pacu.df_mean$temp <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #temp needs to be numeric

Pacu.df_mean <- Pacu.df_mean %>%
  mutate(Temp.Cat = temp) %>%
  select(temp, everything(), -Temp.Cat)

# d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                              lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 10,
                                             lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
                                             upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE,
                                             # include weights here!
                                             modelweights = 1/sd)))


# check they work
d_fit_Pacu$weighted[[1]]


# get predictions using augment
newdata_Pacu <- tibble(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out = 100))
d_preds_Pacu <- d_fit_Pacu %>%
  mutate(., preds = map(weighted, augment, newdata = newdata_Pacu)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pacu) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
  geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Am')+
  #title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 



### Residual re-sampling weighted non-linear regression

#So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
#number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
#is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
#as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM

# fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                data = Pacu.df_mean,
#                                start = coef(d_fit_Pacu$weighted[[1]]),
#                                lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                weights = 1/sd)

fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
                                    data = Pacu.df_mean,
                                    start = coef(d_fit_Pacu$weighted[[1]]),
                                    lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'pawar_2018'),
                                    upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'pawar_2018'),
                                    weights = 1/sd)


#Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
#are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# perform residual bootstrap
boot2_Pacu <- Boot(fit_nlsLM_Pacu, method = 'residual')


# predict over new data

# boot2_preds_Pacu <- boot2_Pacu$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref,e,eh,th, tref = 18))

boot2_preds_Pacu <- boot2_Pacu$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 18))


# calculate bootstrapped confidence intervals
boot2_conf_preds_Pacu <- boot2_preds_Pacu %>%
  group_by(temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs - 95% CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pacu, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
  geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Am',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
#ylim(c(-3, 3.5))

p1


## Calculating confidence intervals of estimated and calculated parameters
#As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
#parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
#confidence intervals (using confint-MASS) and asymptotic confidence intervals 
#(using nlstools::confint2()).

library(nlstools)

# get parameters of fitted model
param_Pacu <- broom::tidy(fit_nlsLM_Pacu) %>%
  select(param = term, estimate)


# calculate confidence intervals of models
ci1_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'profile')
#> Waiting for profiling to be done...

ci2_Pacu <- mutate(ci1_Pacu, method = 'profile',
              conf_lower = NA,
              conf_upper = NA)

# CIs from case resampling
ci3_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')


# CIs from residual resampling
ci4_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


ci_Am_Pacu <- bind_rows(ci1_Pacu, ci2_Pacu, ci4_Pacu) %>%
  full_join(., param_Pacu , by = "param")

#> Joining with `by = join_by(param)`

ggplot(ci_Am_Pacu , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters'
    )


# We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

extra_params <- calc_params(fit_nlsLM_Pacu) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_Am_case_Pacu <- left_join(ci_extra_params, extra_params)

ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params_Am_residual_Pacu <- left_join(ci_extra_params, extra_params)
#> Joining with `by = join_by(param)`

ggplot(ci_extra_params_Am_residual_Pacu, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )



#### Pcomp CURVE FIT - Pmax
Pcom.df <- df_wide %>%  
  filter(species=="Porites compressa") 

Pcom.df$Temp.Cat_numeric <- as.numeric(as.character(Pcom.df$Temp.Cat))

# Calculate mean and sd values for Am, grouped by Temp.Cat

Pcom.df_mean <- Pcom.df %>%
  select(Temp.Cat, Am, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Am, na.rm = TRUE),ave_rate = mean(Am, na.rm = TRUE),
            groups = 'drop')


# View the resulting dataframe
print(Pcom.df_mean)

Pcom.df_mean$Temp.Cat <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Pcom.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          pawar = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = Temp.Cat, r_tref, e, eh, topt, tref = 18),
#                                           data = .x,
#                                           iter = c(4,4,4,4),
#                                           start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018') - 10,
#                                           start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018') + 10,
#                                           lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018'),
#                                           upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'pawar_2018'),
#                                           supp_errors = 'Y',
#                                           convergence_count = FALSE,
#                                           modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pcom.df_mean$Temp.Cat), max(Pcom.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Am',
#        title = 'Am across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model

# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Am',
#        title = 'Am across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


############ Pawar 2018, modified Sharpe Schoolfield 1981 model 
#Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

#### we want to run the model taking into account means and standard deviations of the PI parameters
#https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
#https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


Pcom.df_mean <- Pcom.df %>%
  select(Temp.Cat, Am, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Am),ave_rate = mean(Am),
            groups = 'drop')


Pcom.df_mean$temp <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #temp needs to be numeric

Pcom.df_mean <- Pcom.df_mean %>%
  mutate(Temp.Cat = temp) %>%
  select(temp, everything(), -Temp.Cat)

# d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                              lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 10,
                                             lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
                                             upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE,
                                             # include weights here!
                                             modelweights = 1/sd)))


# check they work
d_fit_Pcom$weighted[[1]]


# get predictions using augment
newdata_Pcom <- tibble(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))
d_preds_Pcom <- d_fit_Pcom %>%
  mutate(., preds = map(weighted, augment, newdata = newdata_Pcom)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pcom) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
  geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Am')+
  #title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 



### Residual re-sampling weighted non-linear regression

#So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
#number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
#is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
#as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM

# fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                data = Pcom.df_mean,
#                                start = coef(d_fit_Pcom$weighted[[1]]),
#                                lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                               weights = 1/sd)

fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
                                    data = Pcom.df_mean,
                                    start = coef(d_fit_Pcom$weighted[[1]]),
                                    lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'pawar_2018'),
                                    upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'pawar_2018'),
                                    weights = 1/sd)

#Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
#are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# perform residual bootstrap
boot2_Pcom <- Boot(fit_nlsLM_Pcom, method = 'residual')


# predict over new data
# boot2_preds_Pcom <- boot2_Pcom$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref,e,eh,th, tref = 18))

boot2_preds_Pcom <- boot2_Pcom$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 18))


# calculate bootstrapped confidence intervals
boot2_conf_preds_Pcom <- boot2_preds_Pcom %>%
  group_by(temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs - 95% CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pcom, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
  geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Am',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
#ylim(c(-3, 3.5))

p1


## Calculating confidence intervals of estimated and calculated parameters
#As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
#parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
#confidence intervals (using confint-MASS) and asymptotic confidence intervals 
#(using nlstools::confint2()).

library(nlstools)

# get parameters of fitted model
param_Pcom <- broom::tidy(fit_nlsLM_Pcom) %>%
  select(param = term, estimate)


# calculate confidence intervals of models
ci1_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'profile')
#> Waiting for profiling to be done...
#> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# profiling method fails
ci2_Pcom <- mutate(ci1_Pcom, method = 'profile',
              conf_lower = NA,
              conf_upper = NA)

# CIs from case resampling
# ci3 <- confint(boot1, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')


# CIs from residual resampling
ci4_Pcom <- confint(boot2_Pcom, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


ci_Am_Pcom <- bind_rows(ci1_Pcom, ci2_Pcom, ci4_Pcom ) %>%
  full_join(., param_Pcom, by = "param")
#> Joining with `by = join_by(param)`

ggplot(ci_Am_Pcom , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters'
  )


# We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

extra_params <- calc_params(fit_nlsLM_Pcom) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_Am_case_Pcom <- left_join(ci_extra_params, extra_params)

ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params_Am_residual_Pcom <- left_join(ci_extra_params, extra_params)
#> Joining with `by = join_by(param)`

ggplot(ci_extra_params_Am_residual_Pcom, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )

#filter to only the most relavent and well characterized parameters 
# All_params <- All_params %>%  
#   filter(!param=="ctmin") %>% 
#   filter(!param=="ctmax") %>% 
#   filter(!param=="skewness") %>% 
#   filter(!param=="thermal_safety_margin") %>% 
#   filter(!param=="thermal_tolerance") %>% 
#   filter(!param=="q10")%>% 
#   filter(!param=="breadth") 



#set plot colors 
cols <- c("Mcap"="green",  
          "Pacu"="cyan", 
          "Pcom"="orange") 

# Ensure Temp.Cat is numeric
df_wide$Temp.Cat <- as.numeric(as.character(df_wide$Temp.Cat))
Mcap.df$Temp.Cat <- as.numeric(as.character(Mcap.df$Temp.Cat))
Pacu.df$Temp.Cat <- as.numeric(as.character(Pacu.df$Temp.Cat))
Pcom.df$Temp.Cat <- as.numeric(as.character(Pcom.df$Temp.Cat))

# plot data and model fit for the 3 species
TPC.plot_Pmax <- ggplot(data=df_wide, aes()) + 
  geom_point(aes(Temp.Cat,Am, color = "Mcap"), data =  Mcap.df, size = 2, alpha = 0.5) + 
  geom_point(aes(Temp.Cat,Am,  color = "Pacu"), data = Pacu.df, size = 2, alpha = 0.5) + 
  geom_point(aes(Temp.Cat,Am, color = "Pcom"), data = Pcom.df, size = 2, alpha = 0.5) + 
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = "green", alpha = 0.3) + 
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'cyan', alpha = 0.3) + 
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'orange', alpha = 0.3) + 
  geom_line(aes(temp, .fitted), data = d_preds_Mcap, color = "green") +
  geom_line(aes(temp, .fitted), data = d_preds_Pacu, color = "cyan") +
  geom_line(aes(temp, .fitted), data = d_preds_Pcom, color = "orange") +
  #xlim(11.5,40.5)+ 
  scale_x_continuous(breaks=c(12,14,16,18,20,22,24,26,28,30,32,34,36))+
  theme_bw(base_size = 12) + 
  scale_colour_manual(name="Species",values=cols)+ 
  scale_fill_manual(name = "Species", values = cols) + 
  theme(legend.position = "top", 
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  labs(x = 'Temperature (ºC)', 
       y = expression("Pmax µmol O2" ~cm^{-2}~h^{-1}))

TPC.plot_Pmax  



#### Plot TPC parameters
# Combine the data frames and add a source column
ci_extra_params_Am_case_Mcap <- ci_extra_params_Am_case_Mcap %>% mutate(source = "Mcap")
ci_extra_params_Am_case_Pacu <- ci_extra_params_Am_case_Pacu %>% mutate(source = "Pacu")
ci_extra_params_Am_case_Pcom <- ci_extra_params_Am_case_Pcom %>% mutate(source = "Pcom")

ci_extra_params_Am_residual_Mcap <- ci_extra_params_Am_residual_Mcap %>% mutate(source = "Mcap")
ci_extra_params_Am_residual_Pacu <- ci_extra_params_Am_residual_Pacu %>% mutate(source = "Pacu")
ci_extra_params_Am_residual_Pcom <- ci_extra_params_Am_residual_Pcom %>% mutate(source = "Pcom")

# Combine all data frames into one
ci_extra_params_Am_case_combined <- bind_rows(ci_extra_params_Am_case_Mcap , ci_extra_params_Am_case_Pacu, ci_extra_params_Am_case_Pcom )

ci_extra_params_Am_residual_combined <- bind_rows(ci_extra_params_Am_residual_Mcap , ci_extra_params_Am_residual_Pacu, ci_extra_params_Am_residual_Pcom )

#set plot colors 
cols <- c("Mcap"="green",  
          "Pacu"="cyan", 
          "Pcom"="orange") 

# Plot the combined data
ggplot(ci_extra_params_Am_residual_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  #facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

# Plot the combined data
ggplot(ci_extra_params_Am_case_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  #facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Add a source column to each data frame
ci_Am_Mcap <- ci_Am_Mcap %>% mutate(source = "Mcap")
ci_Am_Pacu <- ci_Am_Pacu %>% mutate(source = "Pacu")
ci_Am_Pcom <- ci_Am_Pcom %>% mutate(source = "Pcom")

# Combine all data frames into one
ci_combined_Am <- bind_rows(ci_Am_Mcap, ci_Am_Pacu, ci_Am_Pcom )

#set plot colors 
cols <- c("Mcap"="green",  
          "Pacu"="cyan", 
          "Pcom"="orange") 

# Plot the combined data
ggplot(ci_combined_Am, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(ci_combined_Am, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(param ~ method, scales = 'free') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ci_residual_Am_combined <- ci_combined_Am %>%
  filter(method == "residual bootstrap")

ci_residual_Am_combined$param = as.factor(ci_residual_Am_combined$param)

ggplot(ci_residual_Am_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap("param", scales = 'free') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Confidence intervals for model parameters - residual bootstrap',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


######## Statistical analysis of TPC parameters
### we have already calculated the confidence intervals for the parameters, so we can use these to calculate the
### statistical differences between the species.

# Filter rows with "residual bootstrap" in the method column
ci_residual_Am_combined <- ci_combined_Am %>%
  filter(method == "residual bootstrap")

# View the new table
print(ci_residual_Am_combined)

#Formally test whether confidence intervals overlap using a statistical approach called confidence interval comparison. One common method 
#is to calculate the Z-score for the difference between two estimates and determine whether the difference is statistically significant.
#The value 3.92 below is used to approximate the standard error of the confidence interval (CI) bounds. It comes from the fact that a 95% 
#confidence interval corresponds to approximately 1.96 standard deviations above and below the mean in a normal distribution. Since the 
#CI width spans both directions (upper and lower bounds), the total width is 2 * 1.96 = 3.92.

# Load required libraries
library(dplyr)

# Function to test overlap of confidence intervals
test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
                            estimate2, conf_lower2, conf_upper2) {
  # Calculate the Z-score for the difference
  diff <- abs(estimate1 - estimate2)
  pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
                      ((conf_upper2 - conf_lower2) / 3.92)^2)
  z_score <- diff / pooled_se
  
  # Calculate p-value
  p_value <- 2 * (1 - pnorm(z_score))
  
  return(p_value)
}

# Apply the function to compare parameters between species
ci_Am_comparison <- ci_residual_Am_combined %>%
  group_by(param) %>%
  summarise(
    p_value_Mcap_vs_Pacu = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
    ),
    p_value_Mcap_vs_Pcom = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    ),
    p_value_Pacu_vs_Pcom = test_ci_overlap(
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    )
  )

# View the results
print(ci_Am_comparison)

#The 95 % confidence intervals from two subgroups or studies may overlap substantially and yet the test 
#for difference between them may still produce P < 0.05. 
#https://pmc.ncbi.nlm.nih.gov/articles/PMC4877414/
# Confidence Interval Misinterpretation: Confidence intervals are often misinterpreted as a direct test 
# of significance. However, overlapping confidence intervals do not necessarily mean the difference is not 
# significant. This is because the overlap does not account for the combined uncertainty of both estimates.

# Calculate max conf_upper for each parameter
ci_residual_Am_combined <- ci_residual_Am_combined %>%
  group_by(param) %>%
  mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
  ungroup()

# Combine ci_residual_combined with ci_comparison
ci_plot_data_Am <- ci_residual_Am_combined %>%
  left_join(ci_Am_comparison, by = "param")

# Prepare the statistical annotations for ci_plot_data
# Compute y.position values for each parameter and source
y_positions <- ci_plot_data_Am %>%
  group_by(param) %>%
  summarise(
    max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
    .groups = "drop"
  ) %>%
  mutate(
    y_base = max_value  # Add a base offset to the maximum value
  )

# Add y.position to annotations with dynamic spacing
ci_plot_data_Am <- ci_plot_data_Am %>%
  left_join(y_positions, by = "param") %>%
  group_by(param) %>%
  mutate(
    y.position = y_base + (row_number() - 1) * (y_base * 0.5)  # Dynamic offset based on y_base
  ) %>%
  ungroup() %>%
  select(-max_value, -y_base)  # Remove intermediate columns

# # Filter for significant comparisons only
# ci_plot_data_Am <- ci_plot_data_Am %>%
#   filter(p_value_Mcap_vs_Pacu < 0.05 | p_value_Mcap_vs_Pcom < 0.05 | p_value_Pacu_vs_Pcom < 0.05)

# Add significance levels for each comparison
ci_plot_data_Am <- ci_plot_data_Am %>%
  mutate(
    signif_label_Mcap_vs_Pacu = case_when(
      p_value_Mcap_vs_Pacu < 0.001 ~ "***",
      p_value_Mcap_vs_Pacu < 0.01 ~ "**",
      p_value_Mcap_vs_Pacu < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Mcap_vs_Pcom = case_when(
      p_value_Mcap_vs_Pcom < 0.001 ~ "***",
      p_value_Mcap_vs_Pcom < 0.01 ~ "**",
      p_value_Mcap_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Pacu_vs_Pcom = case_when(
      p_value_Pacu_vs_Pcom < 0.001 ~ "***",
      p_value_Pacu_vs_Pcom < 0.01 ~ "**",
      p_value_Pacu_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )


# Compute y.position dynamically for each comparison
ci_plot_data_Am <- ci_plot_data_Am %>%
  group_by(param) %>%
  mutate(
    y_max = max(conf_upper, na.rm = TRUE),
    y_min = min(conf_lower, na.rm = TRUE),
    y_range = y_max - y_min,
    y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
    y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
    y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
  ) %>%
  ungroup()


# Create the plot
ci_data_Am_plot <- ggplot(ci_plot_data_Am, aes(x = source, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  # Add bars and asterisks for significant comparisons
  geom_segment(
    aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
    data = ci_plot_data_Am %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 1.5, y = y_position_Mcap_vs_Pacu + 0.000005, label = signif_label_Mcap_vs_Pacu),
    data = ci_plot_data_Am %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
    data = ci_plot_data_Am %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2, y = y_position_Mcap_vs_Pcom + 0.000005, label = signif_label_Mcap_vs_Pcom),
    data = ci_plot_data_Am %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
    data = ci_plot_data_Am %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.000005, label = signif_label_Pacu_vs_Pcom),
    data = ci_plot_data_Am %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  theme_bw() +
  facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(
    title = 'Confidence intervals for model parameters with significance',
    y = 'Estimate',
    x = 'Source'
  ) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) 

# # Update the plot
# ggplot(ci_plot_data_Am, aes(x = source, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   geom_segment(data = ci_plot_data, 
#                aes(x = 1, xend = 2, 
#                    y = max_conf_upper + 0.1, 
#                    yend = max_conf_upper + 0.1), 
#                inherit.aes = FALSE, col = "black") +
#   geom_text(data = ci_plot_data, 
#             aes(x = 1.5, 
#                 y = max_conf_upper + 0.15, 
#                 label = ifelse(p_value_Mcap_vs_Pacu < 0.05, "*", "")), 
#             inherit.aes = FALSE, col = "black") +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Confidence intervals for model parameters with significance',
#        y = 'Estimate',
#        x = 'Source') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"))



#### Function to test overlap of confidence intervals - EXTRA PARAMETERS
#ci_extra_params_residual_combined

test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
                            estimate2, conf_lower2, conf_upper2) {
  # Calculate the Z-score for the difference
  diff <- abs(estimate1 - estimate2)
  pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
                      ((conf_upper2 - conf_lower2) / 3.92)^2)
  z_score <- diff / pooled_se
  
  # Calculate p-value
  p_value <- 2 * (1 - pnorm(z_score))
  
  return(p_value)
}

# Apply the function to compare parameters between species
ci_comparison_Am_extra_params <- ci_extra_params_Am_residual_combined %>%
  group_by(param) %>%
  summarise(
    p_value_Mcap_vs_Pacu = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
    ),
    p_value_Mcap_vs_Pcom = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    ),
    p_value_Pacu_vs_Pcom = test_ci_overlap(
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    )
  )

# View the results
print(ci_comparison_Am_extra_params)


# Combine ci_extra_params_case_combined with ci_comparison_extra_params
ci_plot_data_Am_extra <- ci_extra_params_Am_residual_combined %>%
  left_join(ci_comparison_Am_extra_params, by = "param")

# # Calculate max_conf_upper for each parameter
# ci_plot_data_Am_extra <- ci_plot_data_Am_extra %>%
#   group_by(param) %>%
#   mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
#   ungroup()

# Compute y.position values for each parameter and source
y_positions <- ci_plot_data_Am_extra %>%
  group_by(param) %>%
  summarise(
    max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
    .groups = "drop"
  ) %>%
  mutate(
    y_base = max_value  # Add a base offset to the maximum value
  )

# Add y.position to annotations with dynamic spacing
ci_plot_data_Am_extra <- ci_plot_data_Am_extra %>%
  left_join(y_positions, by = "param") %>%
  group_by(param) %>%
  mutate(
    y.position = y_base + (row_number() - 1) * (y_base * 0.1)  # Dynamic offset based on y_base
  ) %>%
  ungroup() %>%
  select(-max_value, -y_base)  # Remove intermediate columns

# Filter for significant comparisons only
# ci_plot_data_Am_extra <- ci_plot_data_Am_extra %>%
#   filter(p_value_Mcap_vs_Pacu < 0.05 | p_value_Mcap_vs_Pcom < 0.05 | p_value_Pacu_vs_Pcom < 0.05)


# Add significance levels for each comparison
ci_plot_data_Am_extra <- ci_plot_data_Am_extra %>%
  mutate(
    signif_label_Mcap_vs_Pacu = case_when(
      p_value_Mcap_vs_Pacu < 0.001 ~ "***",
      p_value_Mcap_vs_Pacu < 0.01 ~ "**",
      p_value_Mcap_vs_Pacu < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Mcap_vs_Pcom = case_when(
      p_value_Mcap_vs_Pcom < 0.001 ~ "***",
      p_value_Mcap_vs_Pcom < 0.01 ~ "**",
      p_value_Mcap_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Pacu_vs_Pcom = case_when(
      p_value_Pacu_vs_Pcom < 0.001 ~ "***",
      p_value_Pacu_vs_Pcom < 0.01 ~ "**",
      p_value_Pacu_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )

# Compute y.position dynamically for each comparison
ci_plot_data_Am_extra <- ci_plot_data_Am_extra %>%
  group_by(param) %>%
  mutate(
    y_max = max(conf_upper, na.rm = TRUE),
    y_min = min(conf_lower, na.rm = TRUE),
    y_range = y_max - y_min,
    y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
    y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
    y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
  ) %>%
  ungroup()

# Create the plot
ci_data_Am_extra_plot <- ggplot(ci_plot_data_Am_extra, aes(x = source, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  # Add bars and asterisks for significant comparisons
  geom_segment(
    aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
    data = ci_plot_data_Am_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 1.5, y = y_position_Mcap_vs_Pacu + 0.000005, label = signif_label_Mcap_vs_Pacu),
    data = ci_plot_data_Am_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
    data = ci_plot_data_Am_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2, y = y_position_Mcap_vs_Pcom + 0.000005, label = signif_label_Mcap_vs_Pcom),
    data = ci_plot_data_Am_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
    data = ci_plot_data_Am_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.000005, label = signif_label_Pacu_vs_Pcom),
    data = ci_plot_data_Am_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  theme_bw() +
  facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(
    title = 'Confidence intervals for extra parameters with significance',
    y = 'Estimate',
    x = 'Source'
  ) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )








############ Mcap CURVE FIT - Rd (respiration rate)
Mcap.df <- df_wide %>%  
  filter(species=="Montipora capitata") 

Mcap.df$Temp.Cat_numeric <- as.numeric(as.character(Mcap.df$Temp.Cat))

# Calculate mean and sd values for Rd, grouped by Temp.Cat

# Mcap.df_mean <- Mcap.df %>%
#   group_by(Temp.Cat) %>%
#   summarise(across(7:11, mean, na.rm = TRUE)) %>%
#   mutate(colony_id = "MCAP") %>%
#   select(Temp.Cat, Rd, colony_id)

Mcap.df_mean <- Mcap.df %>%
  dplyr::select(Temp.Cat, Rd, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Rd, na.rm = TRUE), ave_rate = mean(Rd, na.rm = TRUE),
            groups = 'drop')

# View the resulting dataframe
print(Mcap.df_mean)

Mcap.df_mean$Temp.Cat <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Mcap.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          #                                    # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Mcap.df_mean$Temp.Cat), max(Mcap.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   dplyr::select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   dplyr::select(-fit) %>%
#   unnest(info) %>%
#   dplyr::select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


############ Sharpe Schoolfield 1981 model 
#Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

#### we want to run the model taking into account means and standard deviations of the PI parameters
#https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
#https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


Mcap.df_mean <- Mcap.df %>%
  dplyr::select(Temp.Cat, Rd, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Rd),ave_rate = mean(Rd),
            groups = 'drop')

Mcap.df_mean$temp <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #temp needs to be numeric

Mcap.df_mean <- Mcap.df_mean %>%
  mutate(Temp.Cat = temp) %>%
  select(temp, everything(), -Temp.Cat)

# d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 1,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 1,
#                                              lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 10,
                                           #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
                                           #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE,
                                           # include weights here!
                                           modelweights = 1/sd)))

                                           
# check they work
d_fit_Mcap$weighted[[1]]

# get predictions using augment
newdata_Mcap <- tibble(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))
d_preds_Mcap <- d_fit_Mcap %>%
  mutate(., preds = map(weighted, augment, newdata = newdata_Mcap)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Mcap) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
  geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Rd')+
  #title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
#ylim(c(-0.25, 3.5))



### Residual re-sampling weighted non-linear regression

#So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
#number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
#is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
#as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM
# fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                     data = Mcap.df_mean,
#                                     start = coef(d_fit_Mcap$weighted[[1]]),
#                                     lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                     upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                     weights = 1/sd)

fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
                                    data = Mcap.df_mean,
                                    start = coef(d_fit_Mcap$weighted[[1]]),
                                    #lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
                                    #upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
                                    weights = 1/sd)

#Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
#are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# perform residual bootstrap
boot2_Mcap <- Boot(fit_nlsLM_Mcap, method = 'residual')

# predict over new data
# boot2_preds_Mcap <- boot2_Mcap$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref,e,eh,th, tref = 18))

boot2_preds_Mcap <- boot2_Mcap$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp, r_tref,e,el,tl,tref = 25))

# calculate bootstrapped confidence intervals
boot2_conf_preds_Mcap <- boot2_preds_Mcap %>%
  group_by(temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs - 95% CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Mcap, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
  geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Rd',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
#ylim(c(-3, 3.5))

p1



## Calculating confidence intervals of estimated and calculated parameters
#As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
#parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
#confidence intervals (using confint-MASS) and asymptotic confidence intervals 
#(using nlstools::confint2()).

library(nlstools)

# get parameters of fitted model
param_Mcap <- broom::tidy(fit_nlsLM_Mcap) %>%
  select(param = term, estimate)


# calculate confidence intervals of models
ci1_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'profile')
#> Waiting for profiling to be done...
#> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# profiling method fails
ci2_Mcap <- mutate(ci1_Mcap, method = 'profile',
                   conf_lower = NA,
                   conf_upper = NA)


# CIs from residual resampling
ci4_Mcap <- confint(boot2_Mcap, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')
## [1] "All values of t are equal to  0.165850014721406 \n Cannot calculate confidence intervals"

ci_Rd_Mcap <- bind_rows(ci1_Mcap, ci2_Mcap, ci4_Mcap) %>%
  full_join(., param_Mcap, by = "param")
#> Joining with `by = join_by(param)`

ggplot(ci_Rd_Mcap, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'For the chlorella TPC; profile method failes')

#### asymptotic works


# We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

extra_params <- calc_params(fit_nlsLM_Mcap) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 400, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_Rd_case_Mcap <- left_join(ci_extra_params, extra_params)

ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 400, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params_Rd_residual_Mcap <- left_join(ci_extra_params, extra_params)


ggplot(ci_extra_params_Rd_case_Mcap, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )

ggplot(ci_extra_params_Rd_residual_Mcap, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )




#### Pacu CURVE FIT - Rd
Pacu.df <- df_wide %>%  
  filter(species=="Pocillopora acuta") 

Pacu.df$Temp.Cat_numeric <- as.numeric(as.character(Pacu.df$Temp.Cat))

# Calculate mean and sd values for Rd, grouped by Temp.Cat

# Mcap.df_mean <- Mcap.df %>%
#   group_by(Temp.Cat) %>%
#   summarise(across(7:11, mean, na.rm = TRUE)) %>%
#   mutate(colony_id = "MCAP") %>%
#   select(Temp.Cat, Rd, colony_id)

Pacu.df_mean <- Pacu.df %>%
  select(Temp.Cat, Rd, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Rd,na.rm = TRUE),ave_rate = mean(Rd,na.rm = TRUE),
            groups = 'drop')


# View the resulting dataframe
print(Pacu.df_mean)

Pacu.df_mean$Temp.Cat <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits_Pacu <- nest(Pacu.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pacu.df_mean$Temp.Cat), max(Pacu.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


############ Sharpe Schoolfield 1981 model 
#Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

#### we want to run the model taking into account means and standard deviations of the PI parameters
#https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
#https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


Pacu.df_mean <- Pacu.df %>%
  select(Temp.Cat, Rd, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Rd),ave_rate = mean(Rd),
            groups = 'drop')


Pacu.df_mean$temp <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #temp needs to be numeric

Pacu.df_mean <- Pacu.df_mean %>%
  mutate(Temp.Cat = temp) %>%
  select(temp, everything(), -Temp.Cat)

# d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                              lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp,r_tref,e,el,tl, tref = 25),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 1,
                                             start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 1,
                                             #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
                                             #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE,
                                             # include weights here!
                                             modelweights = 1/sd)))

# check they work
d_fit_Pacu$weighted[[1]]


# get predictions using augment
newdata_Pacu <- tibble(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out = 100))
d_preds_Pacu <- d_fit_Pacu %>%
  mutate(., preds = map(weighted, augment, newdata = newdata_Pacu)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pacu) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
  geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Rd')+
  #title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 



### Residual re-sampling weighted non-linear regression

#So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
#number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
#is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
#as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM
# fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                     data = Pacu.df_mean,
#                                     start = coef(d_fit_Pacu$weighted[[1]]),
#                                     lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                     upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                     weights = 1/sd)

fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
                                    data = Pacu.df_mean,
                                    start = coef(d_fit_Pacu$weighted[[1]]),
                                    #lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
                                    #upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
                                    weights = 1/sd)

#Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
#are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# perform residual bootstrap
boot2_Pacu <- Boot(fit_nlsLM_Pacu, method = 'residual')

# predict over new data
# boot2_preds_Pacu <- boot2_Pacu$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref,e,eh,th, tref = 18))

boot2_preds_Pacu <- boot2_Pacu$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp, r_tref,e,el,tl, tref = 25))


# calculate bootstrapped confidence intervals
boot2_conf_preds_Pacu <- boot2_preds_Pacu %>%
  group_by(temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs - 95% CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pacu, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
  geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Rd',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
#ylim(c(-3, 3.5))

p1


## Calculating confidence intervals of estimated and calculated parameters
#As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
#parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
#confidence intervals (using confint-MASS) and asymptotic confidence intervals 
#(using nlstools::confint2()).

library(nlstools)

# get parameters of fitted model
param_Pacu <- broom::tidy(fit_nlsLM_Pacu) %>%
  select(param = term, estimate)


# calculate confidence intervals of models
ci1_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'profile')
#> Waiting for profiling to be done...

ci2_Pacu <- mutate(ci1_Pacu, method = 'profile',
                   conf_lower = NA,
                   conf_upper = NA)

# CIs from case resampling
ci3_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')


# CIs from residual resampling
ci4_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')
# [1] "All values of t are equal to  0.210223059177217 \n Cannot calculate confidence intervals"

ci_Rd_Pacu <- bind_rows(ci1_Pacu,ci2_Pacu,ci3_Pacu, ci4_Pacu) %>%
  full_join(., param_Pacu , by = "param")
ci_Rd_Pacu <- bind_rows(ci2_Pacu,ci3_Pacu, ci4_Pacu) %>%
  full_join(., param_Pacu , by = "param")
#> Joining with `by = join_by(param)`

ggplot(ci_Rd_Pacu , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters'
  )


# We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

extra_params <- calc_params(fit_nlsLM_Pacu) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_Rd_case_Pacu <- left_join(ci_extra_params, extra_params)

ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params_Rd_residual_Pacu <- left_join(ci_extra_params, extra_params)


ggplot(ci_extra_params_Rd_residual_Pacu, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )

ggplot(ci_extra_params_Rd_case_Pacu, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )



#### Pcomp CURVE FIT - Rd
Pcom.df <- df_wide %>%  
  filter(species=="Porites compressa") 

Pcom.df$Temp.Cat_numeric <- as.numeric(as.character(Pcom.df$Temp.Cat))

# Calculate mean and sd values for Rd, grouped by Temp.Cat

# Mcap.df_mean <- Mcap.df %>%
#   group_by(Temp.Cat) %>%
#   summarise(across(7:11, mean, na.rm = TRUE)) %>%
#   mutate(colony_id = "MCAP") %>%
#   select(Temp.Cat, Rd, colony_id)

Pcom.df_mean <- Pcom.df %>%
  select(Temp.Cat, Rd, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Rd, na.rm = TRUE),ave_rate = mean(Rd, na.rm = TRUE),
            groups = 'drop')


# View the resulting dataframe
print(Pcom.df_mean)

Pcom.df_mean$Temp.Cat <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Pcom.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pcom.df_mean$Temp.Cat), max(Pcom.df_mean$Temp.Cat), length.out = 10))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


############ Sharpe Schoolfield 1981 model 
#Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

#### we want to run the model taking into account means and standard deviations of the PI parameters
#https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
#https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


Pcom.df_mean <- Pcom.df %>%
  select(Temp.Cat, Rd, colony_id) %>%
  group_by(Temp.Cat) %>%
  summarise(., sd = sd(Rd),ave_rate = mean(Rd),
            groups = 'drop')


Pcom.df_mean$temp <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #temp needs to be numeric

Pcom.df_mean <- Pcom.df_mean %>%
  mutate(Temp.Cat = temp) %>%
  select(temp, everything(), -Temp.Cat)

# d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                              lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 10,
                                             #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
                                             #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE,
                                             # include weights here!
                                             modelweights = 1/sd)))

# check they work
d_fit_Pcom$weighted[[1]]



# get predictions using augment
newdata_Pcom <- tibble(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))
d_preds_Pcom <- d_fit_Pcom %>%
  mutate(., preds = map(weighted, augment, newdata = newdata_Pcom)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pcom) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
  geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Rd')+
  #title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 



### Residual re-sampling weighted non-linear regression

#So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
#number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
#is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
#as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM
# fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 26.8),
#                                     data = Pcom.df_mean,
#                                     start = coef(d_fit_Pcom$weighted[[1]]),
#                                     lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                     upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                     weights = 1/sd)

# refit model using nlsLM
fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp,r_tref,e,el,tl, tref = 25),
                                    data = Pcom.df_mean,
                                    start = coef(d_fit_Pcom$weighted[[1]]),
                                    #lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
                                    #upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
                                    weights = 1/sd)

#Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
#are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# perform residual bootstrap
boot2_Pcom <- Boot(fit_nlsLM_Pcom, method = 'residual')


# predict over new data
# boot2_preds_Pcom <- boot2_Pcom$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref,e,eh,th, tref = 18))

boot2_preds_Pcom <- boot2_Pcom$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp, r_tref,e,el,tl, tref = 25))

# calculate bootstrapped confidence intervals
boot2_conf_preds_Pcom <- boot2_preds_Pcom %>%
  group_by(temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs - 95% CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds_Pcom, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
  geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Rd',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) 
#ylim(c(-3, 3.5))

p1


## Calculating confidence intervals of estimated and calculated parameters
#As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
#parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
#confidence intervals (using confint-MASS) and asymptotic confidence intervals 
#(using nlstools::confint2()).

library(nlstools)

# get parameters of fitted model
param_Pcom <- broom::tidy(fit_nlsLM_Pcom) %>%
  select(param = term, estimate)


# calculate confidence intervals of models
ci1_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'profile')
#> Waiting for profiling to be done...
#> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# profiling method fails
ci2_Pcom <- mutate(ci1_Pcom, method = 'profile',
                   conf_lower = NA,
                   conf_upper = NA)

# CIs from case resampling
ci3 <- confint(boot2_Pcom, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')


# CIs from residual resampling
ci4_Pcom <- confint(boot2_Pcom, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')
# [1] "All values of t are equal to  0.267937013495708 \n Cannot calculate confidence intervals"


ci_Rd_Pcom <- bind_rows(ci1_Pcom, ci2_Pcom, ci4_Pcom  ) %>%
  full_join(., param_Pcom, by = "param")

#> Joining with `by = join_by(param)`

ggplot(ci_Rd_Pcom , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters'
  )


# We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

extra_params <- calc_params(fit_nlsLM_Pcom) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_Rd_case_Pcom<- left_join(ci_extra_params, extra_params)

ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 100, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params_Rd_residual_Pcom<- left_join(ci_extra_params, extra_params)


ggplot(ci_extra_params_Rd_residual_Pcom, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters'
  )

#filter to only the most relavent and well characterized parameters 
# All_params <- All_params %>%  
#   filter(!param=="ctmin") %>% 
#   filter(!param=="ctmax") %>% 
#   filter(!param=="skewness") %>% 
#   filter(!param=="thermal_safety_margin") %>% 
#   filter(!param=="thermal_tolerance") %>% 
#   filter(!param=="q10")%>% 
#   filter(!param=="breadth") 



#set plot colors 
cols <- c("Mcap"="green",  
          "Pacu"="cyan", 
          "Pcom"="orange") 

# Ensure Temp.Cat is numeric
df_wide$Temp.Cat <- as.numeric(as.character(df_wide$Temp.Cat))
Mcap.df$Temp.Cat <- as.numeric(as.character(Mcap.df$Temp.Cat))
Pacu.df$Temp.Cat <- as.numeric(as.character(Pacu.df$Temp.Cat))
Pcom.df$Temp.Cat <- as.numeric(as.character(Pcom.df$Temp.Cat))

# plot data and model fit for the 3 species
TPC.plot_Rd <- ggplot(data=df_wide, aes()) + 
  geom_point(aes(Temp.Cat,Rd, color = "Mcap"), data =  Mcap.df, size = 2, alpha = 0.5) + 
  geom_point(aes(Temp.Cat,Rd,  color = "Pacu"), data = Pacu.df, size = 2, alpha = 0.5) + 
  geom_point(aes(Temp.Cat,Rd, color = "Pcom"), data = Pcom.df, size = 2, alpha = 0.5) + 
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = "green", alpha = 0.3) + 
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'cyan', alpha = 0.3) + 
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'orange', alpha = 0.3) + 
  geom_line(aes(temp, .fitted), data = d_preds_Mcap, color = "green") +
  geom_line(aes(temp, .fitted), data = d_preds_Pacu, color = "cyan") +
  geom_line(aes(temp, .fitted), data = d_preds_Pcom, color = "orange") +
  xlim(11.5,40.5)+ 
  scale_x_continuous(breaks=c(12,14,16,18,20,22,24,26,28,30,32,34,36))+
  theme_bw(base_size = 12) + 
  scale_colour_manual(name="Species",values=cols)+ 
  scale_fill_manual(name = "Species", values = cols) + 
  theme(legend.position = "top", 
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  labs(x = 'Temperature (ºC)', 
       y = expression("Rd µmol O2" ~cm^{-2}~h^{-1}))

TPC.plot_Rd  



#### Plot TPC parameters
# Combine the data frames and add a source column
ci_extra_params_Rd_case_Mcap <- ci_extra_params_Rd_case_Mcap %>% mutate(source = "Mcap")
ci_extra_params_Rd_case_Pacu <- ci_extra_params_Rd_case_Pacu %>% mutate(source = "Pacu")
ci_extra_params_Rd_case_Pcom <- ci_extra_params_Rd_case_Pcom %>% mutate(source = "Pcom")

ci_extra_params_Rd_residual_Mcap <- ci_extra_params_Rd_residual_Mcap %>% mutate(source = "Mcap")
ci_extra_params_Rd_residual_Pacu <- ci_extra_params_Rd_residual_Pacu %>% mutate(source = "Pacu")
ci_extra_params_Rd_residual_Pcom <- ci_extra_params_Rd_residual_Pcom %>% mutate(source = "Pcom")

# Combine all data frames into one
ci_extra_params_Rd_case_combined <- bind_rows(ci_extra_params_Rd_case_Mcap , ci_extra_params_Rd_case_Pacu, ci_extra_params_Rd_case_Pcom )

ci_extra_params_Rd_residual_combined <- bind_rows(ci_extra_params_Rd_residual_Mcap , ci_extra_params_Rd_residual_Pacu, ci_extra_params_Rd_residual_Pcom )

#set plot colors 
cols <- c("Mcap"="green",  
          "Pacu"="cyan", 
          "Pcom"="orange") 

# Plot the combined data
ggplot(ci_extra_params_Rd_residual_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  #facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

# Plot the combined data
ggplot(ci_extra_params_Rd_case_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  #facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



# Add a source column to each data frame
ci_Rd_Mcap <- ci_Rd_Mcap %>% mutate(source = "Mcap")
ci_Rd_Pacu <- ci_Rd_Pacu %>% mutate(source = "Pacu")
ci_Rd_Pcom <- ci_Rd_Pcom %>% mutate(source = "Pcom")

# Combine all data frames into one
ci_Rd_combined <- bind_rows(ci_Rd_Mcap, ci_Rd_Pacu, ci_Rd_Pcom )

#set plot colors 
cols <- c("Mcap"="green",  
          "Pacu"="cyan", 
          "Pcom"="orange") 

# Plot the combined data
ggplot(ci_Rd_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(ci_Rd_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(param ~ method, scales = 'free') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Filter rows with "residual bootstrap" in the method column
ci_residual_Rd_combined <- ci_Rd_combined %>%
  filter(method == "residual bootstrap")

ci_residual_Rd_combined$param = as.factor(ci_residual_Rd_combined$param)

ggplot(ci_residual_Rd_combined, aes(x = param, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap("param", scales = 'free') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(title = 'Confidence intervals for model parameters - residual bootstrap',
       y = 'Estimate',
       x = 'Parameter') +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




######## Statistical analysis of TPC parameters
### we have already calculated the confidence intervals for the parameters, so we can use these to calculate the
### statistical differences between the species.

# Filter rows with "residual bootstrap" in the method column
ci_residual_Rd_combined <- ci_Rd_combined %>%
  filter(method == "residual bootstrap")

# View the new table
print(ci_residual_Rd_combined)

#Formally test whether confidence intervals overlap using a statistical approach called confidence interval comparison. One common method 
#is to calculate the Z-score for the difference between two estimates and determine whether the difference is statistically significant.
#The value 3.92 below is used to approximate the standard error of the confidence interval (CI) bounds. It comes from the fact that a 95% 
#confidence interval corresponds to approximately 1.96 standard deviations above and below the mean in a normal distribution. Since the 
#CI width spans both directions (upper and lower bounds), the total width is 2 * 1.96 = 3.92.

# Load required libraries
library(dplyr)

# Function to test overlap of confidence intervals
test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
                            estimate2, conf_lower2, conf_upper2) {
  # Calculate the Z-score for the difference
  diff <- abs(estimate1 - estimate2)
  pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
                      ((conf_upper2 - conf_lower2) / 3.92)^2)
  z_score <- diff / pooled_se
  
  # Calculate p-value
  p_value <- 2 * (1 - pnorm(z_score))
  
  return(p_value)
}

# Apply the function to compare parameters between species
ci_Rd_comparison <- ci_residual_Rd_combined %>%
  group_by(param) %>%
  summarise(
    p_value_Mcap_vs_Pacu = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
    ),
    p_value_Mcap_vs_Pcom = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    ),
    p_value_Pacu_vs_Pcom = test_ci_overlap(
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    )
  )

# View the results
print(ci_Rd_comparison)


#The 95 % confidence intervals from two subgroups or studies may overlap substantially and yet the test 
#for difference between them may still produce P < 0.05. 
#https://pmc.ncbi.nlm.nih.gov/articles/PMC4877414/
# Confidence Interval Misinterpretation: Confidence intervals are often misinterpreted as a direct test 
# of significance. However, overlapping confidence intervals do not necessarily mean the difference is not 
# significant. This is because the overlap does not account for the combined uncertainty of both estimates.

# Calculate max conf_upper for each parameter
ci_residual_Rd_combined <- ci_residual_Rd_combined %>%
  group_by(param) %>%
  mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
  ungroup()

# Combine ci_residual_combined with ci_comparison
ci_plot_Rd_data <- ci_residual_Rd_combined %>%
  left_join(ci_Rd_comparison, by = "param")

# Prepare the statistical annotations for ci_plot_data
# Compute y.position values for each parameter and source
y_positions <- ci_plot_Rd_data %>%
  group_by(param) %>%
  summarise(
    max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
    .groups = "drop"
  ) %>%
  mutate(
    y_base = max_value  # Add a base offset to the maximum value
  )

# Add y.position to annotations with dynamic spacing
ci_plot_Rd_data <- ci_plot_Rd_data %>%
  left_join(y_positions, by = "param") %>%
  group_by(param) %>%
  mutate(
    y.position = y_base + (row_number() - 1) * (y_base * 0.2)  # Dynamic offset based on y_base
  ) %>%
  ungroup() %>%
  select(-max_value, -y_base)  # Remove intermediate columns

# Filter for significant comparisons only
# ci_plot_data <- ci_plot_data %>%
#   filter(p_value_Mcap_vs_Pacu < 0.05 | p_value_Mcap_vs_Pcom < 0.05 | p_value_Pacu_vs_Pcom < 0.05)

# Add significance levels for each comparison
ci_plot_Rd_data <- ci_plot_Rd_data %>%
  mutate(
    signif_label_Mcap_vs_Pacu = case_when(
      p_value_Mcap_vs_Pacu < 0.001 ~ "***",
      p_value_Mcap_vs_Pacu < 0.01 ~ "**",
      p_value_Mcap_vs_Pacu < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Mcap_vs_Pcom = case_when(
      p_value_Mcap_vs_Pcom < 0.001 ~ "***",
      p_value_Mcap_vs_Pcom < 0.01 ~ "**",
      p_value_Mcap_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Pacu_vs_Pcom = case_when(
      p_value_Pacu_vs_Pcom < 0.001 ~ "***",
      p_value_Pacu_vs_Pcom < 0.01 ~ "**",
      p_value_Pacu_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )


# Compute y.position dynamically for each comparison
ci_plot_Rd_data <- ci_plot_Rd_data %>%
  group_by(param) %>%
  mutate(
    y_max = max(conf_upper, na.rm = TRUE),
    y_min = min(conf_lower, na.rm = TRUE),
    y_range = y_max - y_min,
    y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
    y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
    y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
  ) %>%
  ungroup()


# Create the plot
ci_Rd_data_plot <- ggplot(ci_plot_Rd_data, aes(x = source, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  # Add bars and asterisks for significant comparisons
  geom_segment(
    aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
    data = ci_plot_Rd_data %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 1.5, y = y_position_Mcap_vs_Pacu + 0.000005, label = signif_label_Mcap_vs_Pacu),
    data = ci_plot_Rd_data %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
    data = ci_plot_Rd_data %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2, y = y_position_Mcap_vs_Pcom + 0.000005, label = signif_label_Mcap_vs_Pcom),
    data = ci_plot_Rd_data %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
    data = ci_plot_Rd_data %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.000005, label = signif_label_Pacu_vs_Pcom),
    data = ci_plot_Rd_data %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  theme_bw() +
  facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(
    title = 'Confidence intervals for model parameters with significance',
    y = 'Estimate',
    x = 'Source'
  ) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )



#### Function to test overlap of confidence intervals - EXTRA PARAMETERS
#ci_extra_params_residual_combined

test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
                            estimate2, conf_lower2, conf_upper2) {
  # Calculate the Z-score for the difference
  diff <- abs(estimate1 - estimate2)
  pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
                      ((conf_upper2 - conf_lower2) / 3.92)^2)
  z_score <- diff / pooled_se
  
  # Calculate p-value
  p_value <- 2 * (1 - pnorm(z_score))
  
  return(p_value)
}

# Apply the function to compare parameters between species
ci_comparison_Rd_extra_params <- ci_extra_params_Rd_residual_combined %>%
  group_by(param) %>%
  summarise(
    p_value_Mcap_vs_Pacu = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
    ),
    p_value_Mcap_vs_Pcom = test_ci_overlap(
      estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    ),
    p_value_Pacu_vs_Pcom = test_ci_overlap(
      estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
      estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
    )
  )

# View the results
print(ci_comparison_Rd_extra_params)


# Combine ci_extra_params_case_combined with ci_comparison_extra_params
ci_plot_data_Rd_extra <- ci_extra_params_Rd_residual_combined %>%
  left_join(ci_comparison_Rd_extra_params, by = "param")

# # Calculate max_conf_upper for each parameter
# ci_plot_data_extra <- ci_plot_data_extra %>%
#   group_by(param) %>%
#   mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
#   ungroup()

# Compute y.position values for each parameter and source
y_positions <- ci_plot_data_Rd_extra %>%
  group_by(param) %>%
  summarise(
    max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
    .groups = "drop"
  ) %>%
  mutate(
    y_base = max_value  # Add a base offset to the maximum value
  )

# Add y.position to annotations with dynamic spacing
ci_plot_data_Rd_extra <- ci_plot_data_Rd_extra %>%
  left_join(y_positions, by = "param") %>%
  group_by(param) %>%
  mutate(
    y.position = y_base + (row_number() - 1) * (y_base * 0.1)  # Dynamic offset based on y_base
  ) %>%
  ungroup() %>%
  select(-max_value, -y_base)  # Remove intermediate columns


# Add significance levels for each comparison
ci_plot_data_Rd_extra <- ci_plot_data_Rd_extra %>%
  mutate(
    signif_label_Mcap_vs_Pacu = case_when(
      p_value_Mcap_vs_Pacu < 0.001 ~ "***",
      p_value_Mcap_vs_Pacu < 0.01 ~ "**",
      p_value_Mcap_vs_Pacu < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Mcap_vs_Pcom = case_when(
      p_value_Mcap_vs_Pcom < 0.001 ~ "***",
      p_value_Mcap_vs_Pcom < 0.01 ~ "**",
      p_value_Mcap_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    ),
    signif_label_Pacu_vs_Pcom = case_when(
      p_value_Pacu_vs_Pcom < 0.001 ~ "***",
      p_value_Pacu_vs_Pcom < 0.01 ~ "**",
      p_value_Pacu_vs_Pcom < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )

# Compute y.position dynamically for each comparison
ci_plot_data_Rd_extra <- ci_plot_data_Rd_extra %>%
  group_by(param) %>%
  mutate(
    y_max = max(conf_upper, na.rm = TRUE),
    y_min = min(conf_lower, na.rm = TRUE),
    y_range = y_max - y_min,
    y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
    y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
    y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
  ) %>%
  ungroup()


# Create the plot
ci_data_Rd_extra_plot <- ggplot(ci_plot_data_Rd_extra, aes(x = source, y = estimate, col = source)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
  # Add bars and asterisks for significant comparisons
  geom_segment(
    aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
    data = ci_plot_data_Rd_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 1.5, y = y_position_Mcap_vs_Pacu + 0.000005, label = signif_label_Mcap_vs_Pacu),
    data = ci_plot_data_Rd_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
    data = ci_plot_data_Rd_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2, y = y_position_Mcap_vs_Pcom +  0.000005, label = signif_label_Mcap_vs_Pcom),
    data = ci_plot_data_Rd_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_segment(
    aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
    data = ci_plot_data_Rd_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  geom_text(
    aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.000005, label = signif_label_Pacu_vs_Pcom),
    data = ci_plot_data_Rd_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
    inherit.aes = FALSE, col = "black"
  ) +
  theme_bw() +
  facet_wrap(~param, scales = 'free_y') +
  scale_x_discrete('') +
  scale_color_manual(values = cols) +
  labs(
    title = 'Confidence intervals for extra parameters with significance',
    y = 'Estimate',
    x = 'Source'
  ) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


######### save z-score test results for Am and Rd estimates for supplementary table
library(writexl)

# Create a comprehensive Excel workbook with all TPC parameter comparisons
tpc_comparison_tables <- list(
  "Am_Model_Parameters" = ci_Am_comparison,
  "Am_Extra_Parameters" = ci_comparison_Am_extra_params,
  "Rd_Model_Parameters" = ci_Rd_comparison,
  "Rd_Extra_Parameters" = ci_comparison_Rd_extra_params
)

# Save to Excel with multiple sheets
write_xlsx(tpc_comparison_tables, "output/TPC_estimates_zScore_pvalues.xlsx")







# ############ Mcap CURVE FIT - AQY (Alpha)
# Mcap.df <- df_wide %>%  
#   filter(species=="Montipora capitata") 

# Mcap.df$Temp.Cat_numeric <- as.numeric(as.character(Mcap.df$Temp.Cat))

# # Calculate mean and sd values for AQY, grouped by Temp.Cat

# Mcap.df_mean <- Mcap.df %>%
#   dplyr::select(Temp.Cat, AQY, colony_id) %>%
#   group_by(Temp.Cat) %>%
#   summarise(., sd = sd(AQY, na.rm = TRUE), ave_rate = mean(AQY, na.rm = TRUE),
#             groups = 'drop')

# # View the resulting dataframe
# print(Mcap.df_mean)

# Mcap.df_mean$Temp.Cat <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Mcap.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          #                                    # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Mcap.df_mean$Temp.Cat), max(Mcap.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   dplyr::select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   dplyr::select(-fit) %>%
#   unnest(info) %>%
#   dplyr::select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


# ############ Sharpe Schoolfield 1981 model 
# #Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

# #### we want to run the model taking into account means and standard deviations of the PI parameters
# #https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
# #https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# # load packages
# library(boot)
# library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
# library(rTPC)
# library(nls.multstart)
# library(broom)
# library(tidyverse)
# library(patchwork)
# library(minpack.lm)


# Mcap.df_mean$temp <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #temp needs to be numeric

# Mcap.df_mean <- Mcap.df_mean %>%
#   mutate(Temp.Cat = temp) %>%
#   select(temp, everything(), -Temp.Cat)

# d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 1,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 1,
#                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

# # d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
# #   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                              data = .x,
# #                                              iter = c(4,4,4,4),
# #                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 10,
# #                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 10,
# #                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              supp_errors = 'Y',
# #                                              convergence_count = FALSE,
# #                                              # include weights here!
# #                                              modelweights = 1/sd)))


# # check they work
# d_fit_Mcap$weighted[[1]]

# # get predictions using augment
# newdata_Mcap <- tibble(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))
# d_preds_Mcap <- d_fit_Mcap %>%
#   mutate(., preds = map(weighted, augment, newdata = newdata_Mcap)) %>%
#   select(-weighted) %>%
#   unnest(preds)

# # plot
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Mcap) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
#   geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd')+
#   #title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-0.25, 3.5))



# ### Residual re-sampling weighted non-linear regression

# #So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
# #number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
# #is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
# #as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# # refit model using nlsLM
# fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
#                                     data = Mcap.df_mean,
#                                     start = coef(d_fit_Mcap$weighted[[1]]),
#                                     #lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     #upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     weights = 1/sd)

# # fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                     data = Mcap.df_mean,
# #                                     start = coef(d_fit_Mcap$weighted[[1]]),
# #                                     #lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     #upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     weights = 1/sd)

# #Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
# #are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# # perform residual bootstrap
# boot2_Mcap <- Boot(fit_nlsLM_Mcap, method = 'residual')

# # predict over new data
# boot2_preds_Mcap <- boot2_Mcap$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 18))

# # boot2_preds_Mcap <- boot2_Mcap$t %>%
# #   as.data.frame() %>%
# #   drop_na() %>%
# #   mutate(iter = 1:n()) %>%
# #   group_by_all() %>%
# #   do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
# #   ungroup() %>%
# #   mutate(pred = sharpeschoollow_1981(temp, r_tref,e,el,tl,tref = 25))

# # calculate bootstrapped confidence intervals
# boot2_conf_preds_Mcap <- boot2_preds_Mcap %>%
#   group_by(temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975),
#             .groups = 'drop')

# # plot bootstrapped CIs - 95% CIs
# p1 <- ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Mcap, col = 'black') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = 'black', alpha = 0.3) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
#   geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 10) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd',
#        title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-3, 3.5))

# p1



# ## Calculating confidence intervals of estimated and calculated parameters
# #As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
# #parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
# #confidence intervals (using confint-MASS) and asymptotic confidence intervals 
# #(using nlstools::confint2()).

# library(nlstools)

# # get parameters of fitted model
# param_Mcap <- broom::tidy(fit_nlsLM_Mcap) %>%
#   select(param = term, estimate)


# # calculate confidence intervals of models
# ci1_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'asymptotic') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'asymptotic')
# ci2_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'profile')
# #> Waiting for profiling to be done...
# #> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# # profiling method fails
# ci2_Mcap <- mutate(ci1_Mcap, method = 'profile',
#                    conf_lower = NA,
#                    conf_upper = NA)


# # CIs from residual resampling
# ci4_Mcap <- confint(boot2_Mcap, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# ## [1] "All values of t are equal to  0.165850014721406 \n Cannot calculate confidence intervals"

# ci_AQY_Mcap <- bind_rows(ci1_Mcap, ci2_Mcap, ci4_Mcap) %>%
#   full_join(., param_Mcap, by = "param")
# #> Joining with `by = join_by(param)`

# ggplot(ci_AQY_Mcap, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
#   geom_point(size = 1) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
#   labs(title = 'Calculation of confidence intervals for model parameters',
#        subtitle = 'For the chlorella TPC; profile method failes')

# #### asymptotic works


# # We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

# extra_params <- calc_params(fit_nlsLM_Mcap) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 400, method = 'case') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')

# ci_extra_params_AQY_case_Mcap <- left_join(ci_extra_params, extra_params)

# ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 400, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')

# ci_extra_params_AQY_residual_Mcap <- left_join(ci_extra_params, extra_params)


# ggplot(ci_extra_params_AQY_case_Mcap, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )

# ggplot(ci_extra_params_AQY_residual_Mcap, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )




# #### Pacu CURVE FIT - AQY
# Pacu.df <- df_wide %>%  
#   filter(species=="Pocillopora acuta") 

# Pacu.df$Temp.Cat_numeric <- as.numeric(as.character(Pacu.df$Temp.Cat))

# # Calculate mean and sd values for AQY, grouped by Temp.Cat

# Pacu.df_mean <- Pacu.df %>%
#   select(Temp.Cat, AQY, colony_id) %>%
#   group_by(Temp.Cat) %>%
#   summarise(., sd = sd(AQY,na.rm = TRUE),ave_rate = mean(AQY,na.rm = TRUE),
#             groups = 'drop')


# # View the resulting dataframe
# print(Pacu.df_mean)

# Pacu.df_mean$Temp.Cat <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits_Pacu <- nest(Pacu.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pacu.df_mean$Temp.Cat), max(Pacu.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


# ############ Sharpe Schoolfield 1981 model 
# #Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

# #### we want to run the model taking into account means and standard deviations of the PI parameters
# #https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
# #https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# # load packages
# library(boot)
# library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
# library(rTPC)
# library(nls.multstart)
# library(broom)
# library(tidyverse)
# library(patchwork)
# library(minpack.lm)


# Pacu.df_mean$temp <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #temp needs to be numeric

# Pacu.df_mean <- Pacu.df_mean %>%
#   mutate(Temp.Cat = temp) %>%
#   select(temp, everything(), -Temp.Cat)

# d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 1,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 1,
#                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

# # d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
# #   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp,r_tref,e,el,tl, tref = 25),
# #                                              data = .x,
# #                                              iter = c(4,4,4,4),
# #                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 1,
# #                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 1,
# #                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              supp_errors = 'Y',
# #                                              convergence_count = FALSE,
# #                                              # include weights here!
# #                                              modelweights = 1/sd)))

# # check they work
# d_fit_Pacu$weighted[[1]]


# # get predictions using augment
# newdata_Pacu <- tibble(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out = 100))
# d_preds_Pacu <- d_fit_Pacu %>%
#   mutate(., preds = map(weighted, augment, newdata = newdata_Pacu)) %>%
#   select(-weighted) %>%
#   unnest(preds)

# # plot
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pacu) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
#   geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd')+
#   #title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 



# ### Residual re-sampling weighted non-linear regression

# #So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
# #number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
# #is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
# #as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# # refit model using nlsLM
# fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
#                                     data = Pacu.df_mean,
#                                     start = coef(d_fit_Pacu$weighted[[1]]),
#                                     #lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     #upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     weights = 1/sd)

# # fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                     data = Pacu.df_mean,
# #                                     start = coef(d_fit_Pacu$weighted[[1]]),
# #                                     #lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     #upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     weights = 1/sd)

# #Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
# #are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# # perform residual bootstrap
# boot2_Pacu <- Boot(fit_nlsLM_Pacu, method = 'residual', R = 1000)

# # predict over new data
# boot2_preds_Pacu <- boot2_Pacu$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 18))

# # boot2_preds_Pacu <- boot2_Pacu$t %>%
# #   as.data.frame() %>%
# #   drop_na() %>%
# #   mutate(iter = 1:n()) %>%
# #   group_by_all() %>%
# #   do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
# #   ungroup() %>%
# #   mutate(pred = sharpeschoollow_1981(temp, r_tref,e,el,tl, tref = 25))


# # calculate bootstrapped confidence intervals
# boot2_conf_preds_Pacu <- boot2_preds_Pacu %>%
#   group_by(temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975),
#             .groups = 'drop')

# # plot bootstrapped CIs - 95% CIs
# p1 <- ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pacu, col = 'black') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'black', alpha = 0.3) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
#   geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 10) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd',
#        title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-3, 3.5))

# p1


# ## Calculating confidence intervals of estimated and calculated parameters
# #As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
# #parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
# #confidence intervals (using confint-MASS) and asymptotic confidence intervals 
# #(using nlstools::confint2()).

# library(nlstools)

# # get parameters of fitted model
# param_Pacu <- broom::tidy(fit_nlsLM_Pacu) %>%
#   select(param = term, estimate)


# # calculate confidence intervals of models
# ci1_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'asymptotic') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'asymptotic')
# ci2_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'profile')
# #> Waiting for profiling to be done...

# ci2_Pacu <- mutate(ci1_Pacu, method = 'profile',
#                    conf_lower = NA,
#                    conf_upper = NA)

# # CIs from case resampling
# ci3_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')


# # CIs from residual resampling
# ci4_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# # [1] "All values of t are equal to  0.210223059177217 \n Cannot calculate confidence intervals"

# ci_Pacu <- bind_rows(ci1_Pacu,ci2_Pacu,ci3_Pacu, ci4_Pacu) %>%
#   full_join(., param_Pacu , by = "param")
# ci_AQY_Pacu <- bind_rows(ci2_Pacu,ci4_Pacu) %>%
#   full_join(., param_Pacu , by = "param")
# #> Joining with `by = join_by(param)`

# ggplot(ci_AQY_Pacu , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
#   geom_point(size = 1) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
#   labs(title = 'Calculation of confidence intervals for model parameters'
#   )


# # We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

# extra_params <- calc_params(fit_nlsLM_Pacu) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 400, method = 'case') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')

# ci_extra_params_AQY_case_Pacu <- left_join(ci_extra_params, extra_params)

# ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 400, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')

# ci_extra_params_AQY_residual_Pacu <- left_join(ci_extra_params, extra_params)


# ggplot(ci_extra_params_AQY_residual_Pacu, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )

# ggplot(ci_extra_params_AQY_case_Pacu, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )



# #### Pcomp CURVE FIT - AQY
# Pcom.df <- df_wide %>%  
#   filter(species=="Porites compressa") 

# Pcom.df$Temp.Cat_numeric <- as.numeric(as.character(Pcom.df$Temp.Cat))

# # Calculate mean and sd values for AQY, grouped by Temp.Cat

# Pcom.df_mean <- Pcom.df %>%
#   select(Temp.Cat, AQY, colony_id) %>%
#   group_by(Temp.Cat) %>%
#   summarise(., sd = sd(AQY, na.rm = TRUE),ave_rate = mean(AQY, na.rm = TRUE),
#             groups = 'drop')


# # View the resulting dataframe
# print(Pcom.df_mean)

# Pcom.df_mean$Temp.Cat <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Pcom.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pcom.df_mean$Temp.Cat), max(Pcom.df_mean$Temp.Cat), length.out = 10))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


# ############ Sharpe Schoolfield 1981 model 
# #Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

# #### we want to run the model taking into account means and standard deviations of the PI parameters
# #https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
# #https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# # load packages
# library(boot)
# library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
# library(rTPC)
# library(nls.multstart)
# library(broom)
# library(tidyverse)
# library(patchwork)
# library(minpack.lm)


# Pcom.df_mean$temp <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #temp needs to be numeric

# Pcom.df_mean <- Pcom.df_mean %>%
#   mutate(Temp.Cat = temp) %>%
#   select(temp, everything(), -Temp.Cat)

# d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 1,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 1,
#                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

# # d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
# #   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                              data = .x,
# #                                              iter = c(4,4,4,4),
# #                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 10,
# #                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 10,
# #                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              supp_errors = 'Y',
# #                                              convergence_count = FALSE,
# #                                              # include weights here!
# #                                              modelweights = 1/sd)))

# # check they work
# d_fit_Pcom$weighted[[1]]



# # get predictions using augment
# newdata_Pcom <- tibble(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))
# d_preds_Pcom <- d_fit_Pcom %>%
#   mutate(., preds = map(weighted, augment, newdata = newdata_Pcom)) %>%
#   select(-weighted) %>%
#   unnest(preds)

# # plot
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pcom) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
#   geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd')+
#   #title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 



# ### Residual re-sampling weighted non-linear regression

# #So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
# #number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
# #is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
# #as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# # refit model using nlsLM
# fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 18),
#                                     data = Pcom.df_mean,
#                                     start = coef(d_fit_Pcom$weighted[[1]]),
#                                     #lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     #upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     weights = 1/sd)

# # refit model using nlsLM
# # fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp,r_tref,e,el,tl, tref = 25),
# #                                     data = Pcom.df_mean,
# #                                     start = coef(d_fit_Pcom$weighted[[1]]),
# #                                     #lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     #upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     weights = 1/sd)

# #Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
# #are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# # perform residual bootstrap
# boot2_Pcom <- Boot(fit_nlsLM_Pcom, method = 'residual')


# # predict over new data
# boot2_preds_Pcom <- boot2_Pcom$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 18))

# # boot2_preds_Pcom <- boot2_Pcom$t %>%
# #   as.data.frame() %>%
# #   drop_na() %>%
# #   mutate(iter = 1:n()) %>%
# #   group_by_all() %>%
# #   do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
# #   ungroup() %>%
# #   mutate(pred = sharpeschoollow_1981(temp, r_tref,e,el,tl, tref = 25))

# # calculate bootstrapped confidence intervals
# boot2_conf_preds_Pcom <- boot2_preds_Pcom %>%
#   group_by(temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975),
#             .groups = 'drop')

# # plot bootstrapped CIs - 95% CIs
# p1 <- ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pcom, col = 'black') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'black', alpha = 0.3) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
#   geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 10) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd',
#        title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-3, 3.5))

# p1


# ## Calculating confidence intervals of estimated and calculated parameters
# #As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
# #parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
# #confidence intervals (using confint-MASS) and asymptotic confidence intervals 
# #(using nlstools::confint2()).

# library(nlstools)

# # get parameters of fitted model
# param_Pcom <- broom::tidy(fit_nlsLM_Pcom) %>%
#   select(param = term, estimate)


# # calculate confidence intervals of models
# ci1_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'asymptotic') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'asymptotic')
# ci2_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'profile')
# #> Waiting for profiling to be done...
# #> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# # profiling method fails
# ci2_Pcom <- mutate(ci1_Pcom, method = 'profile',
#                    conf_lower = NA,
#                    conf_upper = NA)

# # CIs from case resampling
# ci3 <- confint(boot2_Pcom, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')


# # CIs from residual resampling
# ci4_Pcom <- confint(boot2_Pcom, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# # [1] "All values of t are equal to  0.267937013495708 \n Cannot calculate confidence intervals"


# ci_AQY_Pcom <- bind_rows(ci1_Pcom, ci2_Pcom, ci4_Pcom  ) %>%
#   full_join(., param_Pcom, by = "param")

# #> Joining with `by = join_by(param)`

# ggplot(ci_AQY_Pcom , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
#   geom_point(size = 1) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
#   labs(title = 'Calculation of confidence intervals for model parameters'
#   )


# # We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

# extra_params <- calc_params(fit_nlsLM_Pcom) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 400, method = 'case') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')

# ci_extra_params_AQY_case_Pcom<- left_join(ci_extra_params, extra_params)

# ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 400, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')

# ci_extra_params_AQY_residual_Pcom<- left_join(ci_extra_params, extra_params)


# ggplot(ci_extra_params_AQY_residual_Pcom, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )

# #filter to only the most relavent and well characterized parameters 
# # All_params <- All_params %>%  
# #   filter(!param=="ctmin") %>% 
# #   filter(!param=="ctmax") %>% 
# #   filter(!param=="skewness") %>% 
# #   filter(!param=="thermal_safety_margin") %>% 
# #   filter(!param=="thermal_tolerance") %>% 
# #   filter(!param=="q10")%>% 
# #   filter(!param=="breadth") 



# #set plot colors 
# cols <- c("Mcap"="green",  
#           "Pacu"="cyan", 
#           "Pcom"="orange") 

# # Ensure Temp.Cat is numeric
# df_wide$Temp.Cat <- as.numeric(as.character(df_wide$Temp.Cat))
# Mcap.df$Temp.Cat <- as.numeric(as.character(Mcap.df$Temp.Cat))
# Pacu.df$Temp.Cat <- as.numeric(as.character(Pacu.df$Temp.Cat))
# Pcom.df$Temp.Cat <- as.numeric(as.character(Pcom.df$Temp.Cat))

# # plot data and model fit for the 3 species
# TPC.plot_AQY <- ggplot(data=df_wide, aes()) + 
#   geom_point(aes(Temp.Cat,AQY, color = "Mcap"), data =  Mcap.df, size = 2, alpha = 0.5) + 
#   geom_point(aes(Temp.Cat,AQY,  color = "Pacu"), data = Pacu.df, size = 2, alpha = 0.5) + 
#   geom_point(aes(Temp.Cat,AQY, color = "Pcom"), data = Pcom.df, size = 2, alpha = 0.5) + 
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = "green", alpha = 0.3) + 
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'cyan', alpha = 0.3) + 
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'orange', alpha = 0.3) + 
#   geom_line(aes(temp, .fitted), data = d_preds_Mcap, color = "green") +
#   geom_line(aes(temp, .fitted), data = d_preds_Pacu, color = "cyan") +
#   geom_line(aes(temp, .fitted), data = d_preds_Pcom, color = "orange") +
#   xlim(11.5,40.5)+ 
#   scale_x_continuous(breaks=c(12,14,16,18,20,22,24,26,28,30,32,34,36))+
#   theme_bw(base_size = 12) + 
#   scale_colour_manual(name="Species",values=cols)+ 
#   scale_fill_manual(name = "Species", values = cols) + 
#   theme(legend.position = "top", 
#         panel.border = element_blank(), panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
#   labs(x = 'Temperature (ºC)', 
#        y = expression("AQY"))

# TPC.plot_AQY  



# #### Plot TPC parameters
# # Combine the data frames and add a source column
# ci_extra_params_AQY_case_Mcap <- ci_extra_params_AQY_case_Mcap %>% mutate(source = "Mcap")
# ci_extra_params_AQY_case_Pacu <- ci_extra_params_AQY_case_Pacu %>% mutate(source = "Pacu")
# ci_extra_params_AQY_case_Pcom <- ci_extra_params_AQY_case_Pcom %>% mutate(source = "Pcom")

# ci_extra_params_AQY_residual_Mcap <- ci_extra_params_AQY_residual_Mcap %>% mutate(source = "Mcap")
# ci_extra_params_AQY_residual_Pacu <- ci_extra_params_AQY_residual_Pacu %>% mutate(source = "Pacu")
# ci_extra_params_AQY_residual_Pcom <- ci_extra_params_AQY_residual_Pcom %>% mutate(source = "Pcom")

# # Combine all data frames into one
# ci_extra_params_AQY_case_combined <- bind_rows(ci_extra_params_AQY_case_Mcap , ci_extra_params_AQY_case_Pacu, ci_extra_params_AQY_case_Pcom )

# ci_extra_params_AQY_residual_combined <- bind_rows(ci_extra_params_AQY_residual_Mcap , ci_extra_params_AQY_residual_Pacu, ci_extra_params_AQY_residual_Pcom )

# #set plot colors 
# cols <- c("Mcap"="green",  
#           "Pacu"="cyan", 
#           "Pcom"="orange") 

# # Plot the combined data
# ggplot(ci_extra_params_AQY_residual_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 1, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   #facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for extra parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

# # Plot the combined data
# ggplot(ci_extra_params_AQY_case_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 1, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   #facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for extra parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



# # Add a source column to each data frame
# ci_AQY_Mcap <- ci_AQY_Mcap %>% mutate(source = "Mcap")
# ci_AQY_Pacu <- ci_AQY_Pacu %>% mutate(source = "Pacu")
# ci_AQY_Pcom <- ci_AQY_Pcom %>% mutate(source = "Pcom")

# # Combine all data frames into one
# ci_AQY_combined <- bind_rows(ci_AQY_Mcap, ci_AQY_Pacu, ci_AQY_Pcom )

# #set plot colors 
# cols <- c("Mcap"="green",  
#           "Pacu"="cyan", 
#           "Pcom"="orange") 

# # Plot the combined data
# ggplot(ci_AQY_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for model parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ggplot(ci_AQY_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(param ~ method, scales = 'free') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for model parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# # Filter rows with "residual bootstrap" in the method column
# ci_residual_AQY_combined <- ci_AQY_combined %>%
#   filter(method == "residual bootstrap")

# ci_residual_AQY_combined$param = as.factor(ci_residual_AQY_combined$param)

# ggplot(ci_residual_AQY_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap("param", scales = 'free') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Confidence intervals for model parameters - residual bootstrap',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




# ######## Statistical analysis of TPC parameters
# ### we have already calculated the confidence intervals for the parameters, so we can use these to calculate the
# ### statistical differences between the species.

# # Filter rows with "residual bootstrap" in the method column
# ci_residual_AQY_combined <- ci_AQY_combined %>%
#   filter(method == "residual bootstrap")

# # View the new table
# print(ci_residual_AQY_combined)

# #Formally test whether confidence intervals overlap using a statistical approach called confidence interval comparison. One common method 
# #is to calculate the Z-score for the difference between two estimates and determine whether the difference is statistically significant.
# #The value 3.92 below is used to approximate the standard error of the confidence interval (CI) bounds. It comes from the fact that a 95% 
# #confidence interval corresponds to approximately 1.96 standard deviations above and below the mean in a normal distribution. Since the 
# #CI width spans both directions (upper and lower bounds), the total width is 2 * 1.96 = 3.92.

# # Load required libraries
# library(dplyr)

# # Function to test overlap of confidence intervals
# test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
#                             estimate2, conf_lower2, conf_upper2) {
#   # Calculate the Z-score for the difference
#   diff <- abs(estimate1 - estimate2)
#   pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
#                       ((conf_upper2 - conf_lower2) / 3.92)^2)
#   z_score <- diff / pooled_se
  
#   # Calculate p-value
#   p_value <- 2 * (1 - pnorm(z_score))
  
#   return(p_value)
# }

# # Apply the function to compare parameters between species
# ci_AQY_comparison <- ci_residual_AQY_combined %>%
#   group_by(param) %>%
#   summarise(
#     p_value_Mcap_vs_Pacu = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
#     ),
#     p_value_Mcap_vs_Pcom = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     ),
#     p_value_Pacu_vs_Pcom = test_ci_overlap(
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     )
#   )

# # View the results
# print(ci_AQY_comparison)


# #The 95 % confidence intervals from two subgroups or studies may overlap substantially and yet the test 
# #for difference between them may still produce P < 0.05. 
# #https://pmc.ncbi.nlm.nih.gov/articles/PMC4877414/
# # Confidence Interval Misinterpretation: Confidence intervals are often misinterpreted as a direct test 
# # of significance. However, overlapping confidence intervals do not necessarily mean the difference is not 
# # significant. This is because the overlap does not account for the combined uncertainty of both estimates.

# # Calculate max conf_upper for each parameter
# ci_residual_AQY_combined <- ci_residual_AQY_combined %>%
#   group_by(param) %>%
#   mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
#   ungroup()

# # Combine ci_residual_combined with ci_comparison
# ci_plot_AQY_data <- ci_residual_AQY_combined %>%
#   left_join(ci_AQY_comparison, by = "param")

# # Prepare the statistical annotations for ci_plot_data
# # Compute y.position values for each parameter and source
# y_positions <- ci_plot_AQY_data %>%
#   group_by(param) %>%
#   summarise(
#     max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
#     .groups = "drop"
#   ) %>%
#   mutate(
#     y_base = max_value  # Add a base offset to the maximum value
#   )

# # Add y.position to annotations with dynamic spacing
# ci_plot_AQY_data <- ci_plot_AQY_data %>%
#   left_join(y_positions, by = "param") %>%
#   group_by(param) %>%
#   mutate(
#     y.position = y_base + (row_number() - 1) * (y_base * 0.4)  # Dynamic offset based on y_base
#   ) %>%
#   ungroup() %>%
#   select(-max_value, -y_base)  # Remove intermediate columns

# # # Filter for significant comparisons only
# # ci_plot_data <- ci_plot_data %>%
# #   filter(p_value_Mcap_vs_Pacu < 0.05 | p_value_Mcap_vs_Pcom < 0.05 | p_value_Pacu_vs_Pcom < 0.05)

# # Add significance levels for each comparison
# ci_plot_AQY_data <- ci_plot_AQY_data %>%
#   mutate(
#     signif_label_Mcap_vs_Pacu = case_when(
#       p_value_Mcap_vs_Pacu < 0.001 ~ "***",
#       p_value_Mcap_vs_Pacu < 0.01 ~ "**",
#       p_value_Mcap_vs_Pacu < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Mcap_vs_Pcom = case_when(
#       p_value_Mcap_vs_Pcom < 0.001 ~ "***",
#       p_value_Mcap_vs_Pcom < 0.01 ~ "**",
#       p_value_Mcap_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Pacu_vs_Pcom = case_when(
#       p_value_Pacu_vs_Pcom < 0.001 ~ "***",
#       p_value_Pacu_vs_Pcom < 0.01 ~ "**",
#       p_value_Pacu_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     )
#   )


# # Compute y.position dynamically for each comparison
# ci_plot_AQY_data <- ci_plot_AQY_data %>%
#   group_by(param) %>%
#   mutate(
#     y_max = max(conf_upper, na.rm = TRUE),
#     y_min = min(conf_lower, na.rm = TRUE),
#     y_range = y_max - y_min,
#     y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
#     y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
#     y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
#   ) %>%
#   ungroup()

# # Create the plot
# ci_AQY_data_plot <- ggplot(ci_plot_AQY_data, aes(x = source, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   # Add bars and asterisks for significant comparisons
#   geom_segment(
#     aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
#     data = ci_plot_AQY_data %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 1.5, y = y_position_Mcap_vs_Pacu +0.000005, label = signif_label_Mcap_vs_Pacu),
#     data = ci_plot_AQY_data %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
#     data = ci_plot_AQY_data %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2, y = y_position_Mcap_vs_Pcom + 0.000005, label = signif_label_Mcap_vs_Pcom),
#     data = ci_plot_AQY_data %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
#     data = ci_plot_AQY_data %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.000005, label = signif_label_Pacu_vs_Pcom),
#     data = ci_plot_AQY_data %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(
#     title = 'Confidence intervals for model parameters with significance',
#     y = 'Estimate',
#     x = 'Source'
#   ) +
#   theme(
#     legend.position = "top",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black")
#   )



# #### Function to test overlap of confidence intervals - EXTRA PARAMETERS
# #ci_extra_params_residual_combined

# test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
#                             estimate2, conf_lower2, conf_upper2) {
#   # Calculate the Z-score for the difference
#   diff <- abs(estimate1 - estimate2)
#   pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
#                       ((conf_upper2 - conf_lower2) / 3.92)^2)
#   z_score <- diff / pooled_se
  
#   # Calculate p-value
#   p_value <- 2 * (1 - pnorm(z_score))
  
#   return(p_value)
# }

# # Apply the function to compare parameters between species
# ci_comparison_extra_AQY_params <- ci_extra_params_AQY_residual_combined %>%
#   group_by(param) %>%
#   summarise(
#     p_value_Mcap_vs_Pacu = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
#     ),
#     p_value_Mcap_vs_Pcom = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     ),
#     p_value_Pacu_vs_Pcom = test_ci_overlap(
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     )
#   )

# # View the results
# print(ci_comparison_extra_AQY_params)


# # Combine ci_extra_params_case_combined with ci_comparison_extra_params
# ci_plot_data_AQY_extra <- ci_extra_params_AQY_residual_combined %>%
#   left_join(ci_comparison_extra_AQY_params, by = "param")

# # # Calculate max_conf_upper for each parameter
# # ci_plot_data_extra <- ci_plot_data_extra %>%
# #   group_by(param) %>%
# #   mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
# #   ungroup()

# # Compute y.position values for each parameter and source
# y_positions <- ci_plot_data_AQY_extra %>%
#   group_by(param) %>%
#   summarise(
#     max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
#     .groups = "drop"
#   ) %>%
#   mutate(
#     y_base = max_value  # Add a base offset to the maximum value
#   )

# # Add y.position to annotations with dynamic spacing
# ci_plot_data_AQY_extra <- ci_plot_data_AQY_extra %>%
#   left_join(y_positions, by = "param") %>%
#   group_by(param) %>%
#   mutate(
#     y.position = y_base + (row_number() - 1) * (y_base * 0.1)  # Dynamic offset based on y_base
#   ) %>%
#   ungroup() %>%
#   select(-max_value, -y_base)  # Remove intermediate columns

# # Add significance levels for each comparison
# ci_plot_data_AQY_extra <- ci_plot_data_AQY_extra %>%
#   mutate(
#     signif_label_Mcap_vs_Pacu = case_when(
#       p_value_Mcap_vs_Pacu < 0.001 ~ "***",
#       p_value_Mcap_vs_Pacu < 0.01 ~ "**",
#       p_value_Mcap_vs_Pacu < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Mcap_vs_Pcom = case_when(
#       p_value_Mcap_vs_Pcom < 0.001 ~ "***",
#       p_value_Mcap_vs_Pcom < 0.01 ~ "**",
#       p_value_Mcap_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Pacu_vs_Pcom = case_when(
#       p_value_Pacu_vs_Pcom < 0.001 ~ "***",
#       p_value_Pacu_vs_Pcom < 0.01 ~ "**",
#       p_value_Pacu_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     )
#   )

# # Compute y.position dynamically for each comparison
# ci_plot_data_AQY_extra <- ci_plot_data_AQY_extra %>%
#   group_by(param) %>%
#   mutate(
#     y_max = max(conf_upper, na.rm = TRUE),
#     y_min = min(conf_lower, na.rm = TRUE),
#     y_range = y_max - y_min,
#     y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
#     y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
#     y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
#   ) %>%
#   ungroup()

# # Create the plot
# ci_data_AQY_extra_plot <- ggplot(ci_plot_data_AQY_extra, aes(x = source, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   # Add bars and asterisks for significant comparisons
#   geom_segment(
#     aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
#     data = ci_plot_data_AQY_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 1.5, y = y_position_Mcap_vs_Pacu +0.0000005, label = signif_label_Mcap_vs_Pacu),
#     data = ci_plot_data_AQY_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
#     data = ci_plot_data_AQY_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2, y = y_position_Mcap_vs_Pcom +0.0000005, label = signif_label_Mcap_vs_Pcom),
#     data = ci_plot_data_AQY_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
#     data = ci_plot_data_AQY_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2.5, y = y_position_Pacu_vs_Pcom +0.0000005, label = signif_label_Pacu_vs_Pcom),
#     data = ci_plot_data_AQY_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(
#     title = 'Confidence intervals for extra parameters with significance',
#     y = 'Estimate',
#     x = 'Source'
#   ) +
#   theme(
#     legend.position = "top",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black")
#   )






# ############ Mcap CURVE FIT - Ik
# Mcap.df <- df_wide %>%  
#   filter(species=="Montipora capitata") 

# Mcap.df$Temp.Cat_numeric <- as.numeric(as.character(Mcap.df$Temp.Cat))

# # Calculate mean and sd values for Ik, grouped by Temp.Cat

# Mcap.df_mean <- Mcap.df %>%
#   dplyr::select(Temp.Cat, Ik, colony_id) %>%
#   group_by(Temp.Cat) %>%
#   summarise(., sd = sd(Ik, na.rm = TRUE), ave_rate = mean(Ik, na.rm = TRUE),
#             groups = 'drop')

# # View the resulting dataframe
# print(Mcap.df_mean)

# Mcap.df_mean$Temp.Cat <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Mcap.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          #                                    # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Mcap.df_mean$Temp.Cat), max(Mcap.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   dplyr::select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   dplyr::select(-fit) %>%
#   unnest(info) %>%
#   dplyr::select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Mcap.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


# ############ Sharpe Schoolfield 1981 model 
# #Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

# #### we want to run the model taking into account means and standard deviations of the PI parameters
# #https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
# #https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# # load packages
# library(boot)
# library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
# library(rTPC)
# library(nls.multstart)
# library(broom)
# library(tidyverse)
# library(patchwork)
# library(minpack.lm)


# Mcap.df_mean$temp <- as.numeric(as.character(Mcap.df_mean$Temp.Cat)) #temp needs to be numeric

# Mcap.df_mean <- Mcap.df_mean %>%
#   mutate(Temp.Cat = temp) %>%
#   select(temp, everything(), -Temp.Cat)

# d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 12),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 1,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 1,
#                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

# # d_fit_Mcap <- nest(Mcap.df_mean, data = c(ave_rate, temp, sd)) %>%
# #   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                              data = .x,
# #                                              iter = c(4,4,4,4),
# #                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 10,
# #                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 10,
# #                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              supp_errors = 'Y',
# #                                              convergence_count = FALSE,
# #                                              # include weights here!
# #                                              modelweights = 1/sd)))

# # check they work
# d_fit_Mcap$weighted[[1]]

# # get predictions using augment
# newdata_Mcap <- tibble(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))
# d_preds_Mcap <- d_fit_Mcap %>%
#   mutate(., preds = map(weighted, augment, newdata = newdata_Mcap)) %>%
#   select(-weighted) %>%
#   unnest(preds)

# # plot
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Mcap) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
#   geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd')+
#   #title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-0.25, 3.5))



# ### Residual re-sampling weighted non-linear regression

# #So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
# #number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
# #is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
# #as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# # refit model using nlsLM
# fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 12),
#                                     data = Mcap.df_mean,
#                                     start = coef(d_fit_Mcap$weighted[[1]]),
#                                     #lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     #upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     weights = 1/sd)

# # fit_nlsLM_Mcap <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                     data = Mcap.df_mean,
# #                                     start = coef(d_fit_Mcap$weighted[[1]]),
# #                                     #lower = get_lower_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     #upper = get_upper_lims(Mcap.df_mean$temp, Mcap.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     weights = 1/sd)

# #Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
# #are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# # perform residual bootstrap
# boot2_Mcap <- Boot(fit_nlsLM_Mcap, method = 'residual', R = 100)

# # predict over new data
# boot2_preds_Mcap <- boot2_Mcap$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Mcap.df_mean$temp), max(Mcap.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 12))

# # calculate bootstrapped confidence intervals
# boot2_conf_preds_Mcap <- boot2_preds_Mcap %>%
#   group_by(temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975),
#             .groups = 'drop')

# # plot bootstrapped CIs - 95% CIs
# p1 <- ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Mcap, col = 'black') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = 'black', alpha = 0.3) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Mcap.df_mean) +
#   geom_point(aes(temp, ave_rate), Mcap.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 10) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd',
#        title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-3, 3.5))

# p1



# ## Calculating confidence intervals of estimated and calculated parameters
# #As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
# #parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
# #confidence intervals (using confint-MASS) and asymptotic confidence intervals 
# #(using nlstools::confint2()).

# library(nlstools)

# # get parameters of fitted model
# param_Mcap <- broom::tidy(fit_nlsLM_Mcap) %>%
#   select(param = term, estimate)


# # calculate confidence intervals of models
# ci1_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'asymptotic') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'asymptotic')
# ci2_Mcap <- nlstools::confint2(fit_nlsLM_Mcap, method = 'profile')
# #> Waiting for profiling to be done...
# #> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# # profiling method fails
# ci2_Mcap <- mutate(ci1_Mcap, method = 'profile',
#                    conf_lower = NA,
#                    conf_upper = NA)


# # CIs from residual resampling
# ci4_Mcap <- confint(boot2_Mcap, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# ## [1] "All values of t are equal to  0.165850014721406 \n Cannot calculate confidence intervals"

# ci_Ik_Mcap <- bind_rows(ci1_Mcap, ci2_Mcap, ci4_Mcap) %>%
#   full_join(., param_Mcap, by = "param")
# #> Joining with `by = join_by(param)`

# ggplot(ci_Ik_Mcap, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
#   geom_point(size = 1) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
#   labs(title = 'Calculation of confidence intervals for model parameters',
#        subtitle = 'For the chlorella TPC; profile method failes')

# #### asymptotic works


# # We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

# extra_params <- calc_params(fit_nlsLM_Mcap) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 400, method = 'case') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')

# ci_extra_params_Ik_case_Mcap <- left_join(ci_extra_params, extra_params)

# ci_extra_params <- Boot(fit_nlsLM_Mcap, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Mcap)), R = 400, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')

# ci_extra_params_Ik_residual_Mcap <- left_join(ci_extra_params, extra_params)


# ggplot(ci_extra_params_Ik_case_Mcap, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )

# ggplot(ci_extra_params_Ik_residual_Mcap, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )




# #### Pacu CURVE FIT - Ik
# Pacu.df <- df_wide %>%  
#   filter(species=="Pocillopora acuta") 

# Pacu.df$Temp.Cat_numeric <- as.numeric(as.character(Pacu.df$Temp.Cat))

# # Calculate mean and sd values for Ik, grouped by Temp.Cat

# Pacu.df_mean <- Pacu.df %>%
#   select(Temp.Cat, Ik, colony_id) %>%
#   group_by(Temp.Cat) %>%
#   summarise(., sd = sd(Ik,na.rm = TRUE),ave_rate = mean(Ik,na.rm = TRUE),
#             groups = 'drop')


# # View the resulting dataframe
# print(Pacu.df_mean)

# Pacu.df_mean$Temp.Cat <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits_Pacu <- nest(Pacu.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pacu.df_mean$Temp.Cat), max(Pacu.df_mean$Temp.Cat), length.out = 6))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pacu.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


# ############ Sharpe Schoolfield 1981 model 
# #Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

# #### we want to run the model taking into account means and standard deviations of the PI parameters
# #https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
# #https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# # load packages
# library(boot)
# library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
# library(rTPC)
# library(nls.multstart)
# library(broom)
# library(tidyverse)
# library(patchwork)
# library(minpack.lm)


# Pacu.df_mean$temp <- as.numeric(as.character(Pacu.df_mean$Temp.Cat)) #temp needs to be numeric

# Pacu.df_mean <- Pacu.df_mean %>%
#   mutate(Temp.Cat = temp) %>%
#   select(temp, everything(), -Temp.Cat)

# d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 12),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 10,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 10,
#                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

# # d_fit_Pacu <- nest(Pacu.df_mean, data = c(ave_rate, temp, sd)) %>%
# #   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp,r_tref,e,el,tl, tref = 25),
# #                                              data = .x,
# #                                              iter = c(4,4,4,4),
# #                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 1,
# #                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 1,
# #                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              supp_errors = 'Y',
# #                                              convergence_count = FALSE,
# #                                              # include weights here!
# #                                              modelweights = 1/sd)))

# # check they work
# d_fit_Pacu$weighted[[1]]


# # get predictions using augment
# newdata_Pacu <- tibble(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out = 100))
# d_preds_Pacu <- d_fit_Pacu %>%
#   mutate(., preds = map(weighted, augment, newdata = newdata_Pacu)) %>%
#   select(-weighted) %>%
#   unnest(preds)

# # plot
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pacu) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
#   geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd')+
#   #title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 



# ### Residual re-sampling weighted non-linear regression

# #So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
# #number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
# #is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
# #as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# # refit model using nlsLM
# fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 12),
#                                     data = Pacu.df_mean,
#                                     start = coef(d_fit_Pacu$weighted[[1]]),
#                                     #lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     #upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     weights = 1/sd)

# # fit_nlsLM_Pacu <- minpack.lm::nlsLM(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                     data = Pacu.df_mean,
# #                                     start = coef(d_fit_Pacu$weighted[[1]]),
# #                                     #lower = get_lower_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     #upper = get_upper_lims(Pacu.df_mean$temp, Pacu.df_mean$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                     weights = 1/sd)

# #Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
# #are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# # perform residual bootstrap
# boot2_Pacu <- Boot(fit_nlsLM_Pacu, method = 'residual', R = 100)

# # predict over new data
# boot2_preds_Pacu <- boot2_Pacu$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pacu.df_mean$temp), max(Pacu.df_mean$temp), length.out =100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 12))


# # calculate bootstrapped confidence intervals
# boot2_conf_preds_Pacu <- boot2_preds_Pacu %>%
#   group_by(temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975),
#             .groups = 'drop')

# # plot bootstrapped CIs - 95% CIs
# p1 <- ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pacu, col = 'black') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'black', alpha = 0.3) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pacu.df_mean) +
#   geom_point(aes(temp, ave_rate), Pacu.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 10) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd',
#        title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-3, 3.5))

# p1


# ## Calculating confidence intervals of estimated and calculated parameters
# #As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
# #parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
# #confidence intervals (using confint-MASS) and asymptotic confidence intervals 
# #(using nlstools::confint2()).

# library(nlstools)

# # get parameters of fitted model
# param_Pacu <- broom::tidy(fit_nlsLM_Pacu) %>%
#   select(param = term, estimate)


# # calculate confidence intervals of models
# ci1_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'asymptotic') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'asymptotic')
# ci2_Pacu <- nlstools::confint2(fit_nlsLM_Pacu, method = 'profile')
# #> Waiting for profiling to be done...

# ci2_Pacu <- mutate(ci1_Pacu, method = 'profile',
#                    conf_lower = NA,
#                    conf_upper = NA)

# # CIs from case resampling
# ci3_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')


# # CIs from residual resampling
# ci4_Pacu <- confint(boot2_Pacu, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# # [1] "All values of t are equal to  0.210223059177217 \n Cannot calculate confidence intervals"

# ci_Pacu <- bind_rows(ci1_Pacu,ci2_Pacu,ci3_Pacu, ci4_Pacu) %>%
#   full_join(., param_Pacu , by = "param")
# ci_Ik_Pacu <- bind_rows(ci2_Pacu,ci4_Pacu) %>%
#   full_join(., param_Pacu , by = "param")
# #> Joining with `by = join_by(param)`

# ggplot(ci_Ik_Pacu , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
#   geom_point(size = 1) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
#   labs(title = 'Calculation of confidence intervals for model parameters'
#   )


# # We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

# extra_params <- calc_params(fit_nlsLM_Pacu) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 400, method = 'case') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')

# ci_extra_params_Ik_case_Pacu <- left_join(ci_extra_params, extra_params)

# ci_extra_params <- Boot(fit_nlsLM_Pacu, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pacu)), R = 400, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')

# ci_extra_params_Ik_residual_Pacu <- left_join(ci_extra_params, extra_params)


# ggplot(ci_extra_params_Ik_residual_Pacu, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )

# ggplot(ci_extra_params_Ik_case_Pacu, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )



# #### Pcomp CURVE FIT - Ik
# Pcom.df <- df_wide %>%  
#   filter(species=="Porites compressa") 

# Pcom.df$Temp.Cat_numeric <- as.numeric(as.character(Pcom.df$Temp.Cat))

# # Calculate mean and sd values for Ik, grouped by Temp.Cat

# Pcom.df_mean <- Pcom.df %>%
#   select(Temp.Cat, Ik, colony_id) %>%
#   group_by(Temp.Cat) %>%
#   summarise(., sd = sd(Ik, na.rm = TRUE),ave_rate = mean(Ik, na.rm = TRUE),
#             groups = 'drop')


# # View the resulting dataframe
# print(Pcom.df_mean)

# Pcom.df_mean$Temp.Cat <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #Temp.Cat needs to be numeric

# # fit 4 chosen model formulations in rTPC, taking into consideration standard deviations
# #https://padpadpadpad.github.io/rTPC/articles/model_averaging_selection.html
# #models to try, DeLong, sharpeschoolfull_1981, sharpeschoolhigh_1981, deutsch_2008
# d_fits <- nest(Pcom.df_mean, data = c(Temp.Cat, ave_rate, sd)) %>%
#   mutate(deutsch = map(data, ~nls_multstart(ave_rate~deutsch_2008(temp = Temp.Cat, rmax, topt, ctmax, a),
#                                             data = .x,
#                                             iter = c(4,4,4,4),
#                                             start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') - 10,
#                                             start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008') + 10,
#                                             lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'deutsch_2008'),
#                                             supp_errors = 'Y',
#                                             convergence_count = FALSE,
#                                             modelweights = 1/sd)),
#          # include weights here!         
#          # delong = map(data, ~nls_multstart(ave_rate~delong_2017(temp = Temp.Cat, c, eb, ef, tm, ehc),
#          #                                   data = .x,
#          #                                   iter = c(4,4,4,4,4),
#          #                                   start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') - 10,
#          #                                   start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017') + 10,
#          #                                   lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'delong_2017'),
#          #                                   supp_errors = 'Y',
#          #                                   convergence_count = FALSE,
#          #                                   modelweights = 1/sd)),
#          sharpeschoolfull = map(data, ~nls_multstart(ave_rate~sharpeschoolfull_1981(temp = Temp.Cat, r_tref,e,el,tl,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolfull_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)),
#          sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = Temp.Cat, r_tref,e,eh,th, tref = 18),
#                                                      data = .x,
#                                                      iter = c(4,4,4,4),
#                                                      start_lower = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                                                      start_upper = get_start_vals(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                                                      lower = get_lower_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      upper = get_upper_lims(.x$Temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
#                                                      supp_errors = 'Y',
#                                                      convergence_count = FALSE,
#                                                      modelweights = 1/sd)))


# # stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', deutsch:sharpeschoolhigh)

# # get predictions using augment
# newdata <- tibble(temp = seq(min(Pcom.df_mean$Temp.Cat), max(Pcom.df_mean$Temp.Cat), length.out = 10))
# d_preds <- d_stack %>%
#   mutate(., preds = map(fit, augment, newdata = newdata)) %>%
#   select(-fit) %>%
#   unnest(preds)

# # take a random point from each model for labelling
# d_labs <- filter(d_preds, temp < 30) %>%
#   group_by(., model_name) %>%
#   sample_n(., 1) %>%
#   ungroup()

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   scale_color_brewer(type = 'qual', palette = 2)


# # As can be seen in the above plot, there is some variation in how the different model formulations fit to 
# # the data. We can use a information theoretic approach to compare between different models, using measures
# # of relative model fit - such as AIC, BIC, and AICc (AIC correcting for small sample size) 

# library(MuMIn)

# d_ic <- d_stack %>%
#   mutate(., info = map(fit, glance),
#          AICc =  map_dbl(fit, MuMIn::AICc)) %>%
#   select(-fit) %>%
#   unnest(info) %>%
#   select(model_name, sigma, AIC, AICc, BIC, df.residual)


# #Model selection
# #In this instance, we will use AICc score to compare between models. For a model selection approach, the 
# #model with the lowest AICc score is chosen as the model that best supports the data.

# # filter for best model
# best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
# best_model
# #[1] "sharpeschoolhigh"


# # get colour code
# col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# # plot
# ggplot(d_preds, aes(temp, .fitted)) +
#   geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
#   geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
#   geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
#   geom_point(aes(Temp.Cat, ave_rate), Pcom.df_mean) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rd',
#        title = 'Rd across temperatures',
#        subtitle= 'The Sharpe-Schoolfield model is the best model') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 


# ############ Sharpe Schoolfield 1981 model 
# #Schoolfield, R. M., Sharpe, P. J. H., & Magnuson, C. E. (1981). Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of theoretical biology, 88(4), 719-731. https://doi.org/10.1016/0022-5193(81)90246-0 

# #### we want to run the model taking into account means and standard deviations of the PI parameters
# #https://padpadpadpad.github.io/rTPC/articles/weighted_bootstrapping.html
# #https://github.com/fscucchia/PIcurve_TPC_geneExpr_HawaiiBermuda/blob/main/RAnalysis/Hawaii/scripts/weighted_bootstrap_many_curves_DPadfield.R


# # load packages
# library(boot)
# library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
# library(rTPC)
# library(nls.multstart)
# library(broom)
# library(tidyverse)
# library(patchwork)
# library(minpack.lm)


# Pcom.df_mean$temp <- as.numeric(as.character(Pcom.df_mean$Temp.Cat)) #temp needs to be numeric

# Pcom.df_mean <- Pcom.df_mean %>%
#   mutate(Temp.Cat = temp) %>%
#   select(temp, everything(), -Temp.Cat)

# d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
#   mutate(weighted = map(data, ~nls_multstart(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 12),
#                                              data = .x,
#                                              iter = c(4,4,4,4),
#                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') - 1,
#                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'pawar_2018') + 1,
#                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'pawar_2018'),
#                                              supp_errors = 'Y',
#                                              convergence_count = FALSE,
#                                              # include weights here!
#                                              modelweights = 1/sd)))

# # d_fit_Pcom <- nest(Pcom.df_mean, data = c(ave_rate, temp, sd)) %>%
# #   mutate(weighted = map(data, ~nls_multstart(ave_rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 25),
# #                                              data = .x,
# #                                              iter = c(4,4,4,4),
# #                                              start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') - 10,
# #                                              start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981') + 10,
# #                                              #lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              #upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'sharpeschoollow_1981'),
# #                                              supp_errors = 'Y',
# #                                              convergence_count = FALSE,
# #                                              # include weights here!
# #                                              modelweights = 1/sd)))

# # check they work
# d_fit_Pcom$weighted[[1]]



# # get predictions using augment
# newdata_Pcom <- tibble(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))
# d_preds_Pcom <- d_fit_Pcom %>%
#   mutate(., preds = map(weighted, augment, newdata = newdata_Pcom)) %>%
#   select(-weighted) %>%
#   unnest(preds)

# # plot
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pcom) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
#   geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Ik')+
#   #title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 



# ### Residual re-sampling weighted non-linear regression

# #So far we got the best fit to the data. If we want confidence bands around this prediction, we can get those by resampling the data a 
# #number of times. The R package car contains the function Boot() that provides a wrapper for the widely used function boot::boot() that 
# #is tailored to bootstrapping regression models. We refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() 
# #as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# # refit model using nlsLM
# fit_nlsLM_Pcom <- minpack.lm::nlsLM(ave_rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 12),
#                                     data = Pcom.df_mean,
#                                     start = coef(d_fit_Pcom$weighted[[1]]),
#                                     #lower = get_lower_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     #upper = get_upper_lims(Pcom.df_mean$temp, Pcom.df_mean$ave_rate, model_name = 'pawar_2018'),
#                                     weights = 1/sd)


# #Here we will use the weights as a correction to the Pearson residuals of the original model fit. These modified residuals 
# #are then added onto the fitted model predictions, mimicking the methods used for weighted linear regression

# # perform residual bootstrap
# boot2_Pcom <- Boot(fit_nlsLM_Pcom, method = 'residual', R = 100)


# # predict over new data
# boot2_preds_Pcom <- boot2_Pcom$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(Pcom.df_mean$temp), max(Pcom.df_mean$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref,e,eh,topt, tref = 12))

# # calculate bootstrapped confidence intervals
# boot2_conf_preds_Pcom <- boot2_preds_Pcom %>%
#   group_by(temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975),
#             .groups = 'drop')

# # plot bootstrapped CIs - 95% CIs
# p1 <- ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_Pcom, col = 'black') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'black', alpha = 0.3) +
#   geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), Pcom.df_mean) +
#   geom_point(aes(temp, ave_rate), Pcom.df_mean, size = 2, shape = 21, fill = 'green4') +
#   theme_bw(base_size = 10) +
#   theme(legend.position = 'none',
#         strip.text = element_text(hjust = 0),
#         strip.background = element_blank()) +
#   labs(x ='Temperature (ºC)',
#        y = 'Rd',
#        title = 'Photosynthesis rates across temperatures') +
#   geom_hline(aes(yintercept = 0), linetype = 2) 
# #ylim(c(-3, 3.5))

# p1


# ## Calculating confidence intervals of estimated and calculated parameters
# #As for standard non-linear regression, bootstrapping can estimate confidence intervals of the
# #parameters explicitly modelled in the weighted regression. We can compare this approach to profiled 
# #confidence intervals (using confint-MASS) and asymptotic confidence intervals 
# #(using nlstools::confint2()).

# library(nlstools)

# # get parameters of fitted model
# param_Pcom <- broom::tidy(fit_nlsLM_Pcom) %>%
#   select(param = term, estimate)


# # calculate confidence intervals of models
# ci1_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'asymptotic') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'asymptotic')
# ci2_Pcom <- nlstools::confint2(fit_nlsLM_Pcom, method = 'profile')
# #> Waiting for profiling to be done...
# #> Error in prof$getProfile(): number of iterations exceeded maximum of 50
# # profiling method fails
# ci2_Pcom <- mutate(ci1_Pcom, method = 'profile',
#                    conf_lower = NA,
#                    conf_upper = NA)

# # CIs from case resampling
# ci3 <- confint(boot2_Pcom, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')


# # CIs from residual resampling
# ci4_Pcom <- confint(boot2_Pcom, method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# # [1] "All values of t are equal to  0.267937013495708 \n Cannot calculate confidence intervals"


# ci_Pcom <- bind_rows(ci1_Pcom, ci2_Pcom, ci4_Pcom  ) %>%
#   full_join(., param_Pcom, by = "param")
# ci_Ik_Pcom <- bind_rows(ci2_Pcom, ci4_Pcom  ) %>%
#   full_join(., param_Pcom, by = "param")
# #> Joining with `by = join_by(param)`

# ggplot(ci_Ik_Pcom , aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
#   geom_point(size = 1) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
#   labs(title = 'Calculation of confidence intervals for model parameters'
#   )


# # We can also bootstrap confidence intervals for the extra parameters calculated in calc_params().

# extra_params <- calc_params(fit_nlsLM_Pcom) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 200, method = 'case') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'case bootstrap')

# ci_extra_params_Ik_case_Pcom<- left_join(ci_extra_params, extra_params)

# ci_extra_params <- Boot(fit_nlsLM_Pcom, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_Pcom)), R = 100, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')

# ci_extra_params_Ik_residual_Pcom<- left_join(ci_extra_params, extra_params)


# ggplot(ci_extra_params_Ik_residual_Pcom, aes(param, estimate)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   labs(title = 'Calculation of confidence intervals for extra parameters'
#   )

# #filter to only the most relavent and well characterized parameters 
# # All_params <- All_params %>%  
# #   filter(!param=="ctmin") %>% 
# #   filter(!param=="ctmax") %>% 
# #   filter(!param=="skewness") %>% 
# #   filter(!param=="thermal_safety_margin") %>% 
# #   filter(!param=="thermal_tolerance") %>% 
# #   filter(!param=="q10")%>% 
# #   filter(!param=="breadth") 



# #set plot colors 
# cols <- c("Mcap"="green",  
#           "Pacu"="cyan", 
#           "Pcom"="orange") 

# # Ensure Temp.Cat is numeric
# df_wide$Temp.Cat <- as.numeric(as.character(df_wide$Temp.Cat))
# Mcap.df$Temp.Cat <- as.numeric(as.character(Mcap.df$Temp.Cat))
# Pacu.df$Temp.Cat <- as.numeric(as.character(Pacu.df$Temp.Cat))
# Pcom.df$Temp.Cat <- as.numeric(as.character(Pcom.df$Temp.Cat))

# # plot data and model fit for the 3 species
# TPC.plot_Ik <- ggplot(data=df_wide, aes()) + 
#   geom_point(aes(Temp.Cat,Ik, color = "Mcap"), data =  Mcap.df, size = 2, alpha = 0.5) + 
#   geom_point(aes(Temp.Cat,Ik,  color = "Pacu"), data = Pacu.df, size = 2, alpha = 0.5) + 
#   geom_point(aes(Temp.Cat,Ik, color = "Pcom"), data = Pcom.df, size = 2, alpha = 0.5) + 
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Mcap, fill = "green", alpha = 0.3) + 
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pacu, fill = 'cyan', alpha = 0.3) + 
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds_Pcom, fill = 'orange', alpha = 0.3) + 
#   geom_line(aes(temp, .fitted), data = d_preds_Mcap, color = "green") +
#   geom_line(aes(temp, .fitted), data = d_preds_Pacu, color = "cyan") +
#   geom_line(aes(temp, .fitted), data = d_preds_Pcom, color = "orange") +
#   xlim(11.5,40.5)+ 
#   scale_x_continuous(breaks=c(12,14,16,18,20,22,24,26,28,30,32,34,36))+
#   theme_bw(base_size = 12) + 
#   scale_colour_manual(name="Species",values=cols)+ 
#   scale_fill_manual(name = "Species", values = cols) + 
#   theme(legend.position = "top", 
#         panel.border = element_blank(), panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
#   labs(x = 'Temperature (ºC)', 
#        y = expression("Ik"))

# TPC.plot_Ik  



# #### Plot TPC parameters
# # Combine the data frames and add a source column
# ci_extra_params_Ik_case_Mcap <- ci_extra_params_Ik_case_Mcap %>% mutate(source = "Mcap")
# ci_extra_params_Ik_case_Pacu <- ci_extra_params_Ik_case_Pacu %>% mutate(source = "Pacu")
# ci_extra_params_Ik_case_Pcom <- ci_extra_params_Ik_case_Pcom %>% mutate(source = "Pcom")

# ci_extra_params_Ik_residual_Mcap <- ci_extra_params_Ik_residual_Mcap %>% mutate(source = "Mcap")
# ci_extra_params_Ik_residual_Pacu <- ci_extra_params_Ik_residual_Pacu %>% mutate(source = "Pacu")
# ci_extra_params_Ik_residual_Pcom <- ci_extra_params_Ik_residual_Pcom %>% mutate(source = "Pcom")

# # Combine all data frames into one
# ci_extra_params_Ik_case_combined <- bind_rows(ci_extra_params_Ik_case_Mcap , ci_extra_params_Ik_case_Pacu, ci_extra_params_Ik_case_Pcom )

# ci_extra_params_Ik_residual_combined <- bind_rows(ci_extra_params_Ik_residual_Mcap , ci_extra_params_Ik_residual_Pacu, ci_extra_params_Ik_residual_Pcom )

# #set plot colors 
# cols <- c("Mcap"="green",  
#           "Pacu"="cyan", 
#           "Pcom"="orange") 

# # Plot the combined data
# ggplot(ci_extra_params_Ik_residual_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 1, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   #facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for extra parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

# # Plot the combined data
# ggplot(ci_extra_params_Ik_case_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 1, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   #facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for extra parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



# # Add a source column to each data frame
# ci_Ik_Mcap <- ci_Ik_Mcap %>% mutate(source = "Mcap")
# ci_Ik_Pacu <- ci_Ik_Pacu %>% mutate(source = "Pacu")
# ci_Ik_Pcom <- ci_Ik_Pcom %>% mutate(source = "Pcom")

# # Combine all data frames into one
# ci_Ik_combined <- bind_rows(ci_Ik_Mcap, ci_Ik_Pacu, ci_Ik_Pcom )

# #set plot colors 
# cols <- c("Mcap"="green",  
#           "Pacu"="cyan", 
#           "Pcom"="orange") 

# # Plot the combined data
# ggplot(ci_Ik_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for model parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ggplot(ci_Ik_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap(param ~ method, scales = 'free') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Calculation of confidence intervals for model parameters',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# # Filter rows with "residual bootstrap" in the method column
# ci_residual_Ik_combined <- ci_Ik_combined %>%
#   filter(method == "residual bootstrap")

# ci_residual_combined$param = as.factor(ci_residual_combined$param)

# ggplot(ci_residual_Ik_combined, aes(x = param, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   facet_wrap("param", scales = 'free') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(title = 'Confidence intervals for model parameters - residual bootstrap',
#        y = 'Estimate',
#        x = 'Parameter') +
#   theme(legend.position = "top", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




# ######## Statistical analysis of TPC parameters
# ### we have already calculated the confidence intervals for the parameters, so we can use these to calculate the
# ### statistical differences between the species.

# # Filter rows with "residual bootstrap" in the method column
# ci_residual_Ik_combined <- ci_Ik_combined %>%
#   filter(method == "residual bootstrap")

# # View the new table
# print(ci_residual_Ik_combined)

# #Formally test whether confidence intervals overlap using a statistical approach called confidence interval comparison. One common method 
# #is to calculate the Z-score for the difference between two estimates and determine whether the difference is statistically significant.
# #The value 3.92 below is used to approximate the standard error of the confidence interval (CI) bounds. It comes from the fact that a 95% 
# #confidence interval corresponds to approximately 1.96 standard deviations above and below the mean in a normal distribution. Since the 
# #CI width spans both directions (upper and lower bounds), the total width is 2 * 1.96 = 3.92.

# # Load required libraries
# library(dplyr)

# # Function to test overlap of confidence intervals
# test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
#                             estimate2, conf_lower2, conf_upper2) {
#   # Calculate the Z-score for the difference
#   diff <- abs(estimate1 - estimate2)
#   pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
#                       ((conf_upper2 - conf_lower2) / 3.92)^2)
#   z_score <- diff / pooled_se
  
#   # Calculate p-value
#   p_value <- 2 * (1 - pnorm(z_score))
  
#   return(p_value)
# }

# # Apply the function to compare parameters between species
# ci_Ik_comparison <- ci_residual_Ik_combined %>%
#   group_by(param) %>%
#   summarise(
#     p_value_Mcap_vs_Pacu = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
#     ),
#     p_value_Mcap_vs_Pcom = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     ),
#     p_value_Pacu_vs_Pcom = test_ci_overlap(
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     )
#   )

# # View the results
# print(ci_comparison)


# #The 95 % confidence intervals from two subgroups or studies may overlap substantially and yet the test 
# #for difference between them may still produce P < 0.05. 
# #https://pmc.ncbi.nlm.nih.gov/articles/PMC4877414/
# # Confidence Interval Misinterpretation: Confidence intervals are often misinterpreted as a direct test 
# # of significance. However, overlapping confidence intervals do not necessarily mean the difference is not 
# # significant. This is because the overlap does not account for the combined uncertainty of both estimates.

# # Calculate max conf_upper for each parameter
# ci_residual_Ik_combined <- ci_residual_Ik_combined %>%
#   group_by(param) %>%
#   mutate(max_conf_upper = max(conf_upper, na.rm = TRUE)) %>%
#   ungroup()

# # Combine ci_residual_combined with ci_comparison
# ci_plot_Ik_data <- ci_residual_Ik_combined %>%
#   left_join(ci_comparison, by = "param")

# # Prepare the statistical annotations for ci_plot_data
# # Compute y.position values for each parameter and source
# y_positions <- ci_plot_Ik_data %>%
#   group_by(param) %>%
#   summarise(
#     max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
#     .groups = "drop"
#   ) %>%
#   mutate(
#     y_base = max_value  # Add a base offset to the maximum value
#   )

# # Add y.position to annotations with dynamic spacing
# ci_plot_Ik_data <- ci_plot_Ik_data %>%
#   left_join(y_positions, by = "param") %>%
#   group_by(param) %>%
#   mutate(
#     y.position = y_base + (row_number() - 1) * (y_base * 0.4)  # Dynamic offset based on y_base
#   ) %>%
#   ungroup() %>%
#   select(-max_value, -y_base)  # Remove intermediate columns

# # # Filter for significant comparisons only
# # ci_plot_data <- ci_plot_data %>%
# #   filter(p_value_Mcap_vs_Pacu < 0.05 | p_value_Mcap_vs_Pcom < 0.05 | p_value_Pacu_vs_Pcom < 0.05)

# # Add significance levels for each comparison
# ci_plot_Ik_data <- ci_plot_Ik_data %>%
#   mutate(
#     signif_label_Mcap_vs_Pacu = case_when(
#       p_value_Mcap_vs_Pacu < 0.001 ~ "***",
#       p_value_Mcap_vs_Pacu < 0.01 ~ "**",
#       p_value_Mcap_vs_Pacu < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Mcap_vs_Pcom = case_when(
#       p_value_Mcap_vs_Pcom < 0.001 ~ "***",
#       p_value_Mcap_vs_Pcom < 0.01 ~ "**",
#       p_value_Mcap_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Pacu_vs_Pcom = case_when(
#       p_value_Pacu_vs_Pcom < 0.001 ~ "***",
#       p_value_Pacu_vs_Pcom < 0.01 ~ "**",
#       p_value_Pacu_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     )
#   )


# # Compute y.position dynamically for each comparison
# ci_plot_Ik_data <- ci_plot_Ik_data %>%
#   group_by(param) %>%
#   mutate(
#     y_max = max(conf_upper, na.rm = TRUE),
#     y_min = min(conf_lower, na.rm = TRUE),
#     y_range = y_max - y_min,
#     y_position_Mcap_vs_Pacu = y_max + 0.1 * y_range,
#     y_position_Mcap_vs_Pcom = y_max + 0.2 * y_range,
#     y_position_Pacu_vs_Pcom = y_max + 0.3 * y_range
#   ) %>%
#   ungroup()

# # Create the plot
# ci_Ik_data_plot <- ggplot(ci_plot_Ik_data, aes(x = source, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   # Add bars and asterisks for significant comparisons
#   geom_segment(
#     aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
#     data = ci_plot_Ik_data %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 1.5, y = y_position_Mcap_vs_Pacu + 0.000005, label = signif_label_Mcap_vs_Pacu),
#     data = ci_plot_Ik_data %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
#     data = ci_plot_Ik_data %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2, y = y_position_Mcap_vs_Pcom + 0.000005, label = signif_label_Mcap_vs_Pcom),
#     data = ci_plot_Ik_data %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
#     data = ci_plot_Ik_data %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.000005, label = signif_label_Pacu_vs_Pcom),
#     data = ci_plot_Ik_data %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(
#     title = 'Confidence intervals for model parameters with significance',
#     y = 'Estimate',
#     x = 'Source'
#   ) +
#   theme(
#     legend.position = "top",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black")
#   )



# #### Function to test overlap of confidence intervals - EXTRA PARAMETERS
# #ci_extra_params_residual_combined

# test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
#                             estimate2, conf_lower2, conf_upper2) {
#   # Calculate the Z-score for the difference
#   diff <- abs(estimate1 - estimate2)
#   pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
#                       ((conf_upper2 - conf_lower2) / 3.92)^2)
#   z_score <- diff / pooled_se
  
#   # Calculate p-value
#   p_value <- 2 * (1 - pnorm(z_score))
  
#   return(p_value)
# }

# # Apply the function to compare parameters between species
# ci_comparison_extra_Ik_params <- ci_extra_params_Ik_residual_combined %>%
#   group_by(param) %>%
#   summarise(
#     p_value_Mcap_vs_Pacu = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"]
#     ),
#     p_value_Mcap_vs_Pcom = test_ci_overlap(
#       estimate[source == "Mcap"], conf_lower[source == "Mcap"], conf_upper[source == "Mcap"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     ),
#     p_value_Pacu_vs_Pcom = test_ci_overlap(
#       estimate[source == "Pacu"], conf_lower[source == "Pacu"], conf_upper[source == "Pacu"],
#       estimate[source == "Pcom"], conf_lower[source == "Pcom"], conf_upper[source == "Pcom"]
#     )
#   )

# # View the results
# print(ci_comparison_extra_Ik_params)


# # Combine ci_extra_params_case_combined with ci_comparison_extra_params
# ci_plot_data_Ik_extra <- ci_extra_params_Ik_residual_combined %>%
#   left_join(ci_comparison_extra_Ik_params, by = "param")

# # Compute y.position values for each parameter and source
# y_positions <- ci_plot_data_Ik_extra %>%
#   group_by(param) %>%
#   summarise(
#     max_value = max(estimate + conf_upper, na.rm = TRUE),  # Maximum value for the group
#     .groups = "drop"
#   ) %>%
#   mutate(
#     y_base = max_value  # Add a base offset to the maximum value
#   )

# # Add y.position to annotations with dynamic spacing
# ci_plot_data_Ik_extra <- ci_plot_data_Ik_extra %>%
#   left_join(y_positions, by = "param") %>%
#   group_by(param) %>%
#   mutate(
#     y.position = y_base + (row_number() - 1) * (y_base * 0.1)  # Dynamic offset based on y_base
#   ) %>%
#   ungroup() %>%
#   select(-max_value, -y_base)  # Remove intermediate columns



# # Add significance levels for each comparison
# ci_plot_data_Ik_extra <- ci_plot_data_Ik_extra %>%
#   mutate(
#     signif_label_Mcap_vs_Pacu = case_when(
#       p_value_Mcap_vs_Pacu < 0.001 ~ "***",
#       p_value_Mcap_vs_Pacu < 0.01 ~ "**",
#       p_value_Mcap_vs_Pacu < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Mcap_vs_Pcom = case_when(
#       p_value_Mcap_vs_Pcom < 0.001 ~ "***",
#       p_value_Mcap_vs_Pcom < 0.01 ~ "**",
#       p_value_Mcap_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     ),
#     signif_label_Pacu_vs_Pcom = case_when(
#       p_value_Pacu_vs_Pcom < 0.001 ~ "***",
#       p_value_Pacu_vs_Pcom < 0.01 ~ "**",
#       p_value_Pacu_vs_Pcom < 0.05 ~ "*",
#       TRUE ~ NA_character_
#     )
#   )

# # Compute y.position dynamically for each comparison
# ci_plot_data_Ik_extra <- ci_plot_data_Ik_extra %>%
#   group_by(param) %>%
#   mutate(
#     y_max = max(conf_upper, na.rm = TRUE),
#     y_min = min(conf_lower, na.rm = TRUE),
#     y_range = y_max - y_min,
#     y_position_Mcap_vs_Pacu = y_max + 0.2 * y_range,
#     y_position_Mcap_vs_Pcom = y_max + 0.3 * y_range,
#     y_position_Pacu_vs_Pcom = y_max + 0.4 * y_range
#   ) %>%
#   ungroup()

# # Create the plot
# ci_data_Ik_extra_plot <- ggplot(ci_plot_data_Ik_extra, aes(x = source, y = estimate, col = source)) +
#   geom_point(size = 2, position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper), position = position_dodge(width = 0.5)) +
#   # Add bars and asterisks for significant comparisons
#   geom_segment(
#     aes(x = 1, xend = 2, y = y_position_Mcap_vs_Pacu, yend = y_position_Mcap_vs_Pacu),
#     data = ci_plot_data_Ik_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 1.5, y = y_position_Mcap_vs_Pacu + 0.5, label = signif_label_Mcap_vs_Pacu),
#     data = ci_plot_data_Ik_extra %>% filter(!is.na(signif_label_Mcap_vs_Pacu)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 1, xend = 3, y = y_position_Mcap_vs_Pcom, yend = y_position_Mcap_vs_Pcom),
#     data = ci_plot_data_Ik_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2, y = y_position_Mcap_vs_Pcom + 0.5, label = signif_label_Mcap_vs_Pcom),
#     data = ci_plot_data_Ik_extra %>% filter(!is.na(signif_label_Mcap_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_segment(
#     aes(x = 2, xend = 3, y = y_position_Pacu_vs_Pcom, yend = y_position_Pacu_vs_Pcom),
#     data = ci_plot_data_Ik_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   geom_text(
#     aes(x = 2.5, y = y_position_Pacu_vs_Pcom + 0.5, label = signif_label_Pacu_vs_Pcom),
#     data = ci_plot_data_Ik_extra %>% filter(!is.na(signif_label_Pacu_vs_Pcom)),
#     inherit.aes = FALSE, col = "black"
#   ) +
#   theme_bw() +
#   facet_wrap(~param, scales = 'free_y') +
#   scale_x_discrete('') +
#   scale_color_manual(values = cols) +
#   labs(
#     title = 'Confidence intervals for extra parameters with significance',
#     y = 'Estimate',
#     x = 'Source'
#   ) +
#   theme(
#     legend.position = "top",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black")
#   )



########### Plot all TPC estimates ##############

# ci_data_Am_plot
# ci_data_Am_extra_plot
# ci_Rd_data_plot
# ci_data_Rd_extra_plot 
# ci_AQY_data_plot 
# ci_data_AQY_extra_plot
# ci_Ik_data_plot
# ci_data_Ik_extra_plot

# # Load the patchwork package
# library(patchwork)

# # Combine the plots into a single figure
# ci_data_Am_plot <- ci_data_Am_plot + theme(legend.position = "none")
# ci_data_Am_extra_plot <- ci_data_Am_extra_plot + theme(legend.position = "none")
# ci_Rd_data_plot <- ci_Rd_data_plot + theme(legend.position = "none")
# ci_data_Rd_extra_plot <- ci_data_Rd_extra_plot + theme(legend.position = "none")
# ci_AQY_data_plot <- ci_AQY_data_plot + theme(legend.position = "none")
# ci_data_AQY_extra_plot <- ci_data_AQY_extra_plot + theme(legend.position = "none")
# ci_Ik_data_plot <- ci_Ik_data_plot + theme(legend.position = "none")
# #ci_data_Ik_extra_plot <- ci_data_Ik_extra_plot + theme(legend.position = "none")

# combined_plot <- ci_data_Am_plot / ci_data_Am_extra_plot / ci_Rd_data_plot / ci_data_Rd_extra_plot / ci_AQY_data_plot / ci_data_AQY_extra_plot / ci_Ik_data_plot / ci_data_Ik_extra_plot 

# # Display the combined plot
# combined_plot


### Challenge: how to synthesize multiple TPC parameters into a single “thermal tolerance” score for each species, especially 
### when the “best” direction differs by parameter.

# 1. Define the “desirable” direction for each parameter

# 2. Standardize and orient all parameters
# Z-score each parameter across all species:
# z = (value - mean) / sd
# Reverse the sign for parameters where lower is better
# z = -z for those parameters

# 3. Combine into a composite score
# For each parameter, calculate the z-score (with sign flipped as needed).
# Sum or average the z-scores for each species to get a composite “thermal tolerance” score.

library(dplyr)
library(tidyr)
library(fmsb) # for radar plot
library(purrr)


# # Exclude 'breadth' and 'ctmax' from ci_extra_params_AQY_residual_combined
# ci_extra_params_AQY_residual_combined <- ci_extra_params_AQY_residual_combined %>%
#   filter(!param %in% c("breadth", "ctmax"))

# # Exclude 'eh' and 'topt' from ci_residual_AQY_combined
# ci_residual_AQY_combined <- ci_residual_AQY_combined %>%
#   filter(!param %in% c("eh", "topt"))

# Helper function to score a parameter set
score_param_set <- function(df, higher_is_better, lower_is_better) {
  wide <- df %>%
    select(source, param, estimate) %>%
    pivot_wider(names_from = param, values_from = estimate)
  param_cols <- setdiff(colnames(wide), "source")
  wide[param_cols] <- lapply(wide[param_cols], as.numeric)
  for (col in param_cols) {
    if (col %in% higher_is_better) {
      wide[[col]] <- scale(wide[[col]])
    } else if (col %in% lower_is_better) {
      wide[[col]] <- -scale(wide[[col]])
    }
  }
  wide$composite_score <- rowMeans(wide[param_cols], na.rm = TRUE)
  wide
}
# value: The estimate for a given parameter and species (e.g., e for Mcap is 0.737).
# mean: The mean of that parameter's estimates across all species (e.g., mean of e for Mcap, Pacu, and Pcom).
# sd: The standard deviation of that parameter's estimates across all species (e.g., sd of e for Mcap, Pacu, and Pcom).

higher_is_better <- c("eh", "ctmax", "topt", "r_tref", "breadth", "rmax")
lower_is_better <- c("e", "ctmin")
params_to_use <- c("ctmax", "topt", "eh", "e", "r_tref", "breadth", "ctmin", "rmax")

# # --- Am ---
# ci_all_Am_combined <- bind_rows(
#   ci_residual_Am_combined,
#   ci_extra_params_Am_residual_combined
# ) %>%
#   filter(param %in% params_to_use) %>%
#   group_by(source, param) %>%
#   summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop")

# score_Am <- score_param_set(ci_all_Am_combined, higher_is_better, lower_is_better)

# # --- AQY ---
# ci_all_AQY_combined <- bind_rows(
#   ci_residual_AQY_combined,
#   ci_extra_params_AQY_residual_combined
# ) %>%
#   filter(param %in% params_to_use) %>%
#   group_by(source, param) %>%
#   summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop")

# score_AQY <- score_param_set(ci_all_AQY_combined, higher_is_better, lower_is_better)

# # --- Ik ---
# ci_all_Ik_combined <- bind_rows(
#   ci_residual_Ik_combined,
#   ci_extra_params_Ik_residual_combined
# ) %>%
#   filter(param %in% params_to_use) %>%
#   group_by(source, param) %>%
#   summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop")

# score_Ik <- score_param_set(ci_all_Ik_combined, higher_is_better, lower_is_better)

# # --- Combine composite scores ---
# score_summary <- reduce(
#   list(
#     score_Am %>% select(source, Am_score = composite_score),
#     score_AQY %>% select(source, AQY_score = composite_score),
#     score_Ik %>% select(source, Ik_score = composite_score)
#   ),
#   full_join, by = "source"
# ) %>%
#   mutate(overall_score = rowMeans(select(., Am_score, AQY_score, Ik_score), na.rm = TRUE))

# print(score_summary)
# # # A tibble: 3 × 5
# # source Am_score AQY_score Ik_score overall_score
# # <chr>     <dbl>     <dbl>    <dbl>         <dbl>
# #   1 Mcap      0.226     0.472   -0.401        0.0987
# # 2 Pacu     -0.391    -0.110    0.116       -0.128 
# # 3 Pcom      0.165    -0.362    0.285        0.0294


############## repeat for just Am (photosynthesis)

# Define which parameters are higher/lower is better
higher_is_better <- c( "ctmax", "topt", "r_tref", "breadth", "rmax")
lower_is_better <- c("e", "eh", "ctmin")

# Combine all relevant Am parameter tables
ci_all_Am_combined <- bind_rows(
  ci_residual_Am_combined,
  ci_extra_params_Am_residual_combined
)

#Filter to only the parameters you want to use
params_to_use <- c("ctmax", "topt", "eh", "e", "r_tref", "breadth", "ctmin", "rmax")

ci_all_Am_combined <- ci_all_Am_combined %>%
  filter(param %in% params_to_use)
# Summarize to ensure one value per (source, param)
ci_all_Am_combined_unique <- ci_all_Am_combined %>%
  group_by(source, param) %>%
  summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop")

# Score 
score_Am <- score_param_set(ci_all_Am_combined_unique, higher_is_better, lower_is_better)

# The composite_score column in score_Am is your Am-only thermal tolerance score for each species
print(score_Am %>% select(source, Am_composite_score = composite_score))
# # A tibble: 3 × 2
#   source Am_composite_score
#   <chr>               <dbl>
# 1 Mcap               0.0171
# 2 Pacu              -0.459 
# 3 Pcom               0.442 


######## include confidence intervals

# Combine all relevant Am parameter tables
ci_all_Am_combined <- bind_rows(
  ci_residual_Am_combined,
  ci_extra_params_Am_residual_combined
) %>%
  filter(param %in% params_to_use) %>%
  group_by(source, param) %>%
  summarise(
    estimate = mean(estimate, na.rm = TRUE),
    conf_lower = mean(conf_lower, na.rm = TRUE),
    conf_upper = mean(conf_upper, na.rm = TRUE),
    .groups = "drop"
  )

# # Get means and sds for each parameter from the main estimates
# param_stats <- ci_all_Am_combined %>%
#   group_by(param) %>%
#   summarise(
#     mean_est = mean(estimate, na.rm = TRUE),
#     sd_est = sd(estimate, na.rm = TRUE)
#   )
# 
# # Join stats to the main table
# ci_all_Am_combined <- ci_all_Am_combined %>%
#   left_join(param_stats, by = "param") %>%
#   mutate(
#     estimate_z = (estimate - mean_est) / sd_est,
#     conf_lower_z = (conf_lower - mean_est) / sd_est,
#     conf_upper_z = (conf_upper - mean_est) / sd_est
#   )

# Calculate relative CI widths, then apply to z-scored estimates
ci_all_Am_combined <- ci_all_Am_combined %>%
  group_by(param) %>%
  mutate(
    mean_est = mean(estimate, na.rm = TRUE),
    sd_est = sd(estimate, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # Z-score the estimate
    estimate_z = (estimate - mean_est) / sd_est,
    
    # Calculate original CI widths as proportions
    lower_width = estimate - conf_lower,
    upper_width = conf_upper - estimate,
    
    # Scale the widths by the same factor as the estimate
    lower_width_z = lower_width / sd_est,
    upper_width_z = upper_width / sd_est,
    
    # Reconstruct CIs around the z-scored estimate
    conf_lower_z = estimate_z - lower_width_z,
    conf_upper_z = estimate_z + upper_width_z
  )

# Prepare wide tables for z-scored values
wide_estimate <- ci_all_Am_combined %>%
  select(source, param, estimate_z) %>%
  pivot_wider(names_from = param, values_from = estimate_z)

wide_lower <- ci_all_Am_combined %>%
  select(source, param, conf_lower_z) %>%
  pivot_wider(names_from = param, values_from = conf_lower_z)

wide_upper <- ci_all_Am_combined %>%
  select(source, param, conf_upper_z) %>%
  pivot_wider(names_from = param, values_from = conf_upper_z)

# Score function for a wide table
score_param_set_wide_noscale <- function(wide, higher_is_better, lower_is_better) {
  param_cols <- setdiff(colnames(wide), "source")
  for (col in param_cols) {
    if (col %in% lower_is_better) {
      wide[[col]] <- -wide[[col]]
    }
  }
  wide$composite_score <- rowMeans(wide[param_cols], na.rm = TRUE)
  wide
}

#Calculate composite scores for lower and upper bounds
score_main <- score_param_set_wide_noscale(wide_estimate, higher_is_better, lower_is_better)
score_lower <- score_param_set_wide_noscale(wide_lower, higher_is_better, lower_is_better)
score_upper <- score_param_set_wide_noscale(wide_upper, higher_is_better, lower_is_better)

score_main$composite_score_lower <- score_lower$composite_score
score_main$composite_score_upper <- score_upper$composite_score

print(score_main %>% select(source, composite_score, composite_score_lower, composite_score_upper))# # A tibble: 3 × 4
# # A tibble: 3 × 4
# source composite_score composite_score_lower composite_score_upper
# <chr>            <dbl>                 <dbl>                 <dbl>
#   1 Mcap             0.226               -0.170                  1.53 
# 2 Pacu            -0.391               -0.536                  0.182
# 3 Pcom             0.165               -0.0614                 0.915


library(ggplot2)

ggplot(score_main, aes(x = source, y = composite_score, color = source)) +
  geom_pointrange(aes(ymin = composite_score_lower, ymax = composite_score_upper), 
                  position = position_dodge(width = 0.4), size = 1.2) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  theme_bw() +
  labs(
    x = "Species",
    y = "Composite Thermal Tolerance Score",
    title = "Composite Score with Confidence Interval"
  ) +
  scale_color_manual(values = c("Mcap"="green", "Pacu"="cyan", "Pcom"="orange")) +
  scale_x_discrete(expand = c(0.1, 0.1)) +  # <-- This removes extra space at both ends
  theme(
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4), # Optional: adjust margins if needed
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#### statistically compare the composite scores among the three species

# Function to test overlap of confidence intervals
test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
                            estimate2, conf_lower2, conf_upper2) {
  diff <- abs(estimate1 - estimate2)
  pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
                      ((conf_upper2 - conf_lower2) / 3.92)^2)
  z_score <- diff / pooled_se
  p_value <- 2 * (1 - pnorm(z_score))
  return(p_value)
}

# Prepare a summary table for composite scores
score_summary <- score_main %>%
  select(source, composite_score, composite_score_lower, composite_score_upper)

# Compare composite scores between species
composite_comparison <- tibble(
  comparison = c("Mcap_vs_Pacu", "Mcap_vs_Pcom", "Pacu_vs_Pcom"),
  p_value = c(
    test_ci_overlap(
      score_summary$composite_score[score_summary$source == "Mcap"],
      score_summary$composite_score_lower[score_summary$source == "Mcap"],
      score_summary$composite_score_upper[score_summary$source == "Mcap"],
      score_summary$composite_score[score_summary$source == "Pacu"],
      score_summary$composite_score_lower[score_summary$source == "Pacu"],
      score_summary$composite_score_upper[score_summary$source == "Pacu"]
    ),
    test_ci_overlap(
      score_summary$composite_score[score_summary$source == "Mcap"],
      score_summary$composite_score_lower[score_summary$source == "Mcap"],
      score_summary$composite_score_upper[score_summary$source == "Mcap"],
      score_summary$composite_score[score_summary$source == "Pcom"],
      score_summary$composite_score_lower[score_summary$source == "Pcom"],
      score_summary$composite_score_upper[score_summary$source == "Pcom"]
    ),
    test_ci_overlap(
      score_summary$composite_score[score_summary$source == "Pacu"],
      score_summary$composite_score_lower[score_summary$source == "Pacu"],
      score_summary$composite_score_upper[score_summary$source == "Pacu"],
      score_summary$composite_score[score_summary$source == "Pcom"],
      score_summary$composite_score_lower[score_summary$source == "Pcom"],
      score_summary$composite_score_upper[score_summary$source == "Pcom"]
    )
  )
)

print(composite_comparison)
# # A tibble: 3 × 2
# comparison   p_value
# <chr>          <dbl>
#   1 Mcap_vs_Pacu  0.191 
# 2 Mcap_vs_Pcom  0.904 
# 3 Pacu_vs_Pcom  0.0724







############## repeat for just Rd (respiration)

# Define which parameters are higher/lower is better
higher_is_better <- c()
lower_is_better <- c("el", "e", "ctmin")

# Combine all relevant Rd parameter tables
ci_all_Rd_combined <- bind_rows(
  ci_residual_Rd_combined,
  ci_extra_params_Rd_residual_combined
)

#Filter to only the parameters you want to use
params_to_use <- c("el", "e", "ctmin")

ci_all_Rd_combined <- ci_all_Rd_combined %>%
  filter(param %in% params_to_use)
# Summarize to ensure one value per (source, param)
ci_all_Rd_combined_unique <- ci_all_Rd_combined %>%
  group_by(source, param) %>%
  summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop")

# Score 
score_Rd <- score_param_set(ci_all_Rd_combined_unique, higher_is_better, lower_is_better)

# The composite_score column in score_Rd is your Rd-only thermal tolerance score for each species
print(score_Rd %>% select(source, Rd_composite_score = composite_score))
# A tibble: 3 × 2
# source Rd_composite_score
# <chr>               <dbl>
#   1 Mcap              -0.424 
# 2 Pacu               0.498 
# 3 Pcom              -0.0735



######## include confidence intervals

# Combine all relevant Rd parameter tables
ci_all_Rd_combined <- bind_rows(
  ci_residual_Rd_combined,
  ci_extra_params_Rd_residual_combined
) %>%
  filter(param %in% params_to_use) %>%
  group_by(source, param) %>%
  summarise(
    estimate = mean(estimate, na.rm = TRUE),
    conf_lower = mean(conf_lower, na.rm = TRUE),
    conf_upper = mean(conf_upper, na.rm = TRUE),
    .groups = "drop"
  )

# # Get means and sds for each parameter from the main estimates
# param_stats <- ci_all_Rd_combined %>%
#   group_by(param) %>%
#   summarise(
#     mean_est = mean(estimate, na.rm = TRUE),
#     sd_est = sd(estimate, na.rm = TRUE)
#   )
# 
# # Join stats to the main table
# ci_all_Rd_combined <- ci_all_Rd_combined %>%
#   left_join(param_stats, by = "param") %>%
#   mutate(
#     estimate_z = (estimate - mean_est) / sd_est,
#     conf_lower_z = (conf_lower - mean_est) / sd_est,
#     conf_upper_z = (conf_upper - mean_est) / sd_est
#   )

# Calculate relative CI widths, then apply to z-scored estimates
ci_all_Rd_combined <- ci_all_Rd_combined %>%
  group_by(param) %>%
  mutate(
    mean_est = mean(estimate, na.rm = TRUE),
    sd_est = sd(estimate, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # Z-score the estimate
    estimate_z = (estimate - mean_est) / sd_est,
    
    # Calculate original CI widths as proportions
    lower_width = estimate - conf_lower,
    upper_width = conf_upper - estimate,
    
    # Scale the widths by the same factor as the estimate
    lower_width_z = lower_width / sd_est,
    upper_width_z = upper_width / sd_est,
    
    # Reconstruct CIs around the z-scored estimate
    conf_lower_z = estimate_z - lower_width_z,
    conf_upper_z = estimate_z + upper_width_z
  )

# Prepare wide tables for z-scored values
wide_estimate <- ci_all_Rd_combined %>%
  select(source, param, estimate_z) %>%
  pivot_wider(names_from = param, values_from = estimate_z)

wide_lower <- ci_all_Rd_combined %>%
  select(source, param, conf_lower_z) %>%
  pivot_wider(names_from = param, values_from = conf_lower_z)

wide_upper <- ci_all_Rd_combined %>%
  select(source, param, conf_upper_z) %>%
  pivot_wider(names_from = param, values_from = conf_upper_z)

# Score function for a wide table
score_param_set_wide_noscale <- function(wide, higher_is_better, lower_is_better) {
  param_cols <- setdiff(colnames(wide), "source")
  for (col in param_cols) {
    if (col %in% lower_is_better) {
      wide[[col]] <- -wide[[col]]
    }
  }
  wide$composite_score <- rowMeans(wide[param_cols], na.rm = TRUE)
  wide
}

#Calculate composite scores for lower and upper bounds
score_main <- score_param_set_wide_noscale(wide_estimate, higher_is_better, lower_is_better)
score_lower <- score_param_set_wide_noscale(wide_lower, higher_is_better, lower_is_better)
score_upper <- score_param_set_wide_noscale(wide_upper, higher_is_better, lower_is_better)

score_main$composite_score_lower <- score_lower$composite_score
score_main$composite_score_upper <- score_upper$composite_score

print(score_main %>% select(source, composite_score, composite_score_lower, composite_score_upper))# # A tibble: 3 × 4
# # A tibble: 3 × 4
# source composite_score composite_score_lower composite_score_upper
# <chr>            <dbl>                 <dbl>                 <dbl>
#   1 Mcap           -0.424                 -0.254                -0.400
# 2 Pacu            0.498                  1.38                 -0.104
# 3 Pcom           -0.0735                 0.243                -0.668


library(ggplot2)

ggplot(score_main, aes(x = source, y = composite_score, color = source)) +
  geom_pointrange(aes(ymin = composite_score_lower, ymax = composite_score_upper), 
                  position = position_dodge(width = 0.4), size = 1.2) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  theme_bw() +
  labs(
    x = "Species",
    y = "Composite Thermal Tolerance Score",
    title = "Composite Score Rd"
  ) +
  scale_color_manual(values = c("Mcap"="green", "Pacu"="cyan", "Pcom"="orange")) +
  scale_x_discrete(expand = c(0.1, 0.1)) + 
  theme(
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#### statistically compare the composite scores among the three species

# Function to test overlap of confidence intervals
test_ci_overlap <- function(estimate1, conf_lower1, conf_upper1, 
                            estimate2, conf_lower2, conf_upper2) {
  diff <- abs(estimate1 - estimate2)
  pooled_se <- sqrt(((conf_upper1 - conf_lower1) / 3.92)^2 + 
                      ((conf_upper2 - conf_lower2) / 3.92)^2)
  z_score <- diff / pooled_se
  p_value <- 2 * (1 - pnorm(z_score))
  return(p_value)
}

# Prepare a summary table for composite scores
score_summary <- score_main %>%
  select(source, composite_score, composite_score_lower, composite_score_upper)

# Compare composite scores between species
composite_comparison <- tibble(
  comparison = c("Mcap_vs_Pacu", "Mcap_vs_Pcom", "Pacu_vs_Pcom"),
  p_value = c(
    test_ci_overlap(
      score_summary$composite_score[score_summary$source == "Mcap"],
      score_summary$composite_score_lower[score_summary$source == "Mcap"],
      score_summary$composite_score_upper[score_summary$source == "Mcap"],
      score_summary$composite_score[score_summary$source == "Pacu"],
      score_summary$composite_score_lower[score_summary$source == "Pacu"],
      score_summary$composite_score_upper[score_summary$source == "Pacu"]
    ),
    test_ci_overlap(
      score_summary$composite_score[score_summary$source == "Mcap"],
      score_summary$composite_score_lower[score_summary$source == "Mcap"],
      score_summary$composite_score_upper[score_summary$source == "Mcap"],
      score_summary$composite_score[score_summary$source == "Pcom"],
      score_summary$composite_score_lower[score_summary$source == "Pcom"],
      score_summary$composite_score_upper[score_summary$source == "Pcom"]
    ),
    test_ci_overlap(
      score_summary$composite_score[score_summary$source == "Pacu"],
      score_summary$composite_score_lower[score_summary$source == "Pacu"],
      score_summary$composite_score_upper[score_summary$source == "Pacu"],
      score_summary$composite_score[score_summary$source == "Pcom"],
      score_summary$composite_score_lower[score_summary$source == "Pcom"],
      score_summary$composite_score_upper[score_summary$source == "Pcom"]
    )
  )
)

print(composite_comparison)
# # A tibble: 3 × 2
# comparison   p_value
# <chr>          <dbl>
#   1 Mcap_vs_Pacu  0.0151
# 2 Mcap_vs_Pcom  0.136 
# 3 Pacu_vs_Pcom  0.198 




###### Prepare data for radar plot: use z-scored parameters from score - Am
radar_data <- score_Am %>%
  select(-composite_score, -source) # keep only parameter columns

# Add max and min rows required by fmsb
radar_data_plot <- rbind(
  apply(radar_data, 2, max, na.rm = TRUE),
  apply(radar_data, 2, min, na.rm = TRUE),
  radar_data
)

# Set row names for legend
rownames(radar_data_plot) <- c("Max", "Min", score_Am$source)

# Plot
#This will plot each species’ z-score for each Am parameter.
library(fmsb)
radarchart(radar_data_plot, axistype = 1,
           pcol = c("green", "cyan", "orange"),
           plwd = 2, plty = 1,
           title = "Am Parameter Z-scores Radar Plot")
legend("topright", legend = score_Am$source, col = c("green", "cyan", "orange"), lty = 1, lwd = 2)


##### account for higher/lower is better for thermal tolerance

# Use the z-scored, oriented values from score_Am
radar_data <- score_Am %>% select(-composite_score, -source)

# Add max and min rows for fmsb
radar_data_plot <- rbind(
  apply(radar_data, 2, max, na.rm = TRUE),
  apply(radar_data, 2, min, na.rm = TRUE),
  radar_data
)
rownames(radar_data_plot) <- c("Max", "Min", score_Am$source)

# Optionally, rename columns to show direction
colnames(radar_data_plot) <- ifelse(
  colnames(radar_data_plot) %in% higher_is_better,
  paste0(colnames(radar_data_plot), " (+)"),
  paste0(colnames(radar_data_plot), " (−)")
)

# Plot
library(fmsb)
radarchart(radar_data_plot, axistype = 1,
           pcol = c("green", "cyan", "orange"),
           plwd = 2, plty = 1,
           title = "Photosynthetic Thermal Tolerance Radar Plot (Z-scored)")
legend("topright", legend = score_Am$source, col = c("green", "cyan", "orange"), lty = 1, lwd = 2)



###### Prepare data for radar plot: use z-scored parameters from score - Rd
radar_data <- score_Rd %>%
  select(-composite_score, -source) # keep only parameter columns

# Add max and min rows required by fmsb
radar_data_plot <- rbind(
  apply(radar_data, 2, max, na.rm = TRUE),
  apply(radar_data, 2, min, na.rm = TRUE),
  radar_data
)

# Set row names for legend
rownames(radar_data_plot) <- c("Max", "Min", score_Rd$source)

# Plot
#This will plot each species’ z-score for each Rd parameter.
library(fmsb)
radarchart(radar_data_plot, axistype = 1,
           pcol = c("green", "cyan", "orange"),
           plwd = 2, plty = 1,
           title = "Rd Parameter Z-scores Radar Plot")
legend("topright", legend = score_Rd$source, col = c("green", "cyan", "orange"), lty = 1, lwd = 2)


##### account for higher/lower is better for thermal tolerance

### Use the z-scored, oriented values from score - Rd
radar_data <- score_Rd %>% select(-composite_score, -source)

# Add max and min rows for fmsb
radar_data_plot <- rbind(
  apply(radar_data, 2, max, na.rm = TRUE),
  apply(radar_data, 2, min, na.rm = TRUE),
  radar_data
)
rownames(radar_data_plot) <- c("Max", "Min", score_Rd$source)

# Optionally, rename columns to show direction
colnames(radar_data_plot) <- ifelse(
  colnames(radar_data_plot) %in% higher_is_better,
  paste0(colnames(radar_data_plot), " (+)"),
  paste0(colnames(radar_data_plot), " (−)")
)

# Plot
library(fmsb)
radarchart(radar_data_plot, axistype = 1,
           pcol = c("green", "cyan", "orange"),
           plwd = 2, plty = 1,
           title = "Rd Tolerance Radar Plot (Z-scored)")
legend("topright", legend = score_Rd$source, col = c("green", "cyan", "orange"), lty = 1, lwd = 2)




# Make the Rd radar plot have the same axes/parameters as the Am radar plot,
# and show missing parameters for Rd in grey (or as blank spokes).

# List of all Am parameters (from your Am radar plot)
am_params <- colnames(score_Am %>% select(-composite_score, -source))

# Combine Am and Rd parameters, keeping order from Am and adding any extras from Rd
full_params <- unique(c(am_params, setdiff(colnames(score_Rd), c("source", "composite_score"))))

# For Rd: add missing Am parameters as NA, then set NAs to zero
score_Rd_full <- score_Rd %>%
  mutate(!!!setNames(rep(list(NA_real_), length(setdiff(full_params, colnames(score_Rd)))), setdiff(full_params, colnames(score_Rd)))) %>%
  select(source, all_of(full_params))

score_Rd_full[is.na(score_Rd_full)] <- 1

score_Rd_full <- as.data.frame(score_Rd_full)

colnames(score_Rd_full) <- c("source", "breadth", "ctmax", "ctmin", "e", "eh", "r_tref", "rmax", "topt", "el")

# Prepare radar data for Rd
radar_data <- score_Rd_full %>% select(-source)

# Add max and min rows for fmsb
radar_data_plot <- rbind(
  apply(radar_data, 2, function(x) max(x, na.rm = TRUE)),
  apply(radar_data, 2, function(x) min(x, na.rm = TRUE)),
  radar_data
)

# Convert to data.frame for fmsb compatibility
radar_data_plot <- as.data.frame(radar_data_plot)
rownames(radar_data_plot) <- c("Max", "Min", as.character(score_Rd_full$source))


# Plot with all axes present, missing parameters as zero
library(fmsb)
radarchart(radar_data_plot, axistype = 1,
           pcol = c("green", "cyan", "orange"),
           pfcol = c(scales::alpha("green", 0.2), scales::alpha("cyan", 0.2), scales::alpha("orange", 0.2)),
           plwd = 2, plty = 1,
           vlabels = colnames(radar_data_plot),
           title = "Rd Radar Plot (missing parameters set to zero)")
legend("topright", legend = score_Rd_full$source, col = c("green", "cyan", "orange"), lty = 1, lwd = 2)


##### PCA biplot Pmax metrics ###########

library(tidyr)
library(dplyr)
library(tibble)

tpc_pmax_wide <- bind_rows(
  ci_extra_params_Am_residual_combined %>%
    filter(param %in% c("rmax", "ctmin", "ctmax", "breadth")) %>%
    select(species = source, param, estimate),
  ci_residual_Am_combined %>%
    select(species = source, param, estimate)
) %>%
  pivot_wider(names_from = param, values_from = estimate)

# standardize by z-scoring
tpc_pmax_scaled <- tpc_pmax_wide %>%
  column_to_rownames("species") %>%
  scale() %>%
  as.data.frame()       
       
#Run PCA
pca_res <- prcomp(tpc_pmax_scaled, center = TRUE, scale. = TRUE)
       
#Plot the Biplot    
library(ggfortify)
biplot(pca_res, cex = 1.2)

#or
library(factoextra)
fviz_pca_biplot(
  pca_res,
  label = "var", # show variable names
  habillage = rownames(tpc_pmax_scaled), # color by species
  addEllipses = FALSE,
  col.var = "black",
  col.ind = c("green", "cyan", "orange"),
  repel = TRUE,
  title = "PCA Biplot of TPC Metrics (Pmax)"
)


# Bootstrap TPC metrics for each species - take into consideration ci
set.seed(123)
n_boot <- 1000

# Combine both data frames
ci_all <- bind_rows(
  ci_extra_params_Am_residual_combined %>%
    filter(param %in% c("rmax", "ctmin", "ctmax", "breadth")),
  ci_residual_Am_combined
)

# Split by species
boot_samples <- lapply(split(ci_all, ci_all$source), function(df) {
  mat <- mapply(function(mean, lower, upper) {
    rnorm(n_boot, mean, (upper - lower) / (2 * 1.96))
  }, df$estimate, df$conf_lower, df$conf_upper)
  # Set column names to parameter names
  colnames(mat) <- df$param
  as.data.frame(mat)
})
names(boot_samples) <- names(split(ci_all, ci_all$source))

# boot_samples is a list: each element is a n_boot x n_param data.frame for a species
# Add species column and combine
boot_long <- purrr::imap_dfr(boot_samples, ~{
  df <- as.data.frame(.x)
  df$species <- .y
  df
})

# Reorder columns so species is first
boot_long <- boot_long %>% dplyr::select(species, everything())


# Gather into long format for PCA
boot_long <- as.data.frame(boot_long)

# Remove species column, z-score, then add species back
boot_scaled <- boot_long
boot_scaled[,-1] <- scale(boot_long[,-1])

# Run PCA on All Bootstrapped Samples
pca_boot <- prcomp(boot_scaled[,-1], center = TRUE, scale. = TRUE)
# Get PCA scores (PC1, PC2) for each bootstrap sample
scores <- as.data.frame(pca_boot$x)
scores$species <- boot_scaled$species

####### Plot the PCA cloud with loadings
library(ggplot2)

# Get loadings for first two PCs
loadings <- as.data.frame(pca_boot$rotation[,1:2])
loadings$param <- rownames(loadings)

# Plot
ggplot(scores, aes(x = PC1, y = PC2, color = species)) +
  geom_point(alpha = 0.1, size = 1) +  # clouds
  stat_ellipse(aes(fill = species), geom = "polygon", alpha = 0.2, color = NA) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3), 
               arrow = arrow(length = unit(0.2,"cm")), color = "black") +
  geom_text(data = loadings, aes(x = PC1*3.2, y = PC2*3.2, label = param), 
            color = "black", size = 4) +
  scale_color_manual(values = c("Mcap"="green", "Pacu"="cyan", "Pcom"="orange")) +
  scale_fill_manual(values = c("Mcap"="green", "Pacu"="cyan", "Pcom"="orange")) +
  #theme_bw(base_size = 14) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),         # Remove grid
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_blank(),
    legend.position = "top"
  ) +
  labs(title = "PCA Cloud Plot of Bootstrapped TPC Metrics (Pmax)",
       x = "PC1", y = "PC2") +
  theme(legend.position = "top")


# # Prepare PCA scores and loadings
# scores <- as.data.frame(pca_res$x)
# scores$species <- rownames(scores)
# loadings <- as.data.frame(pca_res$rotation[, 1:2])
# loadings$metric <- rownames(loadings)

# Use loadings from the bootstrapped PCA
loadings_boot <- as.data.frame(pca_boot$rotation[, 1:2])
loadings_boot$metric <- rownames(loadings_boot)

arrow_scale <- 2.5

ggplot(scores, aes(x = PC1, y = PC2, color = species)) +
  geom_point(size = 5) +
  # Arrows for loadings from bootstrapped PCA
  geom_segment(
    data = loadings_boot,
    aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
    arrow = arrow(length = unit(0.25, "cm")), color = "black"
  ) +
  # Labels for loadings
  geom_text(
    data = loadings_boot,
    aes(x = PC1 * arrow_scale * 1.1, y = PC2 * arrow_scale * 1.1, label = metric),
    color = "black", size = 5, fontface = "bold"
  ) +
  scale_color_manual(values = c("Mcap" = "green", "Pacu" = "cyan", "Pcom" = "orange")) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA Biplot of TPC Metrics (Pmax) - Bootstrapped Loadings"
  )
       