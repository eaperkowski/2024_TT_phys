# Libraries
library(tidyverse)
library(lubridate)
library(ggpubr)

## Color palettes
gm.colors <- c("#00B2BE", "#F1B700")

# Load photosynthesis model
source("../functions/photosynthesis_model.R")

# Load hourly climate data
weather_hourly <- read.csv("../data/TT24_hourly_weather_data.csv") %>%
  filter(doy > 99 & doy < 250) %>%
  mutate(atm_pressure_pa = atm_pressure_kpa * 1000,
         hour = hour(date)) %>%
  dplyr::select(date_time = date, 
                doy,
                hour,
                par = par_umol_m2_s,
                temp_c = air_temperature_c,
                patm = atm_pressure_pa) %>%
  separate(date_time, c("date", "time"), sep = " ", remove = FALSE)

# Load daily climate data to get rolling means for temps
weather_daily <- read.csv("../data/TT24_daily_weather_summary.csv") %>%
  filter(doy > 99 & doy < 250) %>%
  dplyr::select(doy, tavg10)

# Merge rolling mean with hourly climate
weather <- weather_hourly %>%
  left_join(weather_daily, by = "doy")
head(weather)

# Load bootstrapped model predictions (n = 1000 per day)
model_params <- read.csv("../data/c_budget/TT24_budget_params_without_size.csv") %>%
  filter(doy > 100 & doy < 250)

# Prep hourly climate dataset
weather_with_params <- weather %>%
  full_join(model_params, by = "doy")
head(weather_with_params)

# Create weather_with params such that it just has descriptors
weather_with_params_descript <- weather_with_params %>%
  dplyr::select(date_time:hour, gm.trt:Ci_est)

# Test photosynthesis model
photosynthesis()

# Run photosynthesis model!
model_results <- photosynthesis(temp_c = weather_with_params$temp_c,
                                temp_c_rolling = weather_with_params$tavg10,
                                ci = weather_with_params$Ci_est,
                                par = weather_with_params$par,
                                patm = weather_with_params$patm,
                                vcmax25 = weather_with_params$Vcmax25_est,
                                jmax25 = weather_with_params$Jmax25_est)

# Merge photosynthesis model results with `weather_with_params`
full_model_results <- cbind(weather_with_params_descript, model_results)
head(full_model_results)

# Convert A, Rd, and Anet from umol/m2/s to g CO2/ms/s
# (umol/m2/s * (3600 sec/1 hr) * (1 mol CO2 / 10^6 umol CO2) * (44 g CO2 / 1 mol CO2))
full_model_results$hourly_C_gain <- full_model_results$a * 3600 / 1e6 * 44 # g CO2/m2/hr
full_model_results$hourly_C_lost <- full_model_results$rd * 3600 / 1e6 * 44 # g CO2 mol/m2/hr
full_model_results$hourly_netC_gain <- full_model_results$anet * 3600 / 1e6 * 44 # g CO2/m2/hr
head(full_model_results)

# Ok, so we have 1000 estimates for hourly C gain. We now need to 
# scale these values up to daily C gain. I think the best way to 
# do this is to randomly sample (with replacement) one hourly_netC_gain
# value for each hour of each day, then repeat this process 5000 times
# (i.e., another bootstrap)

# First, lets create a data frame for each gm.trt and spp for easy
# load into for loop
full_model_tri_ambient <- subset(full_model_results, spp == "Tri" & gm.trt == "ambient")
full_model_tri_weeded <- subset(full_model_results, spp == "Tri" & gm.trt == "weeded")
full_model_mai_ambient <- subset(full_model_results, spp == "Mai" & gm.trt == "ambient")
full_model_mai_weeded <- subset(full_model_results, spp == "Mai" & gm.trt == "weeded")

# Let's set up a function that randomly samples one value from
# photosynthesis model results
calculate_daily_c_budget <- function(df) {
  
  # First, sample one value for net C gain, C gain, and C lost for each
  # hour of each doy
  hourly_bootstrap <- df %>%
    group_by(doy, hour) %>%
    summarize(hourly_cGain_resample = sample(hourly_C_gain, 1),
              hourly_cLost_resample = sample(hourly_C_lost, 1),
              hourly_netC_resample = sample(hourly_netC_gain, 1), .groups = "drop")
  
  # Then, use samples generated in hourly_bootstrap to estimate total
  # daily C gain. This is done as a sum of the single resampled value
  # for each hour in a given day
  daily_results <- hourly_bootstrap %>%
    group_by(doy) %>%
    summarize(daily_c_gain = sum(hourly_cGain_resample),
              daily_c_lost = sum(hourly_cLost_resample),
              daily_netC_gain = sum(hourly_netC_resample), .groups = "drop")
  
  # Return the resampled hourly values and the daily budget results in a list
  return(list(hourly_input = hourly_bootstrap, daily_budget = daily_results))
}

# Create vector for number of bootstrap iterations
n_iterations <- 1000

# Alright, it's time. Let's iterate the bootstrap 5000 times and
# start with Trillium weeded
hourly_results_tri_weeded <- vector("list", n_iterations)
daily_results_tri_weeded <- vector("list", n_iterations)

for (i in 1:n_iterations) {
  # Run the modified function
  results <- calculate_daily_c_budget(full_model_tri_weeded)
  
  # Store hourly and daily results, adding iteration number
  hourly_results_tri_weeded[[i]] <- results$hourly_input %>%
    mutate(iteration = i)
  daily_results_tri_weeded[[i]] <- results$daily_budget %>%
    mutate(iteration = i)
}

# Combine bootstrapped results into single data frame- one for
# hourly iterations and one for daily iterations. The hourly
# iteration df is mainly to store the values used to calculate
# the daily means for reproducibility
hourly_results_tri_weeded <- bind_rows(hourly_results_tri_weeded) %>%
  arrange(doy, hour, iteration)
daily_results_tri_weeded <- bind_rows(daily_results_tri_weeded) %>%
  arrange(doy)

daily_c_budget_tri_weeded <- daily_results_tri_weeded %>%
  filter(doy > 104 & doy < 167) %>%
  group_by(doy) %>%
  summarize(n = length(iteration),
            
            daily_netC_gain_mean = mean(daily_netC_gain),
            daily_netC_gain_sd = sd(daily_netC_gain),
            daily_netC_gain_se = daily_netC_gain_sd / sqrt(n),
            
            daily_c_gain_mean = mean(daily_c_gain),
            daily_c_gain_sd = sd(daily_c_gain),
            daily_c_gain_se = daily_c_gain_sd / sqrt(n),
            
            daily_c_lost_mean = mean(daily_c_lost),
            daily_c_lost_sd = sd(daily_c_lost),
            daily_c_lost_se = daily_c_lost_sd / sqrt(n)) %>%
  mutate(daily_netC_gain_LCI = daily_netC_gain_mean - (1.96 * daily_netC_gain_se),
         daily_netC_gain_UCI = daily_netC_gain_mean + (1.96 * daily_netC_gain_se),
         
         daily_C_gain_LCI = daily_c_gain_mean - (1.96 * daily_c_gain_se),
         daily_C_gain_UCI = daily_c_gain_mean + (1.96 * daily_c_gain_se),
         
         daily_C_lost_LCI = daily_c_lost_mean - (1.96 * daily_c_lost_se),
         daily_C_lost_UCI = daily_c_lost_mean + (1.96 * daily_c_lost_se),
         
         cumul_netC_gain = cumsum(daily_netC_gain_mean),
         cumul_c_gain = cumsum(daily_c_gain_mean),
         cumul_c_lost = cumsum(daily_c_lost_mean),
         
         cumul_netC_gain_fraction = cumul_netC_gain / max(cumul_netC_gain),
         cumul_netC_gain_percent = cumsum(cumul_netC_gain_fraction),
         
         gm.trt = "weeded") %>%
  add_row(doy = 103, cumul_netC_gain = 0, cumul_c_gain = 0, cumul_c_lost = 0,
          cumul_netC_gain_fraction = 0, gm.trt = "weeded")

head(daily_c_budget_tri_weeded)
## Everything seems to work fine! Repeat process above with the other unique
## spp-by-gm.trt combinations

# Alright, it's time. Let's iterate the bootstrap 5000 times and
# start with Trillium weeded
hourly_results_tri_ambient <- vector("list", n_iterations)
daily_results_tri_ambient <- vector("list", n_iterations)

for (i in 1:n_iterations) {
  # Run the modified function
  results <- calculate_daily_c_budget(full_model_tri_ambient)
  
  # Store hourly and daily results, adding iteration number
  hourly_results_tri_ambient[[i]] <- results$hourly_input %>%
    mutate(iteration = i)
  daily_results_tri_ambient[[i]] <- results$daily_budget %>%
    mutate(iteration = i)
}

# Combine bootstrapped results into single data frame- one for
# hourly iterations and one for daily iterations. The hourly
# iteration df is mainly to store the values used to calculate
# the daily means for reproducibility
hourly_results_tri_ambient <- bind_rows(hourly_results_tri_ambient) %>%
  arrange(doy, hour, iteration)
daily_results_tri_ambient <- bind_rows(daily_results_tri_ambient) %>%
  arrange(doy)

daily_c_budget_tri_ambient <- daily_results_tri_ambient %>%
  filter(doy > 104 & doy < 167) %>%
  group_by(doy) %>%
  summarize(n = length(iteration),
            
            daily_netC_gain_mean = mean(daily_netC_gain),
            daily_netC_gain_sd = sd(daily_netC_gain),
            daily_netC_gain_se = daily_netC_gain_sd / sqrt(n),
            
            daily_c_gain_mean = mean(daily_c_gain),
            daily_c_gain_sd = sd(daily_c_gain),
            daily_c_gain_se = daily_c_gain_sd / sqrt(n),
            
            daily_c_lost_mean = mean(daily_c_lost),
            daily_c_lost_sd = sd(daily_c_lost),
            daily_c_lost_se = daily_c_lost_sd / sqrt(n)) %>%
  mutate(daily_netC_gain_LCI = daily_netC_gain_mean - (1.96 * daily_netC_gain_se),
         daily_netC_gain_UCI = daily_netC_gain_mean + (1.96 * daily_netC_gain_se),
         
         daily_C_gain_LCI = daily_c_gain_mean - (1.96 * daily_c_gain_se),
         daily_C_gain_UCI = daily_c_gain_mean + (1.96 * daily_c_gain_se),
         
         daily_C_lost_LCI = daily_c_lost_mean - (1.96 * daily_c_lost_se),
         daily_C_lost_UCI = daily_c_lost_mean + (1.96 * daily_c_lost_se),
         
         cumul_netC_gain = cumsum(daily_netC_gain_mean),
         cumul_c_gain = cumsum(daily_c_gain_mean),
         cumul_c_lost = cumsum(daily_c_lost_mean),
         
         cumul_netC_gain_fraction = cumul_netC_gain / max(cumul_netC_gain),
         
         gm.trt = "ambient") %>%
  add_row(doy = 103, cumul_netC_gain = 0, cumul_c_gain = 0, cumul_c_lost = 0,
          cumul_netC_gain_fraction = 0, gm.trt = "ambient")

# Combine weeded and ambient model results
daily_c_budget_tri_combined <- daily_c_budget_tri_weeded %>%
  full_join(daily_c_budget_tri_ambient) %>%
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")),
         spp = "Tri") %>%
  dplyr::select(doy, gm.trt, n, everything())
head(daily_c_budget_tri_combined)

# Add code for facet labels
facet_lab_tri <- "Trillium spp."
names(facet_lab_tri) <- "Tri"
  
# Plots
tri_netC_gain_plot <- ggplot(data = daily_c_budget_tri_combined,
       aes(x = doy, y = daily_netC_gain_mean)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  geom_line(aes(color = gm.trt)) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_x_continuous(limits = c(100, 170), breaks = seq(100, 170, 10)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C gain (g CO"["2"]*" m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

tri_cumul_gain_plot <- ggplot(data = daily_c_budget_tri_combined,
       aes(x = doy, y = cumul_netC_gain_fraction)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_line(aes(color = gm.trt)) +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(limits = c(100, 170), breaks = seq(100, 170, 10)) +
  labs(x = "Day of year",
       y = "Cumulative C gain (%)",
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# Iterate the bootstrap 5000 times for Mai weeded
hourly_results_mai_weeded <- vector("list", n_iterations)
daily_results_mai_weeded <- vector("list", n_iterations)

for (i in 1:n_iterations) {
  # Run the modified function
  results <- calculate_daily_c_budget(full_model_mai_weeded)
  
  # Store hourly and daily results, adding iteration number
  hourly_results_mai_weeded[[i]] <- results$hourly_input %>%
    mutate(iteration = i)
  daily_results_mai_weeded[[i]] <- results$daily_budget %>%
    mutate(iteration = i)
}

# Combine bootstrapped results into single data frame- one for
# hourly iterations and one for daily iterations. The hourly
# iteration df is mainly to store the values used to calculate
# the daily means for reproducibility
hourly_results_mai_weeded <- bind_rows(hourly_results_mai_weeded) %>%
  arrange(doy, hour, iteration)
daily_results_mai_weeded <- bind_rows(daily_results_mai_weeded) %>%
  arrange(doy)

daily_c_budget_mai_weeded <- daily_results_mai_weeded %>%
  filter(doy > 119) %>%
  group_by(doy) %>%
  summarize(n = length(iteration),
            
            daily_netC_gain_mean = mean(daily_netC_gain),
            daily_netC_gain_sd = sd(daily_netC_gain),
            daily_netC_gain_se = daily_netC_gain_sd / sqrt(n),
            
            daily_c_gain_mean = mean(daily_c_gain),
            daily_c_gain_sd = sd(daily_c_gain),
            daily_c_gain_se = daily_c_gain_sd / sqrt(n),
            
            daily_c_lost_mean = mean(daily_c_lost),
            daily_c_lost_sd = sd(daily_c_lost),
            daily_c_lost_se = daily_c_lost_sd / sqrt(n)) %>%
  mutate(daily_netC_gain_LCI = daily_netC_gain_mean - (1.96 * daily_netC_gain_se),
         daily_netC_gain_UCI = daily_netC_gain_mean + (1.96 * daily_netC_gain_se),
         
         daily_C_gain_LCI = daily_c_gain_mean - (1.96 * daily_c_gain_se),
         daily_C_gain_UCI = daily_c_gain_mean + (1.96 * daily_c_gain_se),
         
         daily_C_lost_LCI = daily_c_lost_mean - (1.96 * daily_c_lost_se),
         daily_C_lost_UCI = daily_c_lost_mean + (1.96 * daily_c_lost_se),
         
         cumul_netC_gain = cumsum(daily_netC_gain_mean),
         cumul_c_gain = cumsum(daily_c_gain_mean),
         cumul_c_lost = cumsum(daily_c_lost_mean),
         
         cumul_netC_gain_fraction = cumul_netC_gain / max(cumul_netC_gain),
         
         gm.trt = "weeded") %>%
  add_row(doy = 118, cumul_netC_gain = 0, cumul_c_gain = 0, cumul_c_lost = 0,
          cumul_netC_gain_fraction = 0, gm.trt = "weeded")

# Iterate the bootstrap 5000 times for Mai ambient
hourly_results_mai_ambient <- vector("list", n_iterations)
daily_results_mai_ambient <- vector("list", n_iterations)

for (i in 1:n_iterations) {
  # Run the modified function
  results <- calculate_daily_c_budget(full_model_mai_ambient)
  
  # Store hourly and daily results, adding iteration number
  hourly_results_mai_ambient[[i]] <- results$hourly_input %>%
    mutate(iteration = i)
  daily_results_mai_ambient[[i]] <- results$daily_budget %>%
    mutate(iteration = i)
}

# Combine bootstrapped results into single data frame- one for
# hourly iterations and one for daily iterations. The hourly
# iteration df is mainly to store the values used to calculate
# the daily means for reproducibility
hourly_results_mai_ambient <- bind_rows(hourly_results_mai_ambient) %>%
  arrange(doy, hour, iteration)
daily_results_mai_ambient <- bind_rows(daily_results_mai_ambient) %>%
  arrange(doy)

daily_c_budget_mai_ambient <- daily_results_mai_ambient %>%
  filter(doy > 119) %>%
  group_by(doy) %>%
  summarize(n = length(iteration),
            
            daily_netC_gain_mean = mean(daily_netC_gain),
            daily_netC_gain_sd = sd(daily_netC_gain),
            daily_netC_gain_se = daily_netC_gain_sd / sqrt(n),
            
            daily_c_gain_mean = mean(daily_c_gain),
            daily_c_gain_sd = sd(daily_c_gain),
            daily_c_gain_se = daily_c_gain_sd / sqrt(n),
            
            daily_c_lost_mean = mean(daily_c_lost),
            daily_c_lost_sd = sd(daily_c_lost),
            daily_c_lost_se = daily_c_lost_sd / sqrt(n)) %>%
  mutate(daily_netC_gain_LCI = daily_netC_gain_mean - (1.96 * daily_netC_gain_se),
         daily_netC_gain_UCI = daily_netC_gain_mean + (1.96 * daily_netC_gain_se),
         
         daily_C_gain_LCI = daily_c_gain_mean - (1.96 * daily_c_gain_se),
         daily_C_gain_UCI = daily_c_gain_mean + (1.96 * daily_c_gain_se),
         
         daily_C_lost_LCI = daily_c_lost_mean - (1.96 * daily_c_lost_se),
         daily_C_lost_UCI = daily_c_lost_mean + (1.96 * daily_c_lost_se),
         
         cumul_netC_gain = cumsum(daily_netC_gain_mean),
         cumul_c_gain = cumsum(daily_c_gain_mean),
         cumul_c_lost = cumsum(daily_c_lost_mean),
         
         cumul_netC_gain_fraction = cumul_netC_gain / max(cumul_netC_gain),
         
         gm.trt = "ambient") %>%
  add_row(doy = 118, cumul_netC_gain = 0, cumul_c_gain = 0, cumul_c_lost = 0,
          cumul_netC_gain_fraction = 0, gm.trt = "ambient")

# Combine model results from weeded and ambient plots
daily_c_budget_mai_combined <- daily_c_budget_mai_weeded %>%
  full_join(daily_c_budget_mai_ambient) %>%
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")),
         spp = "Mai") %>%
  dplyr::select(doy, gm.trt, n, everything())
head(daily_c_budget_mai_combined)

# Add code for facet labels
facet_lab_tri <- "M. racemosum"
names(facet_lab_tri) <- "Mai"


mai_netC_gain_plot <- ggplot(data = daily_c_budget_mai_combined,
       aes(x = doy, y = daily_netC_gain_mean)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_line(aes(color = gm.trt)) +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C gain (g CO"["2"]*" m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

mai_cumul_gain_plot <- ggplot(data = daily_c_budget_mai_combined,
       aes(x = doy, y = cumul_netC_gain_fraction)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_line(aes(color = gm.trt)) +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = "Cumulative C gain (%)",
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

## Make figure for Summit
png("../drafts/figs/TT24_cBudget.png", width = 10, height = 10,
    units = "in", res = 600)
ggarrange(tri_netC_gain_plot, mai_netC_gain_plot,
          tri_cumul_gain_plot, mai_cumul_gain_plot,
          nrow = 2, ncol = 2, common.legend = TRUE,
          legend = "bottom", align = "hv",
          labels = c("(a)", "(b)", "(c)", "(d)"),
          font.label = list(size = 18, hjust = 0))
dev.off()




