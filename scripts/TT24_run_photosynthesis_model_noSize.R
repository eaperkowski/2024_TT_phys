# Libraries
library(tidyverse)
library(lubridate)

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

# Next, lets set up the bootstrap with 5000 iterations and
# create a placeholder vector for values to be loaded into

calculate_daily_netC_gain <- function(df) {
  hourly_bootstrap <- df %>%
    group_by(doy, hour) %>%
    summarize(resampled_value = sample(hourly_netC_gain, 1), .groups = "drop")
  
  daily_netC_gain <- hourly_bootstrap %>%
    group_by(doy) %>%
    summarize(daily_netC_gain = sum(resampled_value), .groups = "drop")
  
  return(list(hourly_input = hourly_bootstrap, daily_netC_gain = daily_netC_gain))
}


# Next, lets set up the bootstrap with 5000 iterations and
# create a placeholder vector for values to be loaded into
# Number of iterations
n_iterations <- 5000
hourly_results_list <- vector("list", n_iterations)
daily_results_list <- vector("list", n_iterations)

# Alright, lets iterate the bootstrap
for (i in 1:n_iterations) {
  # Run the modified function
  results <- calculate_daily_netC_gain(full_model_tri_weeded)
  
  # Store hourly and daily results, adding iteration number
  hourly_results_list[[i]] <- results$hourly_input %>%
    mutate(iteration = i)
  daily_results_list[[i]] <- results$daily_netC_gain %>%
    mutate(iteration = i)
}

# Combine results into single data frames
hourly_all_iterations <- bind_rows(hourly_results_list) %>%
  arrange(doy, hour, iteration)
daily_all_iterations <- bind_rows(daily_results_list) %>%
  arrange(doy)

daily_netC_gain_summary <- daily_all_iterations %>%
  group_by(doy) %>%
  summarize(n = length(iteration),
            
            daily_netC_gain_mean = mean(daily_netC_gain),
            daily_netC_gain_sd = sd(daily_netC_gain),
            daily_netC_gain_se = daily_netC_gain_sd / sqrt(n)) %>%
  mutate(daily_netC_gain_LCI = daily_netC_gain_mean - (1.96 * daily_netC_gain_se),
         daily_netC_gain_UCI = daily_netC_gain_mean + (1.96 * daily_netC_gain_se),
         cum_netC_gain = cumsum(daily_netC_gain_mean))
         

ggplot(daily_netC_gain_summary, aes(x = doy, y = daily_netC_gain_mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = daily_netC_gain_LCI,
                    ymax = daily_netC_gain_UCI)) +
  labs(x = "Day of year",
       y = expression("Daily net C gain (g CO"["2"]*" m"^"-2"*" s"^"-1"*")")) +
  geom_line() +
  theme_classic(base_size = 18)

ggplot(daily_netC_gain_summary, aes(x = doy, y = cum_netC_gain)) +
  geom_point(size = 2) +
  labs(x = "Day of year",
       y = expression("Cumulative net C gain (g CO"["2"]*" m"^"-2"*" s"^"-1"*")")) +
  geom_line() +
  theme_classic(base_size = 18)



# Calculate summary statistics for hourly ac, aj, a, rd, and anet
model_daily_summary <- full_model_results %>%
  group_by(date_time, doy, time, gm.trt, spp) %>%
  summarize(n = length(boot_iter),
            
            Anet_mean = mean(anet, na.rm = TRUE),
            Anet_sd = sd(anet, na.rm = TRUE),
            Anet_se = Anet_sd / sqrt(n),
            
            Vcmax_mean = mean(vcmax, na.rm = TRUE),
            Vcmax_sd = sd(vcmax, na.rm = TRUE),
            Vcmax_se = Vcmax_sd / sqrt(n),
            
            Vcmax25_mean = mean(Vcmax25_est, na.rm = TRUE), 
            Vcmax25_sd = sd(Vcmax25_est),
            Vcmax25_se = Vcmax25_sd / sqrt(n),
            
            Jmax_mean = mean(jmax, na.rm = TRUE), 
            Jmax_sd = sd(jmax),
            Jmax_se = Jmax_sd / sqrt(n),
            
            Jmax25_mean = mean(Jmax25_est, na.rm = TRUE), 
            Jmax25_sd = sd(Jmax25_est),
            Jmax25_se = Jmax25_sd / sqrt(n)
            )



ggplot(data = subset(model_daily_summary, doy == 200 & spp == "Tri" & gm.trt == "ambient"), 
       aes(x = time)) +
  geom_point(aes(y = Anet_mean, color = gm.trt, shape = spp),
             size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~spp*gm.trt)


