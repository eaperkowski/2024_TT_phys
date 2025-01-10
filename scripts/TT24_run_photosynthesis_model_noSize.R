# Libraries
library(tidyverse)

# Load photosynthesis model
source("../functions/photosynthesis_model.R")

# Load hourly climate data
weather_hourly <- read.csv("../data/TT24_hourly_weather_data.csv") %>%
  filter(doy > 99 & doy < 250) %>%
  mutate(atm_pressure_pa = atm_pressure_kpa * 1000) %>%
  dplyr::select(date_time = date, 
                doy, 
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

# Load bootstrapped model predictions (n = 1000 per day)
model_params <- read.csv("../data/c_budget/TT24_budget_params_without_size.csv") %>%
  filter(doy > 100 & doy < 250)

# Prep hourly climate dataset
weather_with_params <- weather %>%
  full_join(model_params, by = "doy")
head(weather_with_params)

# Create weather_with params such that it just has descriptors
weather_with_params_descript <- weather_with_params %>%
  dplyr::select(date_time:doy, gm.trt:Ci_est)

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


