#####################################################################
# Load libraries, photosynthesis model, plot aesthetics, and datasets
#####################################################################

# Libraries
library(tidyverse)
library(lubridate)

# Load photosynthesis model
source("../functions/photosynthesis_model.R")

# Read in weather datasets
## Hourly climate data
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

## Daily climate data to get rolling means for temps
weather_daily <- read.csv("../data/TT24_daily_weather_summary.csv") %>%
  filter(doy > 99 & doy < 250) %>%
  dplyr::select(doy, tavg10)

# Merge rolling mean with hourly climate
weather <- weather_hourly %>%
  left_join(weather_daily, by = "doy")
head(weather)

# Read in compiled trait dataset
photo_traits <- read.csv("../data/TT24_photo_traits.csv") %>%
  mutate(doy_photo = doy,
         date_photo = date,
         ci = ifelse(ci < 150, 290, ci)) %>%
  group_by(id) %>%
  mutate(doy_max = max(doy)) %>%
  dplyr::select(-doy, -date)
head(photo_traits)

# How many measurements per ID?
n_measurements <- photo_traits %>%
  group_by(id, spp, plot, subplot, gm.trt) %>%
  summarize(n_measurements = length(id)) %>%
  ungroup() %>%
  dplyr::select(id, n_meas = n_measurements)
head(n_measurements)

# How many IDs have >2 measurements?
count(filter(n_measurements, n_meas > 2))

# Unique IDs
id <- unique(photo_traits$id)

# ID with DOY photo and DOY max
photo_traits_doy <- photo_traits %>%
  dplyr::select(id, doy_photo, doy_max)
head(photo_traits_doy)

#####################################################################
# Create data frame that adds id of all individuals in hourly weather 
# dataset, then adds doy of each photo measurement and the final 
# photo measurement. Then, add photo_traits for each individual, 
# filter dataset to include only the most recent photo measurement 
# per doy of weather data. Final prep for spitting data into 
# photosynthesis model
#####################################################################
weather_photo_merged <- weather %>%
  crossing(id = id) %>%
  left_join(photo_traits_doy, by = "id", relationship = "many-to-many") %>%
  left_join(photo_traits, by = c("id", "doy_photo", "doy_max")) %>%
  mutate(time_diff = abs(doy - doy_photo)) %>%
  group_by(doy, id) %>%
  filter(time_diff == min(time_diff)) %>% # Keep the closest measurement
  ungroup() %>%
  dplyr::select(-time_diff, date_time, time:patm, tavg10_weather = tavg10.x, 
                machine, doy_photo, doy_max, id,spp:Tleaf, tavg10_vcmax25 = tavg10.y,
                init_leaf_length_cm:init_leaf_area, -(vp10:vpd10), -anet, -vcmax, -jmax, -rd) %>%
  left_join(n_measurements, by = "id")
head(weather_photo_merged)

#####################################################################
# Run photosynthesis model!
#####################################################################
model_results <- photosynthesis(temp_c = weather_photo_merged$temp_c,
                                temp_c_rolling = weather_photo_merged$tavg10_weather,
                                ci = weather_photo_merged$ci,
                                par = weather_photo_merged$par,
                                patm = weather_photo_merged$patm,
                                vcmax25 = weather_photo_merged$vcmax25,
                                jmax25 = weather_photo_merged$jmax25)
head(model_results)

# Merge photosynthesis model results with `weather_photo_merged`
hourly_model_results <- cbind(weather_photo_merged, model_results)
head(hourly_model_results)

# Remove repeat cols (Vcmax, Jmax, Ci)
hourly_model_results <- hourly_model_results[, -c(36, 38, 40, 43:44)]

names(hourly_model_results)

#####################################################################
# Convert A, Rd, and Anet from umol/m2/s to g C/hr, filter to only
# include individuals that have greater than 2 gas exchange measurements
#####################################################################
hourly_model_results <- hourly_model_results %>%
  filter(n_meas > 2) %>%
  mutate(hourly_C_assim = a * 3600 / 1e6 * 12, # gC per hour
         hourly_C_resp = rd * 3600 / 1e6 * 12,
         hourly_netC_assim = anet * 3600 / 1e6 * 12) %>%
  mutate(hourly_C_assim_tla = hourly_C_assim * (total_leaf_area_cm2 / 10000),
         hourly_C_resp_tla = hourly_C_resp * (total_leaf_area_cm2 / 10000),
         hourly_netC_assim_tla = hourly_netC_assim * (total_leaf_area_cm2 / 10000))

#####################################################################
# Calculate total daily C estimates
#####################################################################
# Calculate daily carbon budgets
daily_model_results <- hourly_model_results %>%
  group_by(doy, id, spp, plot, subplot, gm.trt, n_meas) %>%
  summarize(daily_c_assim = sum(hourly_C_assim, na.rm = TRUE),
            daily_c_resp = sum(hourly_C_resp, na.rm = TRUE),
            daily_netc_assim = sum(hourly_netC_assim, na.rm = TRUE),
            
            daily_c_assim_tla = sum(hourly_C_assim_tla, na.rm = TRUE),
            daily_c_resp_tla = sum(hourly_C_resp_tla, na.rm = TRUE),
            daily_netc_assim_tla = sum(hourly_netC_assim_tla, na.rm = TRUE))

daily_model_results <- daily_model_results %>%
  arrange(id, doy) %>%
  group_by(id, spp, plot, subplot, gm.trt, n_meas) %>%
  mutate(cumul_c_assim = cumsum(daily_c_assim),
         cumul_c_resp = cumsum(daily_c_resp),
         cumul_netc_assim = cumsum(daily_netc_assim),
         
         cumul_c_assim_tla = cumsum(daily_c_assim_tla),
         cumul_c_resp_tla = cumsum(daily_c_resp_tla),
         cumul_netc_assim_tla = cumsum(daily_netc_assim_tla))

# Calculate total seasonal C budget
total_c_budget <- daily_model_results %>%
  group_by(id, spp, plot, subplot, gm.trt, n_meas) %>%
  summarize(total_c_assim = sum(daily_c_assim),
            total_c_resp = sum(daily_c_resp),
            total_netc_assim = sum(daily_netc_assim),
            
            total_c_assim_tla = sum(daily_c_assim_tla),
            total_c_resp_tla = sum(daily_c_resp_tla),
            total_netc_assim_tla = sum(daily_netc_assim_tla)
            )

# Merge total seasonal C budget with daily model results
full_model_results <- daily_model_results %>%
  full_join(total_c_budget) %>%
  mutate(prop_c_assim = cumul_c_assim / total_c_assim)

# Write total C budget model results to .csv
total_c_budget %>%
  mutate(sample_year = 2024) %>%
  write.csv("../data/c_budget/TT24_seasonal_c_budget.csv", 
            row.names = F)
