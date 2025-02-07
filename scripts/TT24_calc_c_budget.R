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
  mutate(doy_photo = doy) %>%
  group_by(id) %>%
  mutate(doy_max = max(doy)) %>%
  dplyr::select(-doy)
head(photo_traits)

# How many measurements per ID?
n_measurements <- photo_traits %>%
  group_by(id, spp, plot, subplot, gm.trt) %>%
  summarize(n_measurements = length(id)) %>%
  arrange(n_measurements)
head(n_measurements)

id <- unique(photo_traits$id)

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
    dplyr::select(-time_diff, date_time, weather_date = date.x, time:patm, tavg10_weather = tavg10.x, 
           machine, photo_date = date.y, doy_photo, doy_max, id,spp:Tleaf, tavg10_vcmax25 = tavg10.y,
           init_leaf_length_cm:init_leaf_area, -(vp10:vpd10), -anet, -vcmax, -jmax, -rd) 
head(weather_photo_merged)

# check to make sure different leaf areas are included for sample ID
test <- filter(weather_photo_merged, id == "1021")

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
hourly_model_results <- hourly_model_results[, -c(35:37, 39)]

#####################################################################
# Convert A, Rd, and Anet from umol/m2/s to g C/ms/s, then to
# g C / day
#####################################################################
hourly_model_results <- hourly_model_results %>%
  mutate(hourly_C_assim = a * 3600 / 1e6 * 12,
         hourly_C_resp = rd * 3600 / 1e6 * 12,
         hourly_netC_assim = anet * 3600 / 1e6 * 12) %>%
  mutate(hourly_C_assim_tla = hourly_C_assim * (total_leaf_area_cm2 / 10000),
         hourly_C_resp_tla = hourly_C_resp * (total_leaf_area_cm2 / 10000),
         hourly_netC_assim_tla = hourly_netC_assim * (total_leaf_area_cm2 / 10000))

#####################################################################
# Calculate total daily C estimates
#####################################################################
daily_model_results <- hourly_model_results %>%
  group_by(doy, id, spp, plot, subplot, gm.trt) %>%
  summarize(daily_C_assim = sum(hourly_C_assim, na.rm = TRUE),
            daily_C_resp = sum(hourly_C_resp, na.rm = TRUE),
            daily_netC_assim = sum(hourly_netC_assim, na.rm = TRUE),
            
            daily_C_assim_tla = sum(hourly_C_assim_tla, na.rm = TRUE),
            daily_C_resp_tla = sum(hourly_C_resp_tla, na.rm = TRUE),
            daily_netC_assim_tla = sum(hourly_netC_assim_tla, na.rm = TRUE)) %>%
  ungroup(doy) %>%
  mutate(cumul_C_assim = cumsum(daily_C_assim) / max(daily_C_assim),
         cumul_C_resp = cumsum(daily_C_resp) / max(daily_C_resp),
         cumul_netC_assim = cumsum(daily_netC_assim)/ max(daily_netC_assim),
         
         cumul_C_assim_tla = cumsum(daily_C_assim_tla),
         cumul_C_resp_tla = cumsum(daily_C_resp_tla),
         cumul_netC_assim_tla = cumsum(daily_netC_assim_tla)) %>%
  left_join(n_measurements)
head(daily_model_results)

# Figure out maximum C assim for each id
daily_model_results %>%
  group_by(id, spp, plot, subplot, gm.trt) %>%
  summarize(max_cumul_C_assim = max(cumul_C_assim),
            max_cumul_C_resp = max(cumul_C_resp),
            max_cumul_netC_assim = max(cumul_netC_assim),
            
            max_cumul_C_assim_tla = max(cumul_C_assim_tla),
            max_cumul_C_resp_tla = max(cumul_C_resp_tla),
            max_cumul_netC_assim_tla = max(cumul_netC_assim_tla),



# Check that cumsum() worked
filter(daily_model_results, id == 9412)

# Write daily model results to .csv
write.csv(daily_model_results, "../data/c_budget/TT24_cBudget_data.csv", row.names = F)
