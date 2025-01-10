# TT24_compile_datasheet.R
# R script that 

merges A/Ci curve parameters with snapshot
# photosynthesis data and multispeq data. Script also uses weather 
# station data to correct Vcmax and Jmax to temperature-standardized
# estimates (i.e., to a standard leaf temp of 25degC)
# 
# Note: all paths assume that the folder containing this R script
# is the working root directory

#####################################################################
# Libraries and data files
#####################################################################
# Libraries
library(tidyverse)
library(bigleaf)
library(zoo) ## for running mean fxn

# Read .csv for weather station data
weather <- read.csv("../data/TT24_weather_station_data.csv") %>%
  mutate(par_umol_m2_s = Rg.to.PPFD(solar_radiation_wm2)) %>%
  dplyr::select(date, doy, solar_radiation_wm2, par_umol_m2_s, 
                precipitation_mm, 
                wind_direction_degrees:atm_pressure_kpa, 
                rh_sensor_temp, vpd_kpa)

# write.csv(weather, "../data/TT24_hourly_weather_data.csv",
#           row.names = F)

#####################################################################
# Get 10-day weather average leading up to each sampling day to 
# calculate Vcmax25/Jmax25
#####################################################################
weather_dailymean <- weather %>%
  separate(date, into = c("date_only", "time"), 
           sep = " ", remove = F) %>%
  group_by(date_only, doy) %>%
  summarize(precip_total = sum(precipitation_mm, na.rm = TRUE),
            par_mean = mean(par_umol_m2_s, na.rm = TRUE),
            par_max = max(par_umol_m2_s, na.rm = TRUE),
            airtemp_mean = mean(air_temperature_c),
            airtemp_max = max(air_temperature_c, na.rm = TRUE),
            vaporpressure_mean = mean(vapor_pressure_kpa, na.rm = TRUE),
            atmpressure_mean = mean(atm_pressure_kpa, na.rm = TRUE),
            vpd_mean = mean(vpd_kpa, na.rm = TRUE))

# Add 10-day rolling mean (rolling total for precip) for climate data 
weather_dailymean$parmax10 <- rollmean(x = weather_dailymean$par_max, k = 10,
                                     fill = NA, align = "right")
weather_dailymean$paravg10 <- rollmean(x = weather_dailymean$par_mean, k = 10,
                                       fill = NA, align = "right")
weather_dailymean$tavg10 <- rollmean(x = weather_dailymean$airtemp_mean, k = 10,
                                     fill = NA, align = "right")
weather_dailymean$tmax10 <- rollmean(x = weather_dailymean$airtemp_max, k = 10,
                                     fill = NA, align = "right")
weather_dailymean$vp10 <- rollmean(x = weather_dailymean$vaporpressure_mean, k = 10,
                                   fill = NA, align = "right")
weather_dailymean$atmpres10 <- rollmean(x = weather_dailymean$atmpressure_mean, k = 10,
                                        fill = NA, align = "right")
weather_dailymean$vpd10 <- rollmean(x = weather_dailymean$vpd_mean, k = 10,
                                    fill = NA, align = "right")

# Write .csv with daily weather means
write.csv(weather_dailymean, 
          "../data/TT24_daily_weather_summary.csv", 
          row.names = F)

# visualize rolling mean values
hist(weather_dailymean$parmax10)
hist(weather_dailymean$paravg10)
hist(weather_dailymean$tavg10)
hist(weather_dailymean$tmax10)
hist(weather_dailymean$vp10)
hist(weather_dailymean$atmpres10)
hist(weather_dailymean$vpd10)

# Visualize a few patterns
ggplot(data = subset(weather_dailymean, doy > 100 & doy < 250), 
       aes(x = doy, y = par10)) +
  geom_line(linewidth = 2) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 100)) +
  labs(x = "Day of Year", 
       y = expression("PAR"["mean"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")) +
  theme_bw(base_size = 18)
