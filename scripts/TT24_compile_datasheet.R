# TT24_compile_datasheet.R
# R script that merges A/Ci curve parameters with snapshot
# photosynthesis data and multispeq data. Script also uses weather 
# station data to correct Vcmax and Jmax to temperature-standardized
# estimates (i.e., to a standard leaf temp of 25degC)
# 
# Note: all paths assume that the folder containing this R script
# is the working root directory

#####################################################################
# Libraries and custom functions
#####################################################################
# Libraries
library(tidyverse)
library(plantecophys)
library(lubridate)
library(zoo) ## for running mean fxn

# Load custom functions for cleaning LI-6800 files,
# standardizing Vcmax/Jmax/Rd to single temperature,
# and calculating 
R.utils::sourceDirectory("../functions/")

# Create data frame containing subplots and their treatments for
# easy merge with photosynthetic trait data
gm.ambient <- c(4, 5, 6, 10, 11, 12, 16, 17, 18,
                22, 23, 24, 28, 29, 30, 34, 35, 36)

treatments <- data.frame(subplot = seq(1,36, 1)) %>%
  mutate(gm.trt = ifelse(subplot %in% gm.ambient == TRUE, 
                         "ambient", "weeded"))

# Read .csv for A/Ci curve coefficients and snapshot photosynthesis traits
aci_coefs <- read.csv("../data/TT24_curve_fits.csv")
snapshot <- read.csv("../data/TT24_photo_snapshot.csv")

# Read .csv for weather station data
weather <- read.csv("../data/TT24_weather_station_data.csv")

#####################################################################
# Get 10-day weather average leading up to each sampling day to 
# calculate Vcmax25/Jmax25
#####################################################################
weather_dailymean <- weather %>%
  separate(date, into = c("date_only", "time"), sep = " ", remove = F) %>%
  group_by(date_only, doy) %>%
  summarize(precip_total = sum(precipitation_mm, na.rm = T),
            airtemp_mean = mean(air_temperature_c),
            vaporpressure_mean = mean(vapor_pressure_kpa),
            atmpressure_mean = mean(atm_pressure_kpa),
            vpd_mean = mean(vpd_kpa))

# write.csv(weather_dailymean, "../data/TT24_daily_weather_summary.csv", row.names = F)

# Add 10-day rolling mean (rolling total for precip) for climate data 
weather_dailymean$tavg10 <- rollmean(x = weather_dailymean$airtemp_mean, k = 10,
                                     fill = NA, align = "right")
weather_dailymean$vp10 <- rollmean(x = weather_dailymean$vaporpressure_mean, k = 10,
                                   fill = NA, align = "right")
weather_dailymean$atmpres10 <- rollmean(x = weather_dailymean$atmpressure_mean, k = 10,
                                   fill = NA, align = "right")
weather_dailymean$vpd10 <- rollmean(x = weather_dailymean$vpd_mean, k = 10,
                                   fill = NA, align = "right")

# create data frame with only 10-day rolling means
rolling_average_10 <- weather_dailymean %>%
  filter(!is.na(tavg10)) %>%
  select(date_only, doy, tavg10:vpd10)

# visualize rolling mean values
hist(rolling_average_10$tavg10)
hist(rolling_average_10$vp10)
hist(rolling_average_10$atmpres10)
hist(rolling_average_10$vpd10)


#####################################################################
# Compile A/Ci parameter estimates and snapshot measurements
#####################################################################
total_photo <- aci_coefs %>%
  full_join(snapshot, by = c("id", "doy")) %>%
  left_join(treatments, by = c("subplot")) %>%
  left_join(rolling_average_10, by = "doy") %>%
  mutate(vcmax25 = temp_standardize(estimate = Vcmax,
                                    estimate.type = "Vcmax",
                                    standard.to = 25,
                                    tLeaf = Tleaf,
                                    tGrow = tavg10),
         jmax25 = temp_standardize(estimate = Jmax,
                                   estimate.type = "Jmax",
                                   standard.to = 25,
                                   tLeaf = Tleaf,
                                   tGrow = tavg10),
         rd25 = temp_standardize(estimate = Rd,
                                 estimate.type = "Rd",
                                 pft = "C3H",
                                 standard.to = 25,
                                 tLeaf = Tleaf,
                                 tGrow = tavg10)) %>%
  dplyr::select(id, machine, date = date_only, doy, spp:subplot, gm.trt, 
                anet, ci.ca, gsw, iwue, vcmax = Vcmax, vcmax25, jmax = Jmax, 
                jmax25, rd = Rd, rd25, TPU, leaf_length_cm, stem_length_cm, 
                Tleaf, tavg10:vpd10) %>%
  mutate(across(anet:vpd10, \(x) round(x, digits = 4))) %>%
  arrange(doy, plot, subplot)

# write.csv(total_photo, "../data/TT24_photo_traits.csv", 
#           row.names = F)

