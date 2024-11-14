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
climate_10day <- data.frame(doy = unique(aci_coefs$doy))

select_10day_climate <- function(doy) {
  ten_days <- doy - 10
  
  df <- df[]
}


#####################################################################
# Compile A/Ci parameter estimates and snapshot measurements
#####################################################################
total_photo <- aci_coefs %>%
  full_join(snapshot, by = c("id", "doy")) %>%
  full_join(treatments, by = c("subplot"))



#####################################################################
# Compile A/Ci parameter estimates and snapshot measurements
#####################################################################
photo_cleaned_full <- aci_coefs %>%
  mutate(doy = as.numeric(doy),
         subplot = as.numeric(subplot),
         across(Vcmax:TPU, as.numeric)) %>%
  left_join(treatments, by = c("subplot")) %>%
  full_join(snapshot, by = c("id", "doy")) %>%
  mutate(vcmax25 = temp_standardize(estimate = Vcmax,
                                    estimate.type = "Vcmax",
                                    standard.to = 25,
                                    tLeaf = Tleaf,
                                    tGrow = ),
         jmax25 = temp_standardize(estimate = Jmax,
                                   estimate.type = "Jmax",
                                   standard.to = 25,
                                   tLeaf = Tleaf,
                                   tGrow = ),
         rd25 = temp_standardize(estimate = Rd,
                                 estimate.type = "Rd",
                                 pft = "C3H",
                                 standard.to = 25,
                                 tLeaf = Tleaf,
                                 tGrow = )) %>%
  dplyr::select(id, machine, doy, spp:subplot, gm.trt, Tleaf,
                anet, ci.ca, gsw, iwue, vcmax = Vcmax, vcmax25,
                jmax = Jmax, jmax25, rd = Rd, rd25) %>%
  mutate(across(Tleaf:rd25, \(x) round(x, digits = 4))) %>%
  arrange(plot, doy)

write.csv(photo_cleaned_full,
          "../data/TT24_photo_traits_working.csv", 
          row.names = F)
