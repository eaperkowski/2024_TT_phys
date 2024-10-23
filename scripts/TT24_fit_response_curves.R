# TT24_fit_response_curves.R
# R script that merges LI-6800 files from "../cleaned_li6800/" and
# uses merged file to fit CO2 response curves
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
library(purrr)

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

#####################################################################
# Clean licor files and put in `cleaned_li6800` subfolder
#####################################################################
# clean_licor_files(directory_path = "../data/raw_li6800/",
#                   write_directory = "../data/cleaned_li6800/")

#####################################################################
# Merge cleaned LI-6800 files into single file
#####################################################################

# List and set file names within cleaned_li6800 subfolder
files <- list.files(path = "../data/cleaned_li6800",
                    recursive = T,
                    pattern = "\\.csv$",
                    full.names = T)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))

# Read all files and merge into central data frame
li6800_merged <- plyr::rbind.fill(lapply(files, read.csv)) %>%
  mutate(date = lubridate::ymd_hms(date),
         date_only = stringr::word(date, 1)) %>%
  mutate(date_only = ifelse(date_only == "2021-07-22" | date_only == "2023-06-14",
                            "2024-05-01", date_only),
         date = str_c(date_only, " ", hhmmss)) %>%
  filter(!is.na(date)) %>%
  mutate(doy = yday(date_only),
         Qin_cuvette = 2000,

         ## Some corrections
         machine = ifelse(machine == "5024", "albert", machine),

         id = ifelse(doy == 114 & id == 1 & 
                       elapsed > 6840, 2, id),
         id = ifelse(id == "striped", "striped1", id),
         id = ifelse(doy == 119 & id == 5479 &
                       elapsed > 4200 & elapsed < 5400, 4942, id),
         id = ifelse(doy == 119 & id == 5479 &
                       elapsed > 5500, 1666, id),
         id = ifelse(doy == 119 & id == 392 &
                       elapsed > 10500, 2337, id),
         id = ifelse(doy == 119 & id == 4934 &
                       elapsed > 14000, 1157, id),
         id = ifelse(doy == 120 & id == 5229 &
                       elapsed > 8500, 4990, id),
         id = ifelse(doy == 120 & id == 4990 &
                       elapsed > 11000, NA, id),
         id = ifelse(doy == 127 & id == "TT24_101" &
                       elapsed > 5700, "TT24_102", id),
         id = ifelse(doy == 128 & id == 5567 &
                       elapsed > 9500, 5569, id)) %>%
  dplyr::select(obs, time, elapsed, date, date_only, doy, hhmmss:Qin,
                Qin_cuvette, Qabs:SS_r) %>%
  arrange(machine, date, obs)

# Write merged LI-6800 file
# write.csv(li6800_merged, "../data/TT24_li6800_merged.csv", 
#           row.names = F)

#####################################################################
#####################################################################
# Start curve fits
#####################################################################
#####################################################################

# Object naming scheme: Species, tag ID, day of year separated by "_"

# Read .csv file that has all curve fits merged to it
li6800_merged <- read.csv("../data/TT24_li6800_merged.csv")

#####################################################################
# 4/13/24: plot 3
#####################################################################
li6800_merged <- read.csv("../data/TT24_li6800_merged.csv")


## 583
tri_583_104 <- subset(li6800_merged, id == "583" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_583_104)

aci_coefs <- data.frame(id = 583, spp = "Tri", plot = 3, subplot = 23, 
                        doy = 104, t(coef(tri_583_104)),
                        leaf_length_cm = 7.4, stem_length_cm = NA)

## 4934
tri_4934_104 <- subset(li6800_merged, id == "4934" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4934_104)

aci_coefs[2,] <- c(id = 4934, spp = "Tri", plot = 3, subplot = 23,
                   doy = 104, t(coef(tri_4934_104)),
                   leaf_length_cm = 6.8, stem_length_cm = NA)

## 452
tri_452_104 <- subset(li6800_merged, id == "452" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_452_104)

aci_coefs[3,] <- c(id = 452, spp = "Tri", plot = 3, subplot = 22,
                   doy = 104, t(coef(tri_452_104)),
                   leaf_length_cm = 6.9, stem_length_cm = NA)

## 6558
tri_6558_104 <- subset(li6800_merged, id == "6558" & doy == 104 & (A < 7.5 | A > 12)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6558_104)

aci_coefs[4,] <- c(id = 6558, spp = "Tri", plot = 3, subplot = 22,
                   doy = 104, t(coef(tri_6558_104)),
                   leaf_length_cm = 7.9, stem_length_cm = NA)

## 1666
tri_1666_104 <- subset(li6800_merged, id == "1666" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1666_104)

aci_coefs[5,] <- c(id = 1666, spp = "Tri", plot = 3, subplot = 27,
                   doy = 104, t(coef(tri_1666_104)),
                   leaf_length_cm = 6.8, stem_length_cm = NA)

## 4942
tri_4942_104 <- subset(li6800_merged, id == "4942" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4942_104)

aci_coefs[6,] <- c(id = 4942, spp = "Tri", plot = 3, subplot = 27,
                   doy = 104, t(coef(tri_4942_104)),
                   leaf_length_cm = 6.4, stem_length_cm = NA)

## 5479
tri_5479_104 <- subset(li6800_merged, id == "5479" & doy == 104 & (A < 7 | A > 11)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5479_104)

aci_coefs[7,] <- c(id = 5479, spp = "Tri", plot = 3, subplot = 21,
                   doy = 104, t(coef(tri_5479_104)),
                   leaf_length_cm = 6.9, stem_length_cm = NA)

## 774
tri_774_104 <- subset(li6800_merged, id == "774" & doy == 104 & (A < 11.5 | A > 14)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_774_104)

aci_coefs[8,] <- c(id = 774, spp = "Tri", plot = 3, subplot = 20,
                   doy = 104, t(coef(tri_774_104)),
                   leaf_length_cm = 9.5, stem_length_cm = NA)

## 6885
tri_6885_104 <- subset(li6800_merged, id == "6885" & doy == 104 & (A < 8 | A > 15)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6885_104)

aci_coefs[9,] <- c(id = 6885, spp = "Tri", plot = 3, subplot = 20,
                   doy = 104, t(coef(tri_6885_104)),
                   leaf_length_cm = 6.1, stem_length_cm = NA)

## 2329
tri_2329_104 <- subset(li6800_merged, id == "2329" & doy == 104 & (A < 8 | A >14)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2329_104)

aci_coefs[10,] <- c(id = 2329, spp = "Tri", plot = 3, subplot = 20,
                   doy = 104, t(coef(tri_2329_104)),
                   leaf_length_cm = 7.5, stem_length_cm = NA)

## 4714
tri_4714_104 <- subset(li6800_merged, id == "4714" & doy == 104 & (A < 12 | A > 19)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4714_104)

aci_coefs[11,] <- c(id = 4714, spp = "Tri", plot = 3, subplot = 7,
                    doy = 104, t(coef(tri_4714_104)),
                    leaf_length_cm = 8.0, stem_length_cm = NA)

## 902
tri_902_104 <- subset(li6800_merged, id == "902" & doy == 104 & (A < 11 | A > 15)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_902_104)

aci_coefs[12,] <- c(id = 902, spp = "Tri", plot = 3, subplot = 7,
                    doy = 104, t(coef(tri_902_104)),
                    leaf_length_cm = 7.4, stem_length_cm = NA)

## 552
tri_552_104 <- subset(li6800_merged, id == "552" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_552_104)

aci_coefs[13,] <- c(id = 552, spp = "Tri", plot = 3, subplot = 14,
                    doy = 104, t(coef(tri_552_104)),
                    leaf_length_cm = 8.3, stem_length_cm = NA)

## 5436
tri_5436_104 <- subset(li6800_merged, id == "5436" & doy == 104 & (A < 7 | A > 14)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5436_104)

aci_coefs[14,] <- c(id = 5436, spp = "Tri", plot = 3, subplot = 3,
                    doy = 104, t(coef(tri_5436_104)),
                    leaf_length_cm = 9.2, stem_length_cm = NA)

## 3563
tri_3563_104 <- subset(li6800_merged, id == "3563" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3563_104)

aci_coefs[15,] <- c(id = 3563, spp = "Tri", plot = 3, subplot = 4,
                    doy = 104, t(coef(tri_3563_104)),
                    leaf_length_cm = 6.8, stem_length_cm = NA)

## 1926
tri_1926_104 <- subset(li6800_merged, id == "1926" & doy == 104 & (A < 5 | A > 11)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1926_104)

aci_coefs[16,] <- c(id = 1926, spp = "Tri", plot = 3, subplot = 5,
                    doy = 104, t(coef(tri_1926_104)),
                    leaf_length_cm = 6.7, stem_length_cm = NA)

## 425
tri_425_104 <- subset(li6800_merged, id == "425" & doy == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_425_104)

aci_coefs[17,] <- c(id = 425, spp = "Tri", plot = 3, subplot = 6,
                    doy = 104, t(coef(tri_425_104)),
                    leaf_length_cm = 8.2, stem_length_cm = NA)

#####################################################################
# 4/14/24: plot 5
#####################################################################

# 2924
tri_2924_105 <- subset(li6800_merged, id == "2924" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2924_105)

aci_coefs[18,] <- c(id = 2924, spp = "Tri", plot = 5, subplot = 35,
                    doy = 105, t(coef(tri_2924_105)),
                    leaf_length_cm = 8.3, stem_length_cm = NA)

# 6875
tri_6875_105 <- subset(li6800_merged, id == "6875" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6875_105)

aci_coefs[19,] <- c(id = 6875, spp = "Tri", plot = 5, subplot = 29,
                    doy = 105, t(coef(tri_6875_105)),
                    leaf_length_cm = 8.0, stem_length_cm = NA)

# 2988
tri_2988_105 <- subset(li6800_merged, id == "2988" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2988_105)

aci_coefs[20,] <- c(id = 2988, spp = "Tri", plot = 5, subplot = 33,
                    doy = 105, t(coef(tri_2988_105)),
                    leaf_length_cm = 9.5, stem_length_cm = NA)

# 3829
tri_3829_105 <- subset(li6800_merged, id == "3829" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 300)
# plot(tri_3829_105)

aci_coefs[21,] <- c(id = 3829, spp = "Tri", plot = 5, subplot = 33,
                    doy = 105, t(coef(tri_3829_105)),
                    leaf_length_cm = 7.2, stem_length_cm = NA)

# 5877
tri_5877_105 <- subset(li6800_merged, id == "5877" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 800)
# plot(tri_5877_105)

aci_coefs[22,] <- c(id = 5877, spp = "Tri", plot = 5, subplot = 27,
                    doy = 105, t(coef(tri_5877_105)),
                    leaf_length_cm = 7.5, stem_length_cm = NA)

# 4109
tri_4109_105 <- subset(li6800_merged, id == "4109" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4109_105)

aci_coefs[23,] <- c(id = 4109, spp = "Tri", plot = 5, subplot = 32,
                    doy = 105, t(coef(tri_4109_105)),
                    leaf_length_cm = 7.0, stem_length_cm = NA)

# 4431
tri_4431_105 <- subset(li6800_merged, id == "4431" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4431_105)

aci_coefs[24,] <- c(id = 4431, spp = "Tri", plot = 5, subplot = 23,
                    doy = 105, t(coef(tri_4431_105)),
                    leaf_length_cm = 10.2, stem_length_cm = NA)

# 4959
tri_4959_105 <- subset(li6800_merged, id == "4959" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4959_105)

aci_coefs[25,] <- c(id = 4959, spp = "Tri", plot = 5, subplot = 22,
                    doy = 105, t(coef(tri_4959_105)),
                    leaf_length_cm = 9.2, stem_length_cm = NA)

# 5229
tri_5229_105 <- subset(li6800_merged, id == "5229" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5229_105)

aci_coefs[26,] <- c(id = 5229, spp = "Tri", plot = 5, subplot = 22,
                    doy = 105, t(coef(tri_5229_105)),
                    leaf_length_cm = 8.4, stem_length_cm = NA)

# 482
tri_482_105 <- subset(li6800_merged, id == "482" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_482_105)

aci_coefs[27,] <- c(id = 482, spp = "Tri", plot = 5, subplot = 22,
                    doy = 105, t(coef(tri_482_105)),
                    leaf_length_cm = 9.5, stem_length_cm = NA)

# 4990
tri_4990_105 <- subset(li6800_merged, id == "4990" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4990_105)

aci_coefs[28,] <- c(id = 4990, spp = "Tri", plot = 5, subplot = 22,
                    doy = 105, t(coef(tri_4990_105)),
                    leaf_length_cm = 11.2, stem_length_cm = NA)

# 3004
tri_3004_105 <- subset(li6800_merged, id == "3004" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3004_105)

aci_coefs[29,] <- c(id = 3004, spp = "Tri", plot = 5, subplot = 21,
                    doy = 105, t(coef(tri_3004_105)),
                    leaf_length_cm = 9.6, stem_length_cm = NA)

# 4576
tri_4576_105 <- subset(li6800_merged, id == "4576" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4576_105)

aci_coefs[30,] <- c(id = 4576, spp = "Tri", plot = 5, subplot = 21,
                    doy = 105, t(coef(tri_4576_105)),
                    leaf_length_cm = 9.8, stem_length_cm = NA)

# 3077
tri_3077_105 <- subset(li6800_merged, id == "3077" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 650)
# plot(tri_3077_105)

aci_coefs[31,] <- c(id = 3077, spp = "Tri", plot = 5, subplot = 21,
                    doy = 105, t(coef(tri_3077_105)),
                    leaf_length_cm = 8.5, stem_length_cm = NA)

# 4000
tri_4000_105 <- subset(li6800_merged, id == "4000" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 400)
# plot(tri_4000_105)

aci_coefs[32,] <- c(id = 4000, spp = "Tri", plot = 5, subplot = 20,
                    doy = 105, t(coef(tri_4000_105)),
                    leaf_length_cm = 8.1, stem_length_cm = NA)

# 5115
tri_5115_105 <- subset(li6800_merged, id == "5115" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5115_105)

aci_coefs[33,] <- c(id = 5115, spp = "Tri", plot = 5, subplot = 19,
                    doy = 105, t(coef(tri_5115_105)),
                    leaf_length_cm = 9.9, stem_length_cm = NA)

# 5228
tri_5228_105 <- subset(li6800_merged, id == "5228" & doy == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5228_105)

aci_coefs[34,] <- c(id = 5228, spp = "Tri", plot = 5, subplot = 19,
                    doy = 105, t(coef(tri_5228_105)),
                    leaf_length_cm = 8.0, stem_length_cm = NA)

#####################################################################
# 4/15/24: plot 6
#####################################################################
# 5229
tri_5229_106 <- subset(li6800_merged, id == "5229" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 600)
# plot(tri_5229_106)

aci_coefs[35,] <- c(id = 5229, spp = "Tri", plot = 6, subplot = 1,
                    doy = 106, t(coef(tri_5229_106)),
                    leaf_length_cm = 10.1, stem_length_cm = NA)

# 762
tri_762_106 <- subset(li6800_merged, id == "762" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 600)
# plot(tri_762_106)

aci_coefs[36,] <- c(id = 762, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_762_106)),
                    leaf_length_cm = 7.6, stem_length_cm = NA)

# 4373
tri_4373_106 <- subset(li6800_merged, id == "4373" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4373_106)

aci_coefs[37,] <- c(id = 4373, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_4373_106)),
                    leaf_length_cm = 8.5, stem_length_cm = NA)

# 4760
tri_4760_106 <- subset(li6800_merged, id == "4760" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4760_106)

aci_coefs[38,] <- c(id = 4760, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_4760_106)),
                    leaf_length_cm = 6.8, stem_length_cm = NA)

# 4777
tri_4777_106 <- subset(li6800_merged, id == "4777" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4777_106)

aci_coefs[39,] <- c(id = 4777, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_4777_106)),
                    leaf_length_cm = 6.9, stem_length_cm = NA)

# 5865
tri_5865_106 <- subset(li6800_merged, id == "5865" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5865_106)

aci_coefs[40,] <- c(id = 5865, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_5865_106)),
                    leaf_length_cm = 10.0, stem_length_cm = NA)

# 5267
tri_5267_106 <- subset(li6800_merged, id == "5267" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5267_106)

aci_coefs[41,] <- c(id = 5267, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_5267_106)),
                    leaf_length_cm = 8.1, stem_length_cm = NA)

# 684
tri_684_106 <- subset(li6800_merged, id == "684" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_684_106)

aci_coefs[42,] <- c(id = 684, spp = "Tri", plot = 6, subplot = 2,
                    doy = 106, t(coef(tri_684_106)),
                    leaf_length_cm = 10.5, stem_length_cm = NA)

# 5742
tri_5742_106 <- subset(li6800_merged, id == "5742" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5742_106)

aci_coefs[43,] <- c(id = 5742, spp = "Tri", plot = 6, subplot = 3,
                    doy = 106, t(coef(tri_5742_106)),
                    leaf_length_cm = 9.5, stem_length_cm = NA)

# 3305
tri_3305_106 <- subset(li6800_merged, id == "3305" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3305_106)

aci_coefs[44,] <- c(id = 3305, spp = "Tri", plot = 6, subplot = 5,
                    doy = 106, t(coef(tri_3305_106)),
                    leaf_length_cm = 8.3, stem_length_cm = NA)

# 5619
tri_5619_106 <- subset(li6800_merged, id == "5619" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5619_106)

aci_coefs[45,] <- c(id = 5619, spp = "Tri", plot = 6, subplot = 5,
                    doy = 106, t(coef(tri_5619_106)),
                    leaf_length_cm = 9.3, stem_length_cm = NA)

# striped1
tri_striped1_106 <- subset(li6800_merged, id == "striped1" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_striped1_106)

aci_coefs[46,] <- c(id = "striped1", spp = "Tri", plot = 6, subplot = 5,
                    doy = 106, t(coef(tri_striped1_106)),
                    leaf_length_cm = 8.1, stem_length_cm = NA)

# 978
tri_978_106 <- subset(li6800_merged, id == "978" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_978_106)

aci_coefs[47,] <- c(id = 978, spp = "Tri", plot = 6, subplot = 6,
                    doy = 106, t(coef(tri_978_106)),
                    leaf_length_cm = 8.3, stem_length_cm = NA)

# 5626
tri_5626_106 <- subset(li6800_merged, id == "5626" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5626_106)

aci_coefs[48,] <- c(id = 5626, spp = "Tri", plot = 6, subplot = 6,
                    doy = 106, t(coef(tri_5626_106)),
                    leaf_length_cm = 10.0, stem_length_cm = NA)

# 5403
tri_5403_106 <- subset(li6800_merged, id == "5403" & doy == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5403_106)

aci_coefs[49,] <- c(id = 5403, spp = "Tri", plot = 6, subplot = 6,
                    doy = 106, t(coef(tri_5403_106)),
                    leaf_length_cm = 9.0, stem_length_cm = NA)

#####################################################################
# 4/20/24: plot 6
#####################################################################

# 5060
tri_5060_111 <- subset(li6800_merged, id == "5060" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5060_111)

aci_coefs[50,] <- c(id = 5060, spp = "Tri", plot = 6, subplot = 27,
                    doy = 111, t(coef(tri_5060_111)),
                    leaf_length_cm = 10.3, stem_length_cm = NA)

# 5783
tri_5783_111 <- subset(li6800_merged, id == "5783" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5783_111)

aci_coefs[51,] <- c(id = 5783, spp = "Tri", plot = 6, subplot = 20,
                    doy = 111, t(coef(tri_5783_111)),
                    leaf_length_cm = 10.6, stem_length_cm = NA)

# 5031
tri_5031_111 <- subset(li6800_merged, id == "5031" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5031_111)

aci_coefs[52,] <- c(id = 5031, spp = "Tri", plot = 6, subplot = 28,
                    doy = 111, t(coef(tri_5031_111)),
                    leaf_length_cm = 9.4, stem_length_cm = NA)

# 4543
tri_4543_111 <- subset(li6800_merged, id == "4543" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4543_111)

aci_coefs[53,] <- c(id = 4543, spp = "Tri", plot = 6, subplot = 16,
                    doy = 111, t(coef(tri_4543_111)),
                    leaf_length_cm = 8.8, stem_length_cm = NA)

# 2276
tri_2276_111 <- subset(li6800_merged, id == "2276" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2276_111)

aci_coefs[54,] <- c(id = 2276, spp = "Tri", plot = 6, subplot = 22,
                    doy = 111, t(coef(tri_2276_111)),
                    leaf_length_cm = 6.8, stem_length_cm = NA)

# 4511
tri_4511_111 <- subset(li6800_merged, id == "4511" & doy == 111 & (A < 3 | A > 14)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4511_111)

aci_coefs[55,] <- c(id = 4511, spp = "Tri", plot = 6, subplot = 16,
                    doy = 111, t(coef(tri_4511_111)),
                    leaf_length_cm = 11.6, stem_length_cm = NA)

# 3770
tri_3770_111 <- subset(li6800_merged, id == "3770" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3770_111)

aci_coefs[56,] <- c(id = 3770, spp = "Tri", plot = 6, subplot = 10,
                    doy = 111, t(coef(tri_3770_111)),
                    leaf_length_cm = 8.1, stem_length_cm = NA)

# 4582
tri_4582_111 <- subset(li6800_merged, id == "4582" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4582_111)

aci_coefs[57,] <- c(id = 4582, spp = "Tri", plot = 6, subplot = 11,
                    doy = 111, t(coef(tri_4582_111)),
                    leaf_length_cm = 8.8, stem_length_cm = NA)

# 2980
tri_2980_111 <- subset(li6800_merged, id == "2980" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2980_111)

aci_coefs[58,] <- c(id = 2980, spp = "Tri", plot = 6, subplot = 11,
                    doy = 111, t(coef(tri_2980_111)),
                    leaf_length_cm = 9.7, stem_length_cm = NA)

# 3371
tri_3371_111 <- subset(li6800_merged, id == "3371" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3371_111)

aci_coefs[59,] <- c(id = 3371, spp = "Tri", plot = 6, subplot = 11,
                    doy = 111, t(coef(tri_3371_111)),
                    leaf_length_cm = 11.2, stem_length_cm = NA)

# 3379
tri_3379_111 <- subset(li6800_merged, id == "3379" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3379_111)

aci_coefs[60,] <- c(id = 3379, spp = "Tri", plot = 6, subplot = 11,
                    doy = 111, t(coef(tri_3379_111)),
                    leaf_length_cm = 6.2, stem_length_cm = NA)

# 4505
tri_4505_111 <- subset(li6800_merged, id == "4505" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4505_111)

aci_coefs[61,] <- c(id = 4505, spp = "Tri", plot = 6, subplot = 11,
                    doy = 111, t(coef(tri_4505_111)),
                    leaf_length_cm = 12.0, stem_length_cm = NA)

# 4745
tri_4745_111 <- subset(li6800_merged, id == "4745" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4745_111)

aci_coefs[62,] <- c(id = 4745, spp = "Tri", plot = 6, subplot = 12,
                    doy = 111, t(coef(tri_4745_111)),
                    leaf_length_cm = 13.1, stem_length_cm = NA)

# 5722
tri_5722_111 <- subset(li6800_merged, id == "5722" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5722_111)

aci_coefs[63,] <- c(id = 5722, spp = "Tri", plot = 6, subplot = 6,
                    doy = 111, t(coef(tri_5722_111)),
                    leaf_length_cm = 13.2, stem_length_cm = NA)

# 2921
tri_2921_111 <- subset(li6800_merged, id == "2921" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2921_111)

aci_coefs[64,] <- c(id = 2921, spp = "Tri", plot = 6, subplot = 6,
                    doy = 111, t(coef(tri_2921_111)),
                    leaf_length_cm = 11.1, stem_length_cm = NA)

# striped2
tri_striped2_111 <- subset(li6800_merged, id == "striped2" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_striped2_111)

aci_coefs[65,] <- c(id = "striped2", spp = "Tri", plot = 6, subplot = 5,
                    doy = 111, t(coef(tri_striped2_111)),
                    leaf_length_cm = 10.2, stem_length_cm = NA)

# 5642
tri_5642_111 <- subset(li6800_merged, id == "5642" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5642_111)

aci_coefs[66,] <- c(id = 5642, spp = "Tri", plot = 6, subplot = 4,
                    doy = 111, t(coef(tri_5642_111)),
                    leaf_length_cm = 9.9, stem_length_cm = NA)

# 6879
tri_6879_111 <- subset(li6800_merged, id == "6879" & doy == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6879_111)

aci_coefs[67,] <- c(id = 6879, spp = "Tri", plot = 6, subplot = 36,
                    doy = 111, t(coef(tri_6879_111)),
                    leaf_length_cm = 10.3, stem_length_cm = NA)

#####################################################################
# 4/22/24: plot 5
#####################################################################

# 7120
tri_7120_113 <- subset(li6800_merged, id == "7120" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_7120_113)

aci_coefs[68,] <- c(id = 7120, spp = "Tri", plot = 5, subplot = 17,
                    doy = 113, t(coef(tri_7120_113)),
                    leaf_length_cm = 10.1, stem_length_cm = NA)

# 2912
tri_2912_113 <- subset(li6800_merged, id == "2912" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2912_113)

aci_coefs[69,] <- c(id = 2912, spp = "Tri", plot = 5, subplot = 11,
                    doy = 113, t(coef(tri_2912_113)),
                    leaf_length_cm = 13.0, stem_length_cm = NA)

# 614
tri_614_113 <- subset(li6800_merged, id == "614" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_614_113)

aci_coefs[70,] <- c(id = 614, spp = "Tri", plot = 5, subplot = 11,
                    doy = 113, t(coef(tri_614_113)),
                    leaf_length_cm = NA, stem_length_cm = NA)

# 1374
tri_1374_113 <- subset(li6800_merged, id == "1374" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1374_113)

aci_coefs[71,] <- c(id = 1374, spp = "Tri", plot = 5, subplot = 10,
                    doy = 113, t(coef(tri_1374_113)),
                    leaf_length_cm = 11.9, stem_length_cm = NA)

# 6881
tri_6881_113 <- subset(li6800_merged, id == "6881" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6881_113)

aci_coefs[72,] <- c(id = 6881, spp = "Tri", plot = 5, subplot = 10,
                    doy = 113, t(coef(tri_6881_113)),
                    leaf_length_cm = 11.8, stem_length_cm = NA)

# 5381
tri_5381_113 <- subset(li6800_merged, id == "5381" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5381_113)

aci_coefs[73,] <- c(id = 5381, spp = "Tri", plot = 5, subplot = 10,
                    doy = 113, t(coef(tri_5381_113)),
                    leaf_length_cm = 12.7, stem_length_cm = NA)

# 4105
tri_4105_113 <- subset(li6800_merged, id == "4105" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4105_113)

aci_coefs[74,] <- c(id = 4105, spp = "Tri", plot = 5, subplot = 10,
                    doy = 113, t(coef(tri_4105_113)),
                    leaf_length_cm = 11.6, stem_length_cm = NA)

# 4547
tri_4547_113 <- subset(li6800_merged, id == "4547" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4547_113)

aci_coefs[75,] <- c(id = 4547, spp = "Tri", plot = 5, subplot = 10,
                    doy = 113, t(coef(tri_4547_113)),
                    leaf_length_cm = 12.1, stem_length_cm = NA)

# 2388
tri_2388_113 <- subset(li6800_merged, id == "2388" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2388_113)

aci_coefs[76,] <- c(id = 2388, spp = "Tri", plot = 5, subplot = 10,
                    doy = 113, t(coef(tri_2388_113)),
                    leaf_length_cm = 14.3, stem_length_cm = NA)

# 5904
tri_5904_113 <- subset(li6800_merged, id == "5904" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5904_113)

aci_coefs[77,] <- c(id = 5904, spp = "Tri", plot = 5, subplot = 9,
                    doy = 113, t(coef(tri_5904_113)),
                    leaf_length_cm = 11.0, stem_length_cm = NA)

# 7147
tri_7147_113 <- subset(li6800_merged, id == "7147" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_7147_113)

aci_coefs[78,] <- c(id = 7147, spp = "Tri", plot = 5, subplot = 9,
                    doy = 113, t(coef(tri_7147_113)),
                    leaf_length_cm = 12.8, stem_length_cm = NA)

# 43
tri_43_113 <- subset(li6800_merged, id == "43" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_43_113)

aci_coefs[79,] <- c(id = 43, spp = "Tri", plot = 5, subplot = 9,
                    doy = 113, t(coef(tri_43_113)),
                    leaf_length_cm = 12.7, stem_length_cm = NA)

# 86
tri_86_113 <- subset(li6800_merged, id == "86" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_86_113)

aci_coefs[80,] <- c(id = 86, spp = "Tri", plot = 5, subplot = 9,
                    doy = 113, t(coef(tri_86_113)),
                    leaf_length_cm = 11.2, stem_length_cm = NA)

# flag8
tri_flag8_113 <- subset(li6800_merged, id == "flag8" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_flag8_113)

aci_coefs[81,] <- c(id = "flag8", spp = "Tri", plot = 5, subplot = 16,
                    doy = 113, t(coef(tri_flag8_113)),
                    leaf_length_cm = 11.9, stem_length_cm = NA)

# 2563
tri_2563_113 <- subset(li6800_merged, id == "2563" & doy == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2563_113)

aci_coefs[82,] <- c(id = 2563, spp = "Tri", plot = 5, subplot = 9,
                    doy = 113, t(coef(tri_2563_113)),
                    leaf_length_cm = NA, stem_length_cm = NA)

#####################################################################
# 4/23/24: plot 5
#####################################################################

# 4265
tri_4265_114 <- subset(li6800_merged, id == "4265" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4265_114)

aci_coefs[83,] <- c(id = 4265, spp = "Tri", plot = 5, subplot = 31,
                    doy = 114, t(coef(tri_4265_114)),
                    leaf_length_cm = 14.8, stem_length_cm = NA)

# 2573
tri_2573_114 <- subset(li6800_merged, id == "2573" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2573_114)

aci_coefs[84,] <- c(id = 2573, spp = "Tri", plot = 5, subplot = 31,
                    doy = 114, t(coef(tri_2573_114)),
                    leaf_length_cm = 14.1, stem_length_cm = NA)

# 2547
tri_2547_114 <- subset(li6800_merged, id == "2547" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2547_114)

aci_coefs[85,] <- c(id = 2547, spp = "Tri", plot = 5, subplot = 19,
                    doy = 114, t(coef(tri_2547_114)),
                    leaf_length_cm = 14.0, stem_length_cm = NA)

# 4177
tri_4177_114 <- subset(li6800_merged, id == "4177" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4177_114)

aci_coefs[86,] <- c(id = 4177, spp = "Tri", plot = 5, subplot = 14,
                    doy = 114, t(coef(tri_4177_114)),
                    leaf_length_cm = 15.0, stem_length_cm = NA)

# 1795
tri_1795_114 <- subset(li6800_merged, id == "1795" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1795_114)

aci_coefs[87,] <- c(id = 1795, spp = "Tri", plot = 5, subplot = 14,
                    doy = 114, t(coef(tri_1795_114)),
                    leaf_length_cm = 9.8, stem_length_cm = NA)

# 3
tri_flag3_114 <- subset(li6800_merged, id == "3" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_flag3_114)

aci_coefs[88,] <- c(id = "flag3_tri", spp = "Tri", plot = 5, subplot = 14,
                    doy = 114, t(coef(tri_flag3_114)),
                    leaf_length_cm = 12.6, stem_length_cm = NA)

# 1
tri_flag1_114 <- subset(li6800_merged, id == "1" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_flag1_114)

aci_coefs[89,] <- c(id = "flag1_tri", spp = "Tri", plot = 5, subplot = 14,
                    doy = 114, t(coef(tri_flag1_114)),
                    leaf_length_cm = 12.0, stem_length_cm = NA)

# 2
tri_flag2_114 <- subset(li6800_merged, id == "2" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_flag2_114)

aci_coefs[90,] <- c(id = "flag2_tri", spp = "Tri", plot = 5, subplot = 14,
                    doy = 114, t(coef(tri_flag2_114)),
                    leaf_length_cm = 13.3, stem_length_cm = NA)

# 1739
tri_1739_114 <- subset(li6800_merged, id == "1739" & doy == 114 & (A < 5 | A > 14)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1739_114)

aci_coefs[91,] <- c(id = 1739, spp = "Tri", plot = 5, subplot = 14,
                    doy = 114, t(coef(tri_1739_114)),
                    leaf_length_cm = 11.0, stem_length_cm = NA)

# 4149
tri_4149_114 <- subset(li6800_merged, id == "4149" & doy == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4149_114)

aci_coefs[92,] <- c(id = 4149, spp = "Tri", plot = 5, subplot = 13,
                    doy = 114, t(coef(tri_4149_114)),
                    leaf_length_cm = 13.8, stem_length_cm = NA)

#####################################################################
# 4/28/24: plot 3
#####################################################################

# 2329
tri_2329_119 <- subset(li6800_merged, id == "2329" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2329_119)

aci_coefs[93,] <- c(id = 2329, spp = "Tri", plot = 3, subplot = 20,
                    doy = 119, t(coef(tri_2329_119)),
                    leaf_length_cm = 13.6, stem_length_cm = NA)

# 6885
tri_6885_119 <- subset(li6800_merged, id == "6885" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6885_119)

aci_coefs[94,] <- c(id = 6885, spp = "Tri", plot = 3, subplot = 20,
                    doy = 119, t(coef(tri_6885_119)),
                    leaf_length_cm = 11.5, stem_length_cm = NA)

# 774
tri_774_119 <- subset(li6800_merged, id == "774" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_774_119)

aci_coefs[95,] <- c(id = 774, spp = "Tri", plot = 3, subplot = 20,
                    doy = 119, t(coef(tri_774_119)),
                    leaf_length_cm = 11.6, stem_length_cm = NA)

# 5479
tri_5479_119 <- subset(li6800_merged, id == "5479" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5479_119)

aci_coefs[96,] <- c(id = 5479, spp = "Tri", plot = 3, subplot = 21,
                    doy = 119, t(coef(tri_5479_119)),
                    leaf_length_cm = 9.7, stem_length_cm = NA)

# 4942
tri_4942_119 <- subset(li6800_merged, id == "4942" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4942_119)

aci_coefs[97,] <- c(id = 4942, spp = "Tri", plot = 3, subplot = 27,
                    doy = 119, t(coef(tri_4942_119)),
                    leaf_length_cm = 9.7, stem_length_cm = NA)

# 1666
tri_1666_119 <- subset(li6800_merged, id == "1666" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1666_119)

aci_coefs[98,] <- c(id = 1666, spp = "Tri", plot = 3, subplot = 27,
                    doy = 119, t(coef(tri_1666_119)),
                    leaf_length_cm = 10.4, stem_length_cm = NA)

# 6558
tri_6558_119 <- subset(li6800_merged, id == "6558" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6558_119)

aci_coefs[99,] <- c(id = 6558, spp = "Tri", plot = 3, subplot = 22,
                    doy = 119, t(coef(tri_6558_119)),
                    leaf_length_cm = 13.9, stem_length_cm = NA)

# 392
tri_392_119 <- subset(li6800_merged, id == "392" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_392_119)

aci_coefs[100,] <- c(id = 392, spp = "Tri", plot = 3, subplot = 22,
                    doy = 119, t(coef(tri_392_119)),
                    leaf_length_cm = 9.5, stem_length_cm = NA)

# 2337
mai_2337_119 <- subset(li6800_merged, id == "2337" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2337_119)

aci_coefs[101,] <- c(id = 2337, spp = "Mai", plot = 3, subplot = 22,
                     doy = 119, t(coef(mai_2337_119)),
                     leaf_length_cm = 10.8, stem_length_cm = 24.5)

# 583
tri_583_119 <- subset(li6800_merged, id == "583" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_583_119)

aci_coefs[102,] <- c(id = 583, spp = "Tri", plot = 3, subplot = 23,
                     doy = 119, t(coef(tri_583_119)),
                     leaf_length_cm = 10.7, stem_length_cm = NA)

# 4934
tri_4934_119 <- subset(li6800_merged, id == "4934" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4934_119)

aci_coefs[103,] <- c(id = 4934, spp = "Tri", plot = 3, subplot = 23,
                     doy = 119, t(coef(tri_4934_119)),
                     leaf_length_cm = 11.5, stem_length_cm = NA)

# 1157
mai_1157_119 <- subset(li6800_merged, id == "1157" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_1157_119)

aci_coefs[104,] <- c(id = 1157, spp = "Mai", plot = 3, subplot = 30,
                     doy = 119, t(coef(mai_1157_119)),
                     leaf_length_cm = 10.9, stem_length_cm = 22.9)

# 179
mai_179_119 <- subset(li6800_merged, id == "174" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_179_119)

aci_coefs[105,] <- c(id = 179, spp = "Mai", plot = 3, subplot = 30,
                     doy = 119, t(coef(mai_179_119)),
                     leaf_length_cm = 16.8, stem_length_cm = 57.9)

# 4714
tri_4714_119 <- subset(li6800_merged, id == "4714" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4714_119)

aci_coefs[106,] <- c(id = 4714, spp = "Tri", plot = 3, subplot = 7,
                     doy = 119, t(coef(tri_4714_119)),
                     leaf_length_cm = 11.6, stem_length_cm = NA)

# 902
tri_902_119 <- subset(li6800_merged, id == "902" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_902_119)

aci_coefs[107,] <- c(id = 902, spp = "Tri", plot = 3, subplot = 7,
                     doy = 119, t(coef(tri_902_119)),
                     leaf_length_cm = 9.6, stem_length_cm = NA)

# 5105
mai_5105_119 <- subset(li6800_merged, id == "5105" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5105_119)

aci_coefs[108,] <- c(id = 5105, spp = "Mai", plot = 3, subplot = 8,
                     doy = 119, t(coef(mai_5105_119)),
                     leaf_length_cm = 15.5, stem_length_cm = 38.1)

# 5184
mai_5184_119 <- subset(li6800_merged, id == "5184" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5184_119)

aci_coefs[109,] <- c(id = 5184, spp = "Mai", plot = 3, subplot = 13,
                     doy = 119, t(coef(mai_5184_119)),
                     leaf_length_cm = 16.1, stem_length_cm = 42.0)

# 141
mai_141_119 <- subset(li6800_merged, id == "141" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_141_119)

aci_coefs[110,] <- c(id = 141, spp = "Mai", plot = 3, subplot = 13,
                     doy = 119, t(coef(mai_141_119)),
                     leaf_length_cm = 16.5, stem_length_cm = 43.4)

# 552
tri_552_119 <- subset(li6800_merged, id == "552" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_552_119)

aci_coefs[111,] <- c(id = 552, spp = "Tri", plot = 3, subplot = 8,
                     doy = 119, t(coef(tri_552_119)),
                     leaf_length_cm = 12.5, stem_length_cm = NA)

# 4250
mai_4250_119 <- subset(li6800_merged, id == "4250" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_4250_119)

aci_coefs[112,] <- c(id = 4250, spp = "Mai", plot = 3, subplot = 2,
                     doy = 119, t(coef(mai_4250_119)),
                     leaf_length_cm = 13.2, stem_length_cm = 37.7)

# 5488
tri_5488_119 <- subset(li6800_merged, id == "5488" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5488_119)

aci_coefs[113,] <- c(id = 5488, spp = "Tri", plot = 3, subplot = 1,
                     doy = 119, t(coef(tri_5488_119)),
                     leaf_length_cm = 16.5, stem_length_cm = NA)

# 9412
mai_9412_119 <- subset(li6800_merged, id == "9412" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_9412_119)

aci_coefs[114,] <- c(id = 9412, spp = "Mai", plot = 3, subplot = 3,
                     doy = 119, t(coef(mai_9412_119)),
                     leaf_length_cm = 16.1, stem_length_cm = 46.4)

# 5436
tri_5436_119 <- subset(li6800_merged, id == "5436" & doy == 119 & (A < 5 | A > 12)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5436_119)

aci_coefs[115,] <- c(id = 5436, spp = "Tri", plot = 3, subplot = 3,
                     doy = 119, t(coef(tri_5436_119)),
                     leaf_length_cm = 14.6, stem_length_cm = NA)

# 5495
tri_5495_119 <- subset(li6800_merged, id == "5495" & doy == 119 & (A > 0 & (A < 1 | A > 9))) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 400)
# plot(tri_5495_119)

aci_coefs[116,] <- c(id = 5495, spp = "Tri", plot = 3, subplot = 3,
                     doy = 119, t(coef(tri_5495_119)),
                     leaf_length_cm = 8.6, stem_length_cm = NA)

# 3563
tri_3563_119 <- subset(li6800_merged, id == "3563" & doy == 119 & (A < 4 | A > 10)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3563_119)

aci_coefs[117,] <- c(id = 3563, spp = "Tri", plot = 3, subplot = 4,
                     doy = 119, t(coef(tri_3563_119)),
                     leaf_length_cm = 9.5, stem_length_cm = NA)

# 2310
mai_2310_119 <- subset(li6800_merged, id == "2310" & doy == 119 & (A < 4 | A > 11)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2310_119)

aci_coefs[118,] <- c(id = 2310, spp = "Mai", plot = 3, subplot = 4,
                     doy = 119, t(coef(mai_2310_119)),
                     leaf_length_cm = 13.6, stem_length_cm = 36.4)

# 1926
tri_1926_119 <- subset(li6800_merged, id == "1926" & doy == 119 & (A < 4 | A > 10)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1926_119)

aci_coefs[119,] <- c(id = 1926, spp = "Tri", plot = 3, subplot = 5,
                     doy = 119, t(coef(tri_1926_119)),
                     leaf_length_cm = 12.3, stem_length_cm = NA)

# 5500
tri_5500_119 <- subset(li6800_merged, id == "5500" & doy == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5500_119)

aci_coefs[120,] <- c(id = 5500, spp = "Tri", plot = 3, subplot = 5,
                     doy = 119, t(coef(tri_5500_119)),
                     leaf_length_cm = 10.9, stem_length_cm = NA)

# 1686
tri_1686_119 <- subset(li6800_merged, id == "1686" & doy == 119 & (A < 4 | A > 10)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1686_119)

aci_coefs[121,] <- c(id = 1686, spp = "Tri", plot = 3, subplot = 6,
                     doy = 119, t(coef(tri_1686_119)),
                     leaf_length_cm = 15.5, stem_length_cm = NA)

# 5797
tri_5797_119 <- subset(li6800_merged, id == "5797" & doy == 119 & (A < 4 | A > 10)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5797_119)

aci_coefs[122,] <- c(id = 5797, spp = "Tri", plot = 3, subplot = 6,
                     doy = 119, t(coef(tri_5797_119)),
                     leaf_length_cm = 13.6, stem_length_cm = NA)

# 425
tri_425_119 <- subset(li6800_merged, id == "425" & doy == 119 & (A < 4 | A > 11)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_425_119)

aci_coefs[123,] <- c(id = 425, spp = "Tri", plot = 3, subplot = 6,
                     doy = 119, t(coef(tri_425_119)),
                     leaf_length_cm = 12.4, stem_length_cm = NA)

#####################################################################
# 4/29/24: plot 5
#####################################################################

# 4431
tri_4431_120 <- subset(li6800_merged, id == "4431" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4431_120)

aci_coefs[124,] <- c(id = 4431, spp = "Tri", plot = 5, subplot = 23,
                     doy = 120, t(coef(tri_4431_120)),
                     leaf_length_cm = 17.0, stem_length_cm = NA)

# 1476
mai_1476_120 <- subset(li6800_merged, id == "1476" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_1476_120)

aci_coefs[125,] <- c(id = 1476, spp = "Mai", plot = 5, subplot = 23,
                     doy = 120, t(coef(mai_1476_120)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 5657
mai_5657_120 <- subset(li6800_merged, id == "5657" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5657_120)

aci_coefs[126,] <- c(id = 5657, spp = "Mai", plot = 5, subplot = 23,
                     doy = 120, t(coef(mai_5657_120)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 4781
mai_4781_120 <- subset(li6800_merged, id == "4781" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_4781_120)

aci_coefs[127,] <- c(id = 4781, spp = "Mai", plot = 5, subplot = 23,
                     doy = 120, t(coef(mai_4781_120)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 4959
tri_4959_120 <- subset(li6800_merged, id == "4959" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4959_120)

aci_coefs[128,] <- c(id = 4959, spp = "Tri", plot = 5, subplot = 22,
                     doy = 120, t(coef(tri_4959_120)),
                     leaf_length_cm = 15.6, stem_length_cm = NA)

# 4414
tri_4414_120 <- subset(li6800_merged, id == "4414" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4414_120)

aci_coefs[129,] <- c(id = 4414, spp = "Tri", plot = 5, subplot = 22,
                     doy = 120, t(coef(tri_4414_120)),
                     leaf_length_cm = 14.7, stem_length_cm = NA)

# 5229
tri_5229_120 <- subset(li6800_merged, id == "5229" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5229_120)

aci_coefs[130,] <- c(id = 5229, spp = "Tri", plot = 5, subplot = 22,
                     doy = 120, t(coef(tri_5229_120)),
                     leaf_length_cm = 14.0, stem_length_cm = NA)

# 4990
tri_4990_120 <- subset(li6800_merged, id == "5229" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4990_120)

aci_coefs[131,] <- c(id = 4990, spp = "Tri", plot = 5, subplot = 22,
                     doy = 120, t(coef(tri_4990_120)),
                     leaf_length_cm = 18.1, stem_length_cm = NA)

# 482
tri_482_120 <- subset(li6800_merged, id == "482" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_482_120)

aci_coefs[132,] <- c(id = 482, spp = "Tri", plot = 5, subplot = 22,
                     doy = 120, t(coef(tri_482_120)),
                     leaf_length_cm = 14.5, stem_length_cm = NA)

# 3004
tri_3004_120 <- subset(li6800_merged, id == "3004" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3004_120)

aci_coefs[133,] <- c(id = 3004, spp = "Tri", plot = 5, subplot = 21,
                     doy = 120, t(coef(tri_3004_120)),
                     leaf_length_cm = 18.1, stem_length_cm = NA)

# 2508
tri_2508_120 <- subset(li6800_merged, id == "2508" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2508_120)

aci_coefs[134,] <- c(id = 2508, spp = "Tri", plot = 5, subplot = 21,
                     doy = 120, t(coef(tri_2508_120)),
                     leaf_length_cm = 13.6, stem_length_cm = NA)

# 4576
tri_4576_120 <- subset(li6800_merged, id == "4576" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4576_120)

aci_coefs[135,] <- c(id = 4576, spp = "Tri", plot = 5, subplot = 21,
                     doy = 120, t(coef(tri_4576_120)),
                     leaf_length_cm = 14.6, stem_length_cm = NA)

# 3077
tri_3077_120 <- subset(li6800_merged, id == "3077" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3077_120)

aci_coefs[136,] <- c(id = 3077, spp = "Tri", plot = 5, subplot = 21,
                     doy = 120, t(coef(tri_3077_120)),
                     leaf_length_cm = 11.5, stem_length_cm = NA)


# 5052
mai_5052_120 <- subset(li6800_merged, id == "5052" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5052_120)

aci_coefs[137,] <- c(id = 5052, spp = "Mai", plot = 5, subplot = 36,
                     doy = 120, t(coef(mai_5052_120)),
                     leaf_length_cm = 16.9, stem_length_cm = 54.9)

# 2637
mai_2637_120 <- subset(li6800_merged, id == "2637" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2637_120)

aci_coefs[138,] <- c(id = 2637, spp = "Mai", plot = 5, subplot = 36,
                     doy = 120, t(coef(mai_2637_120)),
                     leaf_length_cm = 15.8, stem_length_cm = 40.8)

# 2662
mai_2662_120 <- subset(li6800_merged, id == "2662" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2662_120)

aci_coefs[139,] <- c(id = 2662, spp = "Mai", plot = 5, subplot = 36,
                     doy = 120, t(coef(mai_2616_120)),
                     leaf_length_cm = 17.5, stem_length_cm = 49.4)

# 2616
mai_2616_120 <- subset(li6800_merged, id == "2616" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2616_120)

aci_coefs[140,] <- c(id = 2616, spp = "Mai", plot = 5, subplot = 36,
                     doy = 120, t(coef(mai_2616_120)),
                     leaf_length_cm = 21.1, stem_length_cm = 74.0)

# 6888
tri_6888_120 <- subset(li6800_merged, id == "6888" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6888_120)

aci_coefs[141,] <- c(id = 6888, spp = "Tri", plot = 5, subplot = 35,
                     doy = 120, t(coef(tri_6888_120)),
                     leaf_length_cm = 14.6, stem_length_cm = NA)

# 6875
tri_6875_120 <- subset(li6800_merged, id == "6875" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 450)
# plot(tri_6875_120)

aci_coefs[142,] <- c(id = 6875, spp = "Tri", plot = 5, subplot = 29,
                     doy = 120, t(coef(tri_6875_120)),
                     leaf_length_cm = 14.0, stem_length_cm = NA)

# 4444
tri_4444_120 <- subset(li6800_merged, id == "4444" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4444_120)

aci_coefs[143,] <- c(id = 4444, spp = "Tri", plot = 5, subplot = 33,
                     doy = 120, t(coef(tri_4444_120)),
                     leaf_length_cm = 12.8, stem_length_cm = NA)

# 2988
tri_2988_120 <- subset(li6800_merged, id == "2988" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2988_120)

aci_coefs[144,] <- c(id = 2988, spp = "Tri", plot = 5, subplot = 33,
                     doy = 120, t(coef(tri_2988_120)),
                     leaf_length_cm = 16.5, stem_length_cm = NA)

# 3829
tri_3829_120 <- subset(li6800_merged, id == "3829" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3829_120)

aci_coefs[145,] <- c(id = 3829, spp = "Tri", plot = 5, subplot = 33,
                     doy = 120, t(coef(tri_3829_120)),
                     leaf_length_cm = 13.6, stem_length_cm = NA)

# 5877
tri_5877_120 <- subset(li6800_merged, id == "5877" & doy == 120 & (A < 3.5 | A > 9)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5877_120)

aci_coefs[146,] <- c(id = 5877, spp = "Tri", plot = 5, subplot = 27,
                     doy = 120, t(coef(tri_5877_120)),
                     leaf_length_cm = 16.5, stem_length_cm = NA)

# 4109
tri_4109_120 <- subset(li6800_merged, id == "4109" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4109_120)

aci_coefs[147,] <- c(id = 4109, spp = "Tri", plot = 5, subplot = 32,
                     doy = 120, t(coef(tri_4109_120)),
                     leaf_length_cm = 13.2, stem_length_cm = NA)

# 2573
tri_2573_120 <- subset(li6800_merged, id == "2573" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2573_120)

aci_coefs[148,] <- c(id = 2573, spp = "Tri", plot = 5, subplot = 31,
                     doy = 120, t(coef(tri_2573_120)),
                     leaf_length_cm = 15.7, stem_length_cm = NA)

# 4265
tri_4265_120 <- subset(li6800_merged, id == "4256" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4265_120)

aci_coefs[149,] <- c(id = 4265, spp = "Tri", plot = 5, subplot = 31,
                     doy = 120, t(coef(tri_4265_120)),
                     leaf_length_cm = 17.6, stem_length_cm = NA)

# 2547
tri_2547_120 <- subset(li6800_merged, id == "2547" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2547_120)

aci_coefs[150,] <- c(id = 2547, spp = "Tri", plot = 5, subplot = 19,
                     doy = 120, t(coef(tri_2547_120)),
                     leaf_length_cm = 16.0, stem_length_cm = NA)

# 5228
tri_5228_120 <- subset(li6800_merged, id == "5228" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5228_120)

aci_coefs[151,] <- c(id = 5228, spp = "Tri", plot = 5, subplot = 19,
                     doy = 120, t(coef(tri_5228_120)),
                     leaf_length_cm = 14.9, stem_length_cm = NA)

# 5115
tri_5115_120 <- subset(li6800_merged, id == "5115" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5115_120)

aci_coefs[152,] <- c(id = 5115, spp = "Tri", plot = 5, subplot = 19,
                     doy = 120, t(coef(tri_5115_120)),
                     leaf_length_cm = 13.2, stem_length_cm = NA)

# 1730
tri_1730_120 <- subset(li6800_merged, id == "1730" & doy == 120) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_1730_120)

aci_coefs[153,] <- c(id = 1730, spp = "Tri", plot = 5, subplot = 19,
                     doy = 120, t(coef(tri_1730_120)),
                     leaf_length_cm = 14.7, stem_length_cm = NA)

##stopped here ##


#####################################################################
# 5/1/24: plot 6
#####################################################################
# 3136
mai_3136_122 <- subset(li6800_merged, id == "3136" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_3136_122)

aci_coefs[154,] <- c(id = 3136, spp = "Mai", plot = 6, subplot = 32,
                     doy = 122, t(coef(mai_3136_122)),
                     leaf_length_cm = 11.9, stem_length_cm = 27.3)

# 5178
mai_5178_122 <- subset(li6800_merged, id == "5178" & doy == 122 & (A < 13 | A > 15)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5178_122)

aci_coefs[155,] <- c(id = 5178, spp = "Mai", plot = 6, subplot = 32,
                     doy = 122, t(coef(mai_5178_122)),
                     leaf_length_cm = 17.3, stem_length_cm = 61.3)

# 5600
mai_5600_122 <- subset(li6800_merged, id == "5600" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5600_122)

aci_coefs[156,] <- c(id = 5600, spp = "Mai", plot = 6, subplot = 32,
                     doy = 122, t(coef(mai_5600_122)),
                     leaf_length_cm = 21.6, stem_length_cm = 73.2)

# 2608
mai_2608_122 <- subset(li6800_merged, id == "2608" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 200)
# plot(mai_2608_122)

aci_coefs[157,] <- c(id = 2608, spp = "Mai", plot = 6, subplot = 25,
                     doy = 122, t(coef(mai_2608_122)),
                     leaf_length_cm = 16.2, stem_length_cm = 48.1)

# 6495
mai_6495_122 <- subset(li6800_merged, id == "6495_c" & doy == 122 & (A < 18 | A > 21)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_6495_122)

aci_coefs[158,] <- c(id = 6495, spp = "Mai", plot = 6, subplot = 25,
                     doy = 122, t(coef(mai_6495_122)),
                     leaf_length_cm = 19.9, stem_length_cm = 61.4)

# 6462
mai_6462_122 <- subset(li6800_merged, id == "6462_b" & doy == 122 & (A < 9.5 | A > 12)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_6462_122)

aci_coefs[159,] <- c(id = 6462, spp = "Mai", plot = 6, subplot = 36,
                     doy = 122, t(coef(mai_6462_122)),
                     leaf_length_cm = 15.5, stem_length_cm = 41.5)

# 6448
mai_6448_122 <- subset(li6800_merged, id == "6448" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_6448_122)

aci_coefs[160,] <- c(id = 6448, spp = "Mai", plot = 6, subplot = 36,
                     doy = 122, t(coef(mai_6448_122)),
                     leaf_length_cm = 17.0, stem_length_cm = 80.1)

# 5229
tri_5229_122 <- subset(li6800_merged, id == "5229" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5229_122)

aci_coefs[161,] <- c(id = 5229, spp = "Tri", plot = 6, subplot = 1,
                     doy = 122, t(coef(tri_5229_122)),
                     leaf_length_cm = 16.2, stem_length_cm = NA)

# 4766
tri_4766_122 <- subset(li6800_merged, id == "4766" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4766_122)

aci_coefs[162,] <- c(id = 4766, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_4766_122)),
                     leaf_length_cm = 10.4, stem_length_cm = NA)

# 762
tri_762_122 <- subset(li6800_merged, id == "762" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_762_122)

aci_coefs[163,] <- c(id = 762, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_762_122)),
                     leaf_length_cm = 14.1, stem_length_cm = NA)

# 4373
tri_4373_122 <- subset(li6800_merged, id == "4373" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4373_122)

aci_coefs[164,] <- c(id = 4373, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_4373_122)),
                     leaf_length_cm = 13.6, stem_length_cm = NA)

# 4777
tri_4777_122 <- subset(li6800_merged, id == "4777" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_4777_122)

aci_coefs[165,] <- c(id = 4777, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_4777_122)),
                     leaf_length_cm = 13.9, stem_length_cm = NA)

# 5865
tri_5865_122 <- subset(li6800_merged, id == "5865" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5865_122)

aci_coefs[166,] <- c(id = 5865, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_5865_122)),
                     leaf_length_cm = 13.1, stem_length_cm = NA)

# 5267
tri_5267_122 <- subset(li6800_merged, id == "5267" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5267_122)

aci_coefs[167,] <- c(id = 5267, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_5267_122)),
                     leaf_length_cm = 15.4, stem_length_cm = NA)

# 684
tri_684_122 <- subset(li6800_merged, id == "684" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_684_122)

aci_coefs[168,] <- c(id = 684, spp = "Tri", plot = 6, subplot = 2,
                     doy = 122, t(coef(tri_684_122)),
                     leaf_length_cm = 15.9, stem_length_cm = NA)

# 5742
tri_5742_122 <- subset(li6800_merged, id == "5742" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5742_122)

aci_coefs[169,] <- c(id = 5742, spp = "Tri", plot = 6, subplot = 3,
                     doy = 122, t(coef(tri_5742_122)),
                     leaf_length_cm = 15.3, stem_length_cm = NA)

# 5641
tri_5641_122 <- subset(li6800_merged, id == "5641" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5641_122)

aci_coefs[170,] <- c(id = 5641, spp = "Tri", plot = 6, subplot = 4,
                     doy = 122, t(coef(tri_5641_122)),
                     leaf_length_cm = 13.7, stem_length_cm = NA)

# 3305
tri_3305_122 <- subset(li6800_merged, id == "3305" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_3305_122)

aci_coefs[171,] <- c(id = 3305, spp = "Tri", plot = 6, subplot = 5,
                     doy = 122, t(coef(tri_3305_122)),
                     leaf_length_cm = 11.6, stem_length_cm = NA)

# 5619
tri_5619_122 <- subset(li6800_merged, id == "5619" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5619_122)

aci_coefs[172,] <- c(id = 5619, spp = "Tri", plot = 6, subplot = 5,
                     doy = 122, t(coef(tri_5619_122)),
                     leaf_length_cm = 13.2, stem_length_cm = NA)

# striped2
tri_striped2_122 <- subset(li6800_merged, id == "striped2" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_striped2_122)

aci_coefs[173,] <- c(id = "striped2", spp = "Tri", plot = 6, subplot = 5,
                     doy = 122, t(coef(tri_striped2_122)),
                     leaf_length_cm = 15.5, stem_length_cm = NA)

# striped1
tri_striped1_122 <- subset(li6800_merged, id == "striped1" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_striped1_122)

aci_coefs[174,] <- c(id = "striped1", spp = "Tri", plot = 6, subplot = 5,
                     doy = 122, t(coef(tri_striped1_122)),
                     leaf_length_cm = 13.7, stem_length_cm = NA)

# 978
tri_978_122 <- subset(li6800_merged, id == "978" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_978_122)

aci_coefs[175,] <- c(id = "978", spp = "Tri", plot = 6, subplot = 6,
                     doy = 122, t(coef(tri_978_122)),
                     leaf_length_cm = 14.3, stem_length_cm = NA)

# 5626
tri_5626_122 <- subset(li6800_merged, id == "5626" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5626_122)

aci_coefs[176,] <- c(id = "5626", spp = "Tri", plot = 6, subplot = 6,
                     doy = 122, t(coef(tri_5626_122)),
                     leaf_length_cm = 15.9, stem_length_cm = NA)

# 5403
tri_5403_122 <- subset(li6800_merged, id == "5403" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5403_122)

aci_coefs[177,] <- c(id = "5403", spp = "Tri", plot = 6, subplot = 6,
                     doy = 122, t(coef(tri_5403_122)),
                     leaf_length_cm = 16.1, stem_length_cm = NA)

# 6879
tri_6879_122 <- subset(li6800_merged, id == "6879" & doy == 122) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_6879_122)

aci_coefs[178,] <- c(id = "6879", spp = "Tri", plot = 6, subplot = 36,
                     doy = 122, t(coef(tri_6879_122)),
                     leaf_length_cm = 13.4, stem_length_cm = NA)

write.csv(aci_coefs, "../data/TT24_curve_fits.csv", row.names = F)
#####################################################################
# 5/2/24: plot 6
#####################################################################

# 5579
mai_5579_123 <- subset(li6800_merged, id == "5579" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5579_123)

aci_coefs[179,] <- c(id = "5579", spp = "Mai", plot = 6, subplot = 25,
                     doy = 123, t(coef(mai_5579_123)),
                     leaf_length_cm = 20.6, stem_length_cm = 67.2)

# 5024
mai_5024_123 <- subset(li6800_merged, id == "5024" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5024_123)

aci_coefs[180,] <- c(id = "5024", spp = "Mai", plot = 6, subplot = 27,
                     doy = 123, t(coef(mai_5024_123)),
                     leaf_length_cm = 12.2, stem_length_cm = 19.1)

# 5060
tri_5060_123 <- subset(li6800_merged, id == "5060" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5060_123)

aci_coefs[181,] <- c(id = "5060", spp = "Tri", plot = 6, subplot = 27,
                     doy = 123, t(coef(tri_5060_123)),
                     leaf_length_cm = 15.8, stem_length_cm = NA)

# 6483
mai_6483_123 <- subset(li6800_merged, id == "6483" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_6483_123)

aci_coefs[182,] <- c(id = "6483", spp = "Mai", plot = 6, subplot = 19,
                     doy = 123, t(coef(mai_6483_123)),
                     leaf_length_cm = 16.2, stem_length_cm = 54.3)

# 5783
tri_5783_123 <- subset(li6800_merged, id == "5783" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5783_123)

aci_coefs[183,] <- c(id = "5783", spp = "Tri", plot = 6, subplot = 20,
                     doy = 123, t(coef(tri_5783_123)),
                     leaf_length_cm = 17.1, stem_length_cm = NA)

# 5506
mai_5506_123 <- subset(li6800_merged, id == "5506" & doy == 123 & (A < 6 | A > 12)) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5506_123)

aci_coefs[184,] <- c(id = "5506", spp = "Mai", plot = 6, subplot = 21,
                     doy = 123, t(coef(mai_5506_123)),
                     leaf_length_cm = 13.8, stem_length_cm = 34.5)

# 1941
mai_1941_123 <- subset(li6800_merged, id == "1941" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_1941_123)

aci_coefs[185,] <- c(id = "1941", spp = "Mai", plot = 6, subplot = 19,
                     doy = 123, t(coef(mai_1941_123)),
                     leaf_length_cm = 15.9, stem_length_cm = 50.6)

# 2726
mai_2726_123 <- subset(li6800_merged, id == "2726" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2726_123)

aci_coefs[186,] <- c(id = "2726", spp = "Mai", plot = 6, subplot = 30,
                     doy = 123, t(coef(mai_2726_123)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 5504
mai_5504_123 <- subset(li6800_merged, id == "5504" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_5504_123)

aci_coefs[187,] <- c(id = "5504", spp = "Mai", plot = 6, subplot = 36,
                     doy = 123, t(coef(mai_5504_123)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 2743
mai_2743_123 <- subset(li6800_merged, id == "2743" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_2743_123)

aci_coefs[188,] <- c(id = "2743", spp = "Mai", plot = 6, subplot = 36,
                     doy = 123, t(coef(mai_2743_123)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 1467
mai_1467_123 <- subset(li6800_merged, id == "1467" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(mai_1467_123)

aci_coefs[189,] <- c(id = "1467", spp = "Mai", plot = 6, subplot = 36,
                     doy = 123, t(coef(mai_1467_123)),
                     leaf_length_cm = NA, stem_length_cm = NA)

# 2276
tri_2276_123 <- subset(li6800_merged, id == "2276" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_2276_123)

aci_coefs[190,] <- c(id = "2276", spp = "Tri", plot = 6, subplot = 22,
                     doy = 123, t(coef(tri_2276_123)),
                     leaf_length_cm = 11.9, stem_length_cm = NA)

# 5031
tri_5031_123 <- subset(li6800_merged, id == "5031" & doy == 123) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
# plot(tri_5031_123)

aci_coefs[191,] <- c(id = "5031", spp = "Tri", plot = 6, subplot = 22,
                     doy = 123, t(coef(tri_5031_123)),
                     leaf_length_cm = 11.0, stem_length_cm = NA)

## stopped here 10/23/24 ##

#####################################################################
# 5/6/24: plot 3
#####################################################################

# 511
mai_511_127 <- subset(li6800_merged, id == "511" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_511_127)

aci_coefs[192,] <- c(id = "511", spp = "Mai", plot = 3, subplot = 30,
                     doy = 127, t(coef(mai_511_127)))

# 1021
mai_1021_127 <- subset(li6800_merged, id == "1021" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1021_127)

aci_coefs[193,] <- c(id = "1021", spp = "Mai", plot = 3, subplot = 24,
                     doy = 127, t(coef(mai_1021_127)))

# 543
mai_543_127 <- subset(li6800_merged, id == "543" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_543_127)

aci_coefs[194,] <- c(id = "543", spp = "Mai", plot = 3, subplot = 30,
                     doy = 127, t(coef(mai_543_127)))

# TT24_202
mai_TT24_202_127 <- subset(li6800_merged, id == "TT24_202" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_202_127)

aci_coefs[195,] <- c(id = "TT24_202", spp = "Mai", plot = 3, subplot = 36,
                     doy = 127, t(coef(mai_TT24_202_127)))

# 2268
mai_2268_127 <- subset(li6800_merged, id == "2268" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2268_127)

aci_coefs[196,] <- c(id = "2268", spp = "Mai", plot = 3, subplot = 18,
                     doy = 127, t(coef(mai_2268_127)))

# 5069
mai_5069_127 <- subset(li6800_merged, id == "5069" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5069_127)

aci_coefs[197,] <- c(id = "5069", spp = "Mai", plot = 3, subplot = 24,
                     doy = 127, t(coef(mai_5069_127)))

# 5030
mai_5030_127 <- subset(li6800_merged, id == "5030" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5030_127)

aci_coefs[198,] <- c(id = "5030", spp = "Mai", plot = 3, subplot = 24,
                     doy = 127, t(coef(mai_5030_127)))

# 6876
mai_6876_127 <- subset(li6800_merged, id == "6876" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_6876_127)

aci_coefs[199,] <- c(id = "6876", spp = "Mai", plot = 3, subplot = 18,
                     doy = 127, t(coef(mai_6876_127)))

# 129
mai_129_127 <- subset(li6800_merged, id == "129" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_129_127)

aci_coefs[200,] <- c(id = "129", spp = "Mai", plot = 3, subplot = 18,
                     doy = 127, t(coef(mai_129_127)))

# 3643
mai_3643_127 <- subset(li6800_merged, id == "3643" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3643_127)

aci_coefs[201,] <- c(id = "3643", spp = "Mai", plot = 3, subplot = 12,
                     doy = 127, t(coef(mai_3643_127)))

# 2416
tri_2416_127 <- subset(li6800_merged, id == "2416" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2416_127)

aci_coefs[202,] <- c(id = "2416", spp = "Tri", plot = 3, subplot = 17,
                     doy = 127, t(coef(tri_2416_127)))

# TT24_201
mai_TT24_201_127 <- subset(li6800_merged, id == "TT24_201" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_201_127)

aci_coefs[203,] <- c(id = "TT24_201", spp = "Mai", plot = 3, subplot = 35,
                     doy = 127, t(coef(mai_TT24_201_127)))

# TT24_203
mai_TT24_203_127 <- subset(li6800_merged, id == "TT24_203" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_203_127)

aci_coefs[204,] <- c(id = "TT24_203", spp = "Mai", plot = 3, subplot = 34,
                     doy = 127, t(coef(mai_TT24_203_127)))

# TT24_204
mai_TT24_204_127 <- subset(li6800_merged, id == "TT24_204" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_204_127)

aci_coefs[205,] <- c(id = "TT24_204", spp = "Mai", plot = 3, subplot = 35,
                     doy = 127, t(coef(mai_TT24_204_127)))

# TT24_103
mai_TT24_103_127 <- subset(li6800_merged, id == "TT24_103" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_103_127)

aci_coefs[206,] <- c(id = "TT24_103", spp = "Mai", plot = 3, subplot = 32,
                     doy = 127, t(coef(mai_TT24_103_127)))

# TT24_101
mai_TT24_101_127 <- subset(li6800_merged, id == "TT24_101" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_101_127)

aci_coefs[207,] <- c(id = "TT24_101", spp = "Mai", plot = 3, subplot = 33,
                     doy = 127, t(coef(mai_TT24_101_127)))

# TT24_102
mai_TT24_102_127 <- subset(li6800_merged, id == "TT24_102" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_102_127)

aci_coefs[208,] <- c(id = "TT24_102", spp = "Mai", plot = 3, subplot = 32,
                     doy = 127, t(coef(mai_TT24_102_127)))

# 452
tri_452_127 <- subset(li6800_merged, id == "452" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_452_127)

aci_coefs[209,] <- c(id = 452, spp = "Tri", plot = 3, subplot = 22,
                     doy = 127, t(coef(tri_452_127)))

# 54
mai_54_127 <- subset(li6800_merged, id == "54" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_54_127)

aci_coefs[210,] <- c(id = 54, spp = "Mai", plot = 3, subplot = 16,
                     doy = 127, t(coef(mai_54_127)))

# 916
tri_916_127 <- subset(li6800_merged, id == "916" & doy == 127) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_916_127)

aci_coefs[211,] <- c(id = 916, spp = "Tri", plot = 3, subplot = 16,
                     doy = 127, t(coef(tri_916_127)))

#####################################################################
# 5/7/24: plot 6
#####################################################################

# TT24_104
mai_TT24_104_128 <- subset(li6800_merged, id == "TT24_104" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_104_128)

aci_coefs[212,] <- c(id = "TT24_104", spp = "Mai", plot = 6, subplot = 33,
                     doy = 128, t(coef(mai_TT24_104_128)))

# TT24_105
mai_TT24_105_128 <- subset(li6800_merged, id == "TT24_105" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_105_128)

aci_coefs[213,] <- c(id = "TT24_105", spp = "Mai", plot = 6, subplot = 27,
                     doy = 128, t(coef(mai_TT24_105_128)))

# TT24_106
mai_TT24_106_128 <- subset(li6800_merged, id == "TT24_106" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_106_128)

aci_coefs[214,] <- c(id = "TT24_106", spp = "Mai", plot = 6, subplot = 27,
                     doy = 128, t(coef(mai_TT24_106_128)))

# 2801
mai_2801_128 <- subset(li6800_merged, id == "2801" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2801_128)

aci_coefs[215,] <- c(id = "2801", spp = "Mai", plot = 6, subplot = 20,
                     doy = 128, t(coef(mai_2801_128)))

# 2692
mai_2692_128 <- subset(li6800_merged, id == "2692" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2692_128)

aci_coefs[216,] <- c(id = "2692", spp = "Mai", plot = 6, subplot = 21,
                     doy = 128, t(coef(mai_2692_128)))

# 2639
mai_2639_128 <- subset(li6800_merged, id == "2639" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2639_128)

aci_coefs[217,] <- c(id = "2639", spp = "Mai", plot = 6, subplot = 21,
                     doy = 128, t(coef(mai_2639_128)))

# 2799
mai_2799_128 <- subset(li6800_merged, id == "279" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2799_128)

aci_coefs[218,] <- c(id = "2799", spp = "Mai", plot = 6, subplot = 14,
                     doy = 128, t(coef(mai_2799_128)))

# 5714
tri_5714_128 <- subset(li6800_merged, id == "5714" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5714_128)

aci_coefs[219,] <- c(id = "5714", spp = "Tri", plot = 6, subplot = 8,
                     doy = 128, t(coef(tri_5714_128)))

# 5596
mai_5596_128 <- subset(li6800_merged, id == "5596" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5596_128)

aci_coefs[220,] <- c(id = "5596", spp = "Mai", plot = 6, subplot = 8,
                     doy = 128, t(coef(mai_5596_128)))

# 5852
tri_5852_128 <- subset(li6800_merged, id == "5852" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5852_128)

aci_coefs[221,] <- c(id = "5852", spp = "Tri", plot = 6, subplot = 15,
                     doy = 128, t(coef(tri_5852_128)))

# 5607
mai_5607_128 <- subset(li6800_merged, id == "5607" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5607_128)

aci_coefs[222,] <- c(id = "5607", spp = "Mai", plot = 6, subplot = 9,
                     doy = 128, t(coef(mai_5607_128)))

# 5664
mai_5664_128 <- subset(li6800_merged, id == "5664" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5664_128)

aci_coefs[223,] <- c(id = "5664", spp = "Mai", plot = 6, subplot = 9,
                     doy = 128, t(coef(mai_5664_128)))

# 5122
mai_5122_128 <- subset(li6800_merged, id == "5122" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5122_128)

aci_coefs[224,] <- c(id = "5122", spp = "Mai", plot = 6, subplot = 2,
                     doy = 128, t(coef(mai_5122_128)))

# 2770
tri_2770_128 <- subset(li6800_merged, id == "2770" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2770_128)

aci_coefs[225,] <- c(id = "2770", spp = "Mai", plot = 6, subplot = 9,
                     doy = 128, t(coef(tri_2770_128)))

# TT24_205
mai_TT24_205_128 <- subset(li6800_merged, id == "TT24_205" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_205_128)

aci_coefs[226,] <- c(id = "TT24_205", spp = "Mai", plot = 6, subplot = 34,
                     doy = 128, t(coef(mai_TT24_205_128)))

# TT24_206
mai_TT24_206_128 <- subset(li6800_merged, id == "TT24_206" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_206_128)

aci_coefs[227,] <- c(id = "TT24_206", spp = "Mai", plot = 6, subplot = 22,
                     doy = 128, t(coef(mai_TT24_206_128)))

# 3175
mai_3175_128 <- subset(li6800_merged, id == "3175" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3175_128)

aci_coefs[228,] <- c(id = "3175", spp = "Mai", plot = 6, subplot = 22,
                     doy = 128, t(coef(mai_3175_128)))

# 5567
mai_5567_128 <- subset(li6800_merged, id == "5567" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5567_128)

aci_coefs[229,] <- c(id = "5567", spp = "Mai", plot = 6, subplot = 22,
                     doy = 128, t(coef(mai_5567_128)))

# 6463
mai_6463_128 <- subset(li6800_merged, id == "6463" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_6463_128)

aci_coefs[230,] <- c(id = "6463", spp = "Mai", plot = 6, subplot = NA,
                     doy = 128, t(coef(mai_6463_128)))

# 4511
tri_4511_128 <- subset(li6800_merged, id == "4511" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4511_128)

aci_coefs[231,] <- c(id = "4511", spp = "Tri", plot = 6, subplot = 16,
                     doy = 128, t(coef(tri_4511_128)))

# 4216 (note ID is correct - switched order in user ID)
mai_4216_128 <- subset(li6800_merged, id == "3438" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4216_128)

aci_coefs[232,] <- c(id = "6463", spp = "Mai", plot = 6, subplot = 11,
                     doy = 128, t(coef(mai_4216_128)))

# 3438 (note ID is correct - switched order in user ID)
mai_3438_128 <- subset(li6800_merged, id == "4216" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3438_128)

aci_coefs[233,] <- c(id = "3438", spp = "Mai", plot = 6, subplot = 11,
                     doy = 128, t(coef(mai_3438_128)))

# 5025
mai_5025_128 <- subset(li6800_merged, id == "5025" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5025_128)

aci_coefs[234,] <- c(id = "5025", spp = "Mai", plot = 6, subplot = 11,
                     doy = 128, t(coef(mai_5025_128)))

# 4582
tri_4582_128 <- subset(li6800_merged, id == "4582" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4582_128)

aci_coefs[235,] <- c(id = "4582", spp = "Tri", plot = 6, subplot = 11,
                     doy = 128, t(coef(tri_4582_128)))

# 2980
tri_2980_128 <- subset(li6800_merged, id == "2980" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2980_128)

aci_coefs[236,] <- c(id = "2980", spp = "Tri", plot = 6, subplot = 11,
                     doy = 128, t(coef(tri_2980_128)))

# 3371
tri_3371_128 <- subset(li6800_merged, id == "3371" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3371_128)

aci_coefs[237,] <- c(id = "3371", spp = "Tri", plot = 6, subplot = 11,
                     doy = 128, t(coef(tri_3371_128)))

# 3379
tri_3379_128 <- subset(li6800_merged, id == "3379" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3379_128)

aci_coefs[238,] <- c(id = "3379", spp = "Tri", plot = 6, subplot = 11,
                     doy = 128, t(coef(tri_3379_128)))

# 4505
tri_4505_128 <- subset(li6800_merged, id == "4505" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4505_128)

aci_coefs[239,] <- c(id = "4505", spp = "Tri", plot = 6, subplot = 11,
                     doy = 128, t(coef(tri_4505_128)))

# 4745
tri_4745_128 <- subset(li6800_merged, id == "4745" & doy == 128) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4745_128)

aci_coefs[240,] <- c(id = "4745", spp = "Tri", plot = 6, subplot = 12,
                     doy = 128, t(coef(tri_4745_128)))

#####################################################################
# 5/8/24: plot 5
#####################################################################
# 4000
tri_4000_129 <- subset(li6800_merged, id == "4000" & doy == 129) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4000_129)

aci_coefs[241,] <- c(id = "4000", spp = "Tri", plot = 5, subplot = 20,
                     doy = 129, t(coef(tri_4000_129)))

#####################################################################
# 5/13/24: plot 5
#####################################################################
# 7120
tri_7120_134 <- subset(li6800_merged, id == "7120" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 900)
plot(tri_7120_134)
summary(tri_7120_134)

aci_coefs[242,] <- c(id = "7120", spp = "Tri", plot = 5, subplot = 17,
                     doy = 134, t(coef(tri_7120_134)))

# 2803
mai_2803_134 <- subset(li6800_merged, id == "2803" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2803_134)

aci_coefs[243,] <- c(id = "2803", spp = "Mai", plot = 5, subplot = 18,
                     doy = 134, t(coef(mai_2803_134)))

# 1499
mai_1499_134 <- subset(li6800_merged, id == "1499" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1499_134)

aci_coefs[244,] <- c(id = "1499", spp = "Mai", plot = 5, subplot = 17,
                     doy = 134, t(coef(mai_1499_134)))

# 1432
mai_1432_134 <- subset(li6800_merged, id == "1432" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1432_134)

aci_coefs[245,] <- c(id = "1432", spp = "Mai", plot = 5, subplot = 17,
                     doy = 134, t(coef(mai_1432_134)))

# 2912
tri_2912_134 <- subset(li6800_merged, id == "2912" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2912_134)

aci_coefs[246,] <- c(id = "2912", spp = "Tri", plot = 5, subplot = 11,
                     doy = 134, t(coef(tri_2912_134)))

# 614
tri_614_134 <- subset(li6800_merged, id == "614" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_614_134)

aci_coefs[247,] <- c(id = "614", spp = "Tri", plot = 5, subplot = 11,
                     doy = 134, t(coef(tri_614_134)))

# 6881
tri_6881_134 <- subset(li6800_merged, id == "6881" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6881_134)

aci_coefs[248,] <- c(id = "6881", spp = "Tri", plot = 5, subplot = 10,
                     doy = 134, t(coef(tri_6881_134)))

# 1374
tri_1374_134 <- subset(li6800_merged, id == "1374" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1374_134)

aci_coefs[249,] <- c(id = "1374", spp = "Tri", plot = 5, subplot = 10,
                     doy = 134, t(coef(tri_1374_134)))

# 1826
tri_1826_134 <- subset(li6800_merged, id == "1826" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1826_134)

aci_coefs[250,] <- c(id = "1826", spp = "Tri", plot = 5, subplot = 10,
                     doy = 134, t(coef(tri_1826_134)))

# 4105
tri_4105_134 <- subset(li6800_merged, id == "4105" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4105_134)

aci_coefs[251,] <- c(id = "4105", spp = "Tri", plot = 5, subplot = 10,
                     doy = 134, t(coef(tri_4105_134)))

# 5381
tri_5381_134 <- subset(li6800_merged, id == "5381" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5381_134)

aci_coefs[252,] <- c(id = "5381", spp = "Tri", plot = 5, subplot = 10,
                     doy = 134, t(coef(tri_5381_134)))

# 4574
tri_4574_134 <- subset(li6800_merged, id == "4574" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4574_134)

aci_coefs[253,] <- c(id = "4574", spp = "Tri", plot = 5, subplot = 10,
                     doy = 134, t(coef(tri_4574_134)))

# flag8
tri_flag8_134 <- subset(li6800_merged, id == "fla8" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag8_134)

aci_coefs[254,] <- c(id = "flag8", spp = "Tri", plot = 5, subplot = 16,
                     doy = 134, t(coef(tri_flag8_134)))

# 2563
tri_2563_134 <- subset(li6800_merged, id == "2563" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2563_134)

aci_coefs[255,] <- c(id = "2563", spp = "Tri", plot = 5, subplot = 16,
                     doy = 134, t(coef(tri_2563_134)))

# 4177
tri_4177_134 <- subset(li6800_merged, id == "4177" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4177_134)

aci_coefs[256,] <- c(id = "4177", spp = "Tri", plot = 5, subplot = 14,
                     doy = 134, t(coef(tri_4177_134)))

# 3
tri_flag3_134 <- subset(li6800_merged, id == "3" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag3_134)

aci_coefs[257,] <- c(id = "flag3_tri", spp = "Tri", plot = 5, subplot = 14,
                     doy = 134, t(coef(tri_flag3_134)))

# 1
tri_flag1_134 <- subset(li6800_merged, id == "1" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag1_134)

aci_coefs[258,] <- c(id = "flag1_tri", spp = "Tri", plot = 5, subplot = 14,
                     doy = 134, t(coef(tri_flag1_134)))

# 2
tri_flag2_134 <- subset(li6800_merged, id == "2" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag2_134)

aci_coefs[259,] <- c(id = "flag2_tri", spp = "Tri", plot = 5, subplot = 14,
                     doy = 134, t(coef(tri_flag2_134)))

# 4149
tri_4149_134 <- subset(li6800_merged, id == "4149" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4149_134)

aci_coefs[260,] <- c(id = "4149", spp = "Tri", plot = 5, subplot = 13,
                     doy = 134, t(coef(tri_4149_134)))

# 1795
tri_1795_134 <- subset(li6800_merged, id == "1795" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1795_134)

aci_coefs[261,] <- c(id = "1795", spp = "Tri", plot = 5, subplot = 14,
                     doy = 134, t(coef(tri_1795_134)))

# 6887
tri_6887_134 <- subset(li6800_merged, id == "6887" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6887_134)

aci_coefs[262,] <- c(id = "6887", spp = "Tri", plot = 5, subplot = 8,
                     doy = 134, t(coef(tri_6887_134)))

# 2
mai_flag2_134 <- subset(li6800_merged, id == "2mai" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_flag2_134)

aci_coefs[263,] <- c(id = "flag2_mai", spp = "Mai", plot = 5, subplot = 3,
                     doy = 134, t(coef(mai_flag2_134)))

# 631
mai_631_134 <- subset(li6800_merged, id == "631" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_631_134)

aci_coefs[264,] <- c(id = "631", spp = "Mai", plot = 5, subplot = 3,
                     doy = 134, t(coef(mai_631_134)))

# 3960
mai_3960_134 <- subset(li6800_merged, id == "3960" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3960_134)

aci_coefs[265,] <- c(id = "3960", spp = "Mai", plot = 5, subplot = 3,
                     doy = 134, t(coef(mai_3960_134)))

# 3993
mai_3993_134 <- subset(li6800_merged, id == "3993" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3993_134)

aci_coefs[266,] <- c(id = "3993", spp = "Mai", plot = 5, subplot = 3,
                     doy = 134, t(coef(mai_3993_134)))

# 86
tri_86_134 <- subset(li6800_merged, id == "86" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_86_134)

aci_coefs[267,] <- c(id = "86", spp = "Tri", plot = 5, subplot = 9,
                     doy = 134, t(coef(tri_86_134)))

# 7147
tri_7147_134 <- subset(li6800_merged, id == "7147" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_7147_134)

aci_coefs[268,] <- c(id = "7147", spp = "Tri", plot = 5, subplot = 9,
                     doy = 134, t(coef(tri_7147_134)))

# 43
tri_43_134 <- subset(li6800_merged, id == "43" & doy == 134) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_43_134)

aci_coefs[269,] <- c(id = "43", spp = "Tri", plot = 5, subplot = 9,
                     doy = 134, t(coef(tri_43_134)))


#####################################################################
# 5/14/24: plot 6
#####################################################################
unique(subset(li6800_merged, doy == 135)$id)

# 5722 (correct incorrect ID)
tri_5722_135 <- subset(li6800_merged, id == "4722" & doy == 135) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5722_135)

aci_coefs[270,] <- c(id = "5722", spp = "Tri", plot = 6, subplot = 6,
                     doy = 135, t(coef(tri_5722_135)))

# 2941
tri_2941_135 <- subset(li6800_merged, id == "2941" & doy == 135) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2941_135)

aci_coefs[271,] <- c(id = "2941", spp = "Tri", plot = 6, subplot = 6,
                     doy = 135, t(coef(tri_2941_135)))



#####################################################################
# Snapshot measurements (Anet, gsw, Ci:Ca)
#####################################################################
snapshot <- li6800_merged %>%
  group_by(id, doy) %>%
  filter(row_number() == 1) %>%
  dplyr::select(id, machine, doy, anet = A, ci = Ci, 
                ca = Ca, co2_ref = CO2_r, gsw, Tleaf) %>%
  mutate(doy = as.numeric(doy),
         id = as.character(id),
         ci.ca = ci / ca,
         iwue = anet / gsw,
         id = ifelse(id == "", NA, id)) %>%
  filter(!is.na(id)) %>%
  mutate(id = ifelse(id == "1" & doy == 114, "flag1_tri", id),
         id = ifelse(id == "2" & doy == 114, "flag2_tri", id),
         id = ifelse(id == "3" & doy == 114, "flag3_tri", id),
         id = ifelse(id == "174" & doy == 119, "179", id))

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
                                    tGrow = 20),
         jmax25 = temp_standardize(estimate = Jmax,
                                    estimate.type = "Jmax",
                                    standard.to = 25,
                                    tLeaf = Tleaf,
                                    tGrow = 20),
         rd25 = temp_standardize(estimate = Rd,
                                 estimate.type = "Rd",
                                 pft = "C3H",
                                 standard.to = 25,
                                 tLeaf = Tleaf,
                                 tGrow = 20)) %>%
  dplyr::select(id, machine, doy, spp:subplot, gm.trt, Tleaf,
                anet, ci.ca, gsw, iwue, vcmax = Vcmax, vcmax25,
                jmax = Jmax, jmax25, rd = Rd, rd25) %>%
  mutate(across(Tleaf:rd25, \(x) round(x, digits = 4))) %>%
  arrange(plot, doy)

write.csv(photo_cleaned_full,
          "../data/TT24_photo_traits_working.csv", 
          row.names = F)
