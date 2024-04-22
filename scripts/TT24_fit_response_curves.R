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

# Load custom functions for cleaning LI-6800 files,
# standardizing Vcmax/Jmax/Rd to single temperature,
# and calculating 
R.utils::sourceDirectory("../functions/")

# Create data frame containing subplots and their treatments for
# easy merge with photosynthetic trait data
gm.ambient <- c(4, 5, 6, 10, 11, 12, 16, 17, 18,
                22, 23, 24, 28, 29, 30, 34, 35, 36)

treatments <- data.frame(subplot = seq(1,36, 1)) %>%
  mutate(gm.trt = ifelse(subplot %in% gm.ambient == TRUE, "ambient", "weeded"))


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

# Read all files and merge into central dataframe
li6800_merged <- lapply(files, read.csv) %>%
  reshape::merge_all() %>%
  mutate(date = lubridate::ymd_hms(date),
         date_only = stringr::word(date, 1),
         julian_date = yday(date_only),
         Qin_cuvette = 2000) %>%
  dplyr::select(obs, time, elapsed, date, date_only, julian_date, hhmmss:Qin,
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

# Object naming scheme: Species, tag ID, julian date separated by "_"

#####################################################################
# 4/13/24: plot 3
#####################################################################

## 583
tri_583_104 <- subset(li6800_merged, id == "583" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_583_104)

aci_coefs <- data.frame(id = 583, spp = "Tri", plot = 3, subplot = 23, 
                        julian_date = 104, t(coef(tri_583_104)))


## 4934
tri_4934_104 <- subset(li6800_merged, id == "4934" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4934_104)

aci_coefs[2,] <- c(id = 4934, spp = "Tri", plot = 3, subplot = 23,
                   julian_date = 104, t(coef(tri_4934_104)))


## 452
tri_452_104 <- subset(li6800_merged, id == "452" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_452_104)

aci_coefs[3,] <- c(id = 452, spp = "Tri", plot = 3, subplot = 22,
                   julian_date = 104, t(coef(tri_452_104)))

## 6558
tri_6558_104 <- subset(li6800_merged, id == "6558" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6558_104)

aci_coefs[4,] <- c(id = 6558, spp = "Tri", plot = 3, subplot = 22,
                   julian_date = 104, t(coef(tri_6558_104)))

## 1666
tri_1666_104 <- subset(li6800_merged, id == "1666" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1666_104)

aci_coefs[5,] <- c(id = 1666, spp = "Tri", plot = 3, subplot = 27,
                   julian_date = 104, t(coef(tri_1666_104)))

## 4942
tri_4942_104 <- subset(li6800_merged, id == "4942" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4942_104)

aci_coefs[6,] <- c(id = 4942, spp = "Tri", plot = 3, subplot = 27,
                   julian_date = 104, t(coef(tri_4942_104)))

## 5479
tri_5479_104 <- subset(li6800_merged, id == "5479" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5479_104)

aci_coefs[7,] <- c(id = 5479, spp = "Tri", plot = 3, subplot = 21,
                   julian_date = 104, t(coef(tri_5479_104)))

## 774
tri_774_104 <- subset(li6800_merged, id == "774" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_774_104)

aci_coefs[8,] <- c(id = 774, spp = "Tri", plot = 3, subplot = 20,
                   julian_date = 104, t(coef(tri_774_104)))

## 6885
tri_6885_104 <- subset(li6800_merged, id == "6885" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6885_104)

aci_coefs[9,] <- c(id = 6885, spp = "Tri", plot = 3, subplot = 20,
                   julian_date = 104, t(coef(tri_6885_104)))

## 2329
tri_2329_104 <- subset(li6800_merged, id == "2329" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2329_104)

aci_coefs[10,] <- c(id = 2329, spp = "Tri", plot = 3, subplot = 20,
                   julian_date = 104, t(coef(tri_2329_104)))

## 4714
tri_4714_104 <- subset(li6800_merged, id == "4714" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4714_104)

aci_coefs[11,] <- c(id = 4714, spp = "Tri", plot = 3, subplot = 7,
                    julian_date = 104, t(coef(tri_4714_104)))

## 902
tri_902_104 <- subset(li6800_merged, id == "902" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_902_104)

aci_coefs[12,] <- c(id = 902, spp = "Tri", plot = 3, subplot = 7,
                    julian_date = 104, t(coef(tri_902_104)))

## 552
tri_552_104 <- subset(li6800_merged, id == "552" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_552_104)

aci_coefs[13,] <- c(id = 552, spp = "Tri", plot = 3, subplot = 14,
                    julian_date = 104, t(coef(tri_552_104)))

## 5436
tri_5436_104 <- subset(li6800_merged, id == "5436" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5436_104)

aci_coefs[14,] <- c(id = 5436, spp = "Tri", plot = 3, subplot = 3,
                    julian_date = 104, t(coef(tri_5436_104)))

## 3563
tri_3563_104 <- subset(li6800_merged, id == "3563" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3563_104)

aci_coefs[15,] <- c(id = 3563, spp = "Tri", plot = 3, subplot = 4,
                    julian_date = 104, t(coef(tri_3563_104)))

## 1926
tri_1926_104 <- subset(li6800_merged, id == "1926" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1926_104)

aci_coefs[16,] <- c(id = 1926, spp = "Tri", plot = 3, subplot = 5,
                    julian_date = 104, t(coef(tri_1926_104)))

## 425
tri_425_104 <- subset(li6800_merged, id == "425" & julian_date == 104) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_425_104)

aci_coefs[17,] <- c(id = 425, spp = "Tri", plot = 3, subplot = 6,
                    julian_date = 104, t(coef(tri_425_104)))


#####################################################################
# 4/14/24: plot 5
#####################################################################

# 2924
tri_2924_105 <- subset(li6800_merged, id == "2924" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2924_105)

aci_coefs[18,] <- c(id = 2924, spp = "Tri", plot = 5, subplot = 35,
                    julian_date = 105, t(coef(tri_2924_105)))

# 6875
tri_6875_105 <- subset(li6800_merged, id == "6875" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6875_105)

aci_coefs[19,] <- c(id = 6875, spp = "Tri", plot = 5, subplot = 29,
                    julian_date = 105, t(coef(tri_6875_105)))

# 2988
tri_2988_105 <- subset(li6800_merged, id == "2988" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2988_105)

aci_coefs[20,] <- c(id = 2988, spp = "Tri", plot = 5, subplot = 33,
                    julian_date = 105, t(coef(tri_2988_105)))

# 3829
tri_3829_105 <- subset(li6800_merged, id == "3829" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3829_105)

aci_coefs[21,] <- c(id = 3829, spp = "Tri", plot = 5, subplot = 33,
                    julian_date = 105, t(coef(tri_3829_105)))

# 5877
tri_5877_105 <- subset(li6800_merged, id == "5877" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5877_105)

aci_coefs[22,] <- c(id = 5877, spp = "Tri", plot = 5, subplot = 27,
                    julian_date = 105, t(coef(tri_5877_105)))

# 4109
tri_4109_105 <- subset(li6800_merged, id == "4109" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4109_105)

aci_coefs[23,] <- c(id = 4109, spp = "Tri", plot = 5, subplot = 32,
                    julian_date = 105, t(coef(tri_4109_105)))

# 4431
tri_4431_105 <- subset(li6800_merged, id == "4431" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4431_105)

aci_coefs[24,] <- c(id = 4431, spp = "Tri", plot = 5, subplot = 23,
                    julian_date = 105, t(coef(tri_4431_105)))

# 4959
tri_4959_105 <- subset(li6800_merged, id == "4959" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4959_105)

aci_coefs[25,] <- c(id = 4959, spp = "Tri", plot = 5, subplot = 22,
                    julian_date = 105, t(coef(tri_4959_105)))

# 5229
tri_5229_105 <- subset(li6800_merged, id == "5229" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5229_105)

aci_coefs[26,] <- c(id = 5229, spp = "Tri", plot = 5, subplot = 22,
                    julian_date = 105, t(coef(tri_5229_105)))

# 482
tri_482_105 <- subset(li6800_merged, id == "482" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_482_105)

aci_coefs[27,] <- c(id = 482, spp = "Tri", plot = 5, subplot = 22,
                    julian_date = 105, t(coef(tri_482_105)))

# 4990
tri_4990_105 <- subset(li6800_merged, id == "4990" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4990_105)

aci_coefs[28,] <- c(id = 4990, spp = "Tri", plot = 5, subplot = 22,
                    julian_date = 105, t(coef(tri_4990_105)))

# 3004
tri_3004_105 <- subset(li6800_merged, id == "3004" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3004_105)

aci_coefs[29,] <- c(id = 3004, spp = "Tri", plot = 5, subplot = 21,
                    julian_date = 105, t(coef(tri_3004_105)))

# 4576
tri_4576_105 <- subset(li6800_merged, id == "4576" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4576_105)

aci_coefs[30,] <- c(id = 4576, spp = "Tri", plot = 5, subplot = 21,
                    julian_date = 105, t(coef(tri_4576_105)))

# 3077
tri_3077_105 <- subset(li6800_merged, id == "3077" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3077_105)

aci_coefs[31,] <- c(id = 3077, spp = "Tri", plot = 5, subplot = 21,
                    julian_date = 105, t(coef(tri_3077_105)))

# 4000
tri_4000_105 <- subset(li6800_merged, id == "4000" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4000_105)

aci_coefs[32,] <- c(id = 4000, spp = "Tri", plot = 5, subplot = 20,
                    julian_date = 105, t(coef(tri_4000_105)))

# 5115
tri_5115_105 <- subset(li6800_merged, id == "5115" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5115_105)

aci_coefs[33,] <- c(id = 5115, spp = "Tri", plot = 5, subplot = 19,
                    julian_date = 105, t(coef(tri_5115_105)))

# 5228
tri_5228_105 <- subset(li6800_merged, id == "5228" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5228_105)

aci_coefs[34,] <- c(id = 5228, spp = "Tri", plot = 5, subplot = 19,
                    julian_date = 105, t(coef(tri_5228_105)))

#####################################################################
# 4/15/24: plot 6
#####################################################################
# 5229
tri_5229_106 <- subset(li6800_merged, id == "5229" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5229_106)

aci_coefs[35,] <- c(id = 5229, spp = "Tri", plot = 6, subplot = 1,
                    julian_date = 106, t(coef(tri_5229_106)))

# 762
tri_762_106 <- subset(li6800_merged, id == "762" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_762_106)

aci_coefs[36,] <- c(id = 762, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_762_106)))

# 4373
tri_4373_106 <- subset(li6800_merged, id == "4373" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4373_106)

aci_coefs[37,] <- c(id = 4373, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_4373_106)))

# 4760
tri_4760_106 <- subset(li6800_merged, id == "4760" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4760_106)

aci_coefs[38,] <- c(id = 4760, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_4760_106)))

# 4777
tri_4777_106 <- subset(li6800_merged, id == "4777" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4777_106)

aci_coefs[39,] <- c(id = 4777, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_4777_106)))

# 5865
tri_5865_106 <- subset(li6800_merged, id == "5865" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5865_106)

aci_coefs[40,] <- c(id = 5865, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_5865_106)))

# 5267
tri_5267_106 <- subset(li6800_merged, id == "5267" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5267_106)

aci_coefs[41,] <- c(id = 5267, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_5267_106)))

# 684
tri_684_106 <- subset(li6800_merged, id == "684" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_684_106)

aci_coefs[42,] <- c(id = 684, spp = "Tri", plot = 6, subplot = 2,
                    julian_date = 106, t(coef(tri_684_106)))

# 5742
tri_5742_106 <- subset(li6800_merged, id == "5742" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5742_106)

aci_coefs[43,] <- c(id = 5742, spp = "Tri", plot = 6, subplot = 3,
                    julian_date = 106, t(coef(tri_5742_106)))

# 3305
tri_3305_106 <- subset(li6800_merged, id == "3305" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3305_106)

aci_coefs[44,] <- c(id = 3305, spp = "Tri", plot = 6, subplot = 5,
                    julian_date = 106, t(coef(tri_3305_106)))

# 5619
tri_5619_106 <- subset(li6800_merged, id == "5619" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5619_106)

aci_coefs[45,] <- c(id = 5619, spp = "Tri", plot = 6, subplot = 5,
                    julian_date = 106, t(coef(tri_5619_106)))

# striped
tri_striped_106 <- subset(li6800_merged, id == "striped" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_striped_106)

aci_coefs[46,] <- c(id = "striped", spp = "Tri", plot = 6, subplot = 5,
                    julian_date = 106, t(coef(tri_striped_106)))

# 978
tri_978_106 <- subset(li6800_merged, id == "978" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_978_106)

aci_coefs[47,] <- c(id = 978, spp = "Tri", plot = 6, subplot = 6,
                    julian_date = 106, t(coef(tri_978_106)))

# 5626
tri_5626_106 <- subset(li6800_merged, id == "5626" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5626_106)

aci_coefs[48,] <- c(id = 5626, spp = "Tri", plot = 6, subplot = 6,
                    julian_date = 106, t(coef(tri_5626_106)))

# 5403
tri_5403_106 <- subset(li6800_merged, id == "5403" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5403_106)

aci_coefs[49,] <- c(id = 5403, spp = "Tri", plot = 6, subplot = 6,
                    julian_date = 106, t(coef(tri_5403_106)))

## End point (tack on additional curves here)

#####################################################################
# Snapshot measurements (Anet, gsw, Ci:Ca)
#####################################################################
snapshot <- li6800_merged %>%
  group_by(id, julian_date) %>%
  filter(row_number() == 1) %>%
  dplyr::select(id, machine, julian_date, anet = A, ci = Ci, 
                ca = Ca, co2_ref = CO2_r, gsw, Tleaf) %>%
  mutate(julian_date = as.numeric(julian_date),
         ci.ca = ci / ca,
         iwue = anet / gsw)
  
#####################################################################
# Compile A/Ci parameter estimates and snapshot measurements
#####################################################################
photo_cleaned_full <- aci_coefs %>%
  mutate(julian_date = as.numeric(julian_date),
         subplot = as.numeric(subplot),
         across(Vcmax:TPU, as.numeric)) %>%
  left_join(treatments, by = c("subplot")) %>%
  full_join(snapshot, by = c("id", "julian_date")) %>%
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
  dplyr::select(id, machine, julian_date, spp:subplot, gm.trt, Tleaf,
                anet, ci.ca, gsw, iwue, vcmax = Vcmax, vcmax25,
                jmax = Jmax, jmax25, rd = Rd, rd25) %>%
  mutate(across(Tleaf:rd25, round, 4)) %>%
  arrange(plot, julian_date)


write.csv(photo_cleaned_full,
          "../data/TT24_photo_traits_working.csv", 
          row.names = F)

