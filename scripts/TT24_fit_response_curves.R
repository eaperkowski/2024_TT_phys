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
         Qin_cuvette = 2000,
         
         ## Some corrections
         id = ifelse(id == "striped", "striped1", id),
         id = ifelse(julian_date == 114 & id == 1 & 
                       elapsed > 6840, 2, id),
         id = ifelse(julian_date == 119 & id == 5479 & 
                       elapsed > 4200 & elapsed < 5400, 4942, id),
         id = ifelse(julian_date == 119 & id == 5479 & 
                       elapsed > 5500, 1666, id),
         id = ifelse(julian_date == 119 & id == 392 & 
                       elapsed > 10500, 2337, id),
         id = ifelse(julian_date == 119 & id == 4934 & 
                       elapsed > 14000, 1157, id)) %>%
  dplyr::select(obs, time, elapsed, date, date_only, julian_date, hhmmss:Qin,
                Qin_cuvette, Qabs:SS_r) %>%
  arrange(machine, date, obs)

# Write merged LI-6800 file
write.csv(li6800_merged, "../data/TT24_li6800_merged.csv", 
          row.names = F)

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

# striped1
tri_striped1_106 <- subset(li6800_merged, id == "striped1" & julian_date == 106) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_striped_106)

aci_coefs[46,] <- c(id = "striped1", spp = "Tri", plot = 6, subplot = 5,
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

#####################################################################
# 4/20/24: plot 6
#####################################################################

# 5060
tri_5060_111 <- subset(li6800_merged, id == "5060" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5060_111)

aci_coefs[50,] <- c(id = 5060, spp = "Tri", plot = 6, subplot = 27,
                    julian_date = 111, t(coef(tri_5060_111)))

# 5783
tri_5783_111 <- subset(li6800_merged, id == "5783" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5783_111)

aci_coefs[51,] <- c(id = 5783, spp = "Tri", plot = 6, subplot = 20,
                    julian_date = 111, t(coef(tri_5783_111)))

# 5031
tri_5031_111 <- subset(li6800_merged, id == "5031" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5031_111)

aci_coefs[52,] <- c(id = 5031, spp = "Tri", plot = 6, subplot = 28,
                    julian_date = 111, t(coef(tri_5031_111)))

# 4543
tri_4543_111 <- subset(li6800_merged, id == "4543" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4543_111)

aci_coefs[53,] <- c(id = 4543, spp = "Tri", plot = 6, subplot = 16,
                    julian_date = 111, t(coef(tri_4543_111)))


# 2276
tri_2276_111 <- subset(li6800_merged, id == "2276" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2276_111)

aci_coefs[54,] <- c(id = 2276, spp = "Tri", plot = 6, subplot = 22,
                    julian_date = 111, t(coef(tri_2276_111)))

# 4511
tri_4511_111 <- subset(li6800_merged, id == "4511" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4511_111)

aci_coefs[55,] <- c(id = 4511, spp = "Tri", plot = 6, subplot = 16,
                    julian_date = 111, t(coef(tri_4511_111)))

# 3770
tri_3770_111 <- subset(li6800_merged, id == "3770" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3770_111)

aci_coefs[56,] <- c(id = 3770, spp = "Tri", plot = 6, subplot = 10,
                    julian_date = 111, t(coef(tri_3770_111)))

# 4582
tri_4582_111 <- subset(li6800_merged, id == "4582" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4582_111)

aci_coefs[57,] <- c(id = 4582, spp = "Tri", plot = 6, subplot = 11,
                    julian_date = 111, t(coef(tri_4582_111)))

# 2980
tri_2980_111 <- subset(li6800_merged, id == "2980" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2980_111)

aci_coefs[58,] <- c(id = 2980, spp = "Tri", plot = 6, subplot = 11,
                    julian_date = 111, t(coef(tri_2980_111)))

# 3371
tri_3371_111 <- subset(li6800_merged, id == "3371" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3371_111)

aci_coefs[59,] <- c(id = 3371, spp = "Tri", plot = 6, subplot = 11,
                    julian_date = 111, t(coef(tri_3371_111)))

# 3379
tri_3379_111 <- subset(li6800_merged, id == "3379" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3379_111)

aci_coefs[60,] <- c(id = 3379, spp = "Tri", plot = 6, subplot = 11,
                    julian_date = 111, t(coef(tri_3379_111)))

# 4505
tri_4505_111 <- subset(li6800_merged, id == "4505" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4505_111)

aci_coefs[61,] <- c(id = 4505, spp = "Tri", plot = 6, subplot = 11,
                    julian_date = 111, t(coef(tri_4505_111)))

# 4745
tri_4745_111 <- subset(li6800_merged, id == "4745" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4745_111)

aci_coefs[62,] <- c(id = 4745, spp = "Tri", plot = 6, subplot = 12,
                    julian_date = 111, t(coef(tri_4745_111)))

# 5722
tri_5722_111 <- subset(li6800_merged, id == "5722" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5722_111)

aci_coefs[63,] <- c(id = 5722, spp = "Tri", plot = 6, subplot = 6,
                    julian_date = 111, t(coef(tri_5722_111)))

# 2921
tri_2921_111 <- subset(li6800_merged, id == "2921" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2921_111)

aci_coefs[64,] <- c(id = 2921, spp = "Tri", plot = 6, subplot = 6,
                    julian_date = 111, t(coef(tri_2921_111)))

# striped2
tri_striped2_111 <- subset(li6800_merged, id == "striped2" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_striped2_111)

aci_coefs[65,] <- c(id = "striped2", spp = "Tri", plot = 6, subplot = 5,
                    julian_date = 111, t(coef(tri_striped2_111)))

# 5642
tri_5642_111 <- subset(li6800_merged, id == "5642" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5642_111)

aci_coefs[66,] <- c(id = 5642, spp = "Tri", plot = 6, subplot = 4,
                    julian_date = 111, t(coef(tri_5642_111)))

# 6879
tri_6879_111 <- subset(li6800_merged, id == "6879" & julian_date == 111) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6879_111)

aci_coefs[67,] <- c(id = 6879, spp = "Tri", plot = 6, subplot = 36,
                    julian_date = 111, t(coef(tri_6879_111)))

#####################################################################
# 4/22/24: plot 5
#####################################################################

# 7120
tri_7120_113 <- subset(li6800_merged, id == "7120" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_7120_113)

aci_coefs[68,] <- c(id = 7120, spp = "Tri", plot = 5, subplot = 17,
                    julian_date = 113, t(coef(tri_7120_113)))

# 2912
tri_2912_113 <- subset(li6800_merged, id == "2912" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2912_113)

aci_coefs[69,] <- c(id = 2912, spp = "Tri", plot = 5, subplot = 11,
                    julian_date = 113, t(coef(tri_2912_113)))

# 614
tri_614_113 <- subset(li6800_merged, id == "614" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_614_113)

aci_coefs[70,] <- c(id = 614, spp = "Tri", plot = 5, subplot = 11,
                    julian_date = 113, t(coef(tri_614_113)))

# 1374
tri_1374_113 <- subset(li6800_merged, id == "1374" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1374_113)

aci_coefs[71,] <- c(id = 1374, spp = "Tri", plot = 5, subplot = 10,
                    julian_date = 113, t(coef(tri_1374_113)))

# 6881
tri_6881_113 <- subset(li6800_merged, id == "6881" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6881_113)

aci_coefs[72,] <- c(id = 6881, spp = "Tri", plot = 5, subplot = 10,
                    julian_date = 113, t(coef(tri_6881_113)))

# 5381
tri_5381_113 <- subset(li6800_merged, id == "5381" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5381_113)

aci_coefs[73,] <- c(id = 5381, spp = "Tri", plot = 5, subplot = 10,
                    julian_date = 113, t(coef(tri_5381_113)))

# 4105
tri_4105_113 <- subset(li6800_merged, id == "4105" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4105_113)

aci_coefs[74,] <- c(id = 4105, spp = "Tri", plot = 5, subplot = 10,
                    julian_date = 113, t(coef(tri_4105_113)))

# 4547
tri_4547_113 <- subset(li6800_merged, id == "4547" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4547_113)

aci_coefs[75,] <- c(id = 4547, spp = "Tri", plot = 5, subplot = 10,
                    julian_date = 113, t(coef(tri_4547_113)))

# 2388
tri_2388_113 <- subset(li6800_merged, id == "2388" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2388_113)

aci_coefs[76,] <- c(id = 2388, spp = "Tri", plot = 5, subplot = 10,
                    julian_date = 113, t(coef(tri_2388_113)))

# 5904
tri_5904_113 <- subset(li6800_merged, id == "5904" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5904_113)

aci_coefs[77,] <- c(id = 5904, spp = "Tri", plot = 5, subplot = 9,
                    julian_date = 113, t(coef(tri_5904_113)))

# 7147
tri_7147_113 <- subset(li6800_merged, id == "7147" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_7147_113)

aci_coefs[78,] <- c(id = 7147, spp = "Tri", plot = 5, subplot = 9,
                    julian_date = 113, t(coef(tri_7147_113)))

# 43
tri_43_113 <- subset(li6800_merged, id == "43" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_43_113)

aci_coefs[79,] <- c(id = 43, spp = "Tri", plot = 5, subplot = 9,
                    julian_date = 113, t(coef(tri_43_113)))

# 86
tri_86_113 <- subset(li6800_merged, id == "86" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_86_113)

aci_coefs[80,] <- c(id = 86, spp = "Tri", plot = 5, subplot = 9,
                    julian_date = 113, t(coef(tri_86_113)))

# flag8
tri_flag8_113 <- subset(li6800_merged, id == "flag8" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_43_113)

aci_coefs[81,] <- c(id = "flag8", spp = "Tri", plot = 5, subplot = 16,
                    julian_date = 113, t(coef(tri_flag8_113)))

# 2563
tri_2563_113 <- subset(li6800_merged, id == "2563" & julian_date == 113) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2563_113)

aci_coefs[82,] <- c(id = 2563, spp = "Tri", plot = 5, subplot = 9,
                    julian_date = 113, t(coef(tri_2563_113)))

#####################################################################
# 4/23/24: plot 5
#####################################################################

# 4265
tri_4265_114 <- subset(li6800_merged, id == "4265" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4265_114)

aci_coefs[83,] <- c(id = 4265, spp = "Tri", plot = 5, subplot = 31,
                    julian_date = 114, t(coef(tri_4265_114)))

# 2573
tri_2573_114 <- subset(li6800_merged, id == "2573" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2573_114)

aci_coefs[84,] <- c(id = 2573, spp = "Tri", plot = 5, subplot = 31,
                    julian_date = 114, t(coef(tri_2573_114)))

# 2547
tri_2547_114 <- subset(li6800_merged, id == "2547" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2547_114)

aci_coefs[85,] <- c(id = 2547, spp = "Tri", plot = 5, subplot = 19,
                    julian_date = 114, t(coef(tri_2547_114)))

# 4177
tri_4177_114 <- subset(li6800_merged, id == "4177" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4177_114)

aci_coefs[86,] <- c(id = 4177, spp = "Tri", plot = 5, subplot = 14,
                    julian_date = 114, t(coef(tri_4177_114)))

# 1795
tri_1795_114 <- subset(li6800_merged, id == "1795" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1795_114)

aci_coefs[87,] <- c(id = 1795, spp = "Tri", plot = 5, subplot = 14,
                    julian_date = 114, t(coef(tri_1795_114)))


aci_coefs <- read.csv("../data/TT24_photo_traits_working.csv") %>%
  select(id, spp, plot, subplot, julian_date, Vcmax = vcmax, 
         Jmax = jmax, Rd = rd) %>%
  mutate(TPU = NA)


# 3
tri_flag3_114 <- subset(li6800_merged, id == "3" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag3_114)

aci_coefs[88,] <- c(id = "flag3_tri", spp = "Tri", plot = 5, subplot = 14,
                    julian_date = 114, t(coef(tri_flag3_114)))

# 1
tri_flag1_114 <- subset(li6800_merged, id == "1" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag1_114)

aci_coefs[89,] <- c(id = "flag1_tri", spp = "Tri", plot = 5, subplot = 14,
                    julian_date = 114, t(coef(tri_flag1_114)))

# 2
tri_flag2_114 <- subset(li6800_merged, id == "2" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_flag2_114)

aci_coefs[90,] <- c(id = "flag2_tri", spp = "Tri", plot = 5, subplot = 14,
                    julian_date = 114, t(coef(tri_flag2_114)))

# 1739
tri_1739_114 <- subset(li6800_merged, id == "1739" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1739_114)

aci_coefs[91,] <- c(id = 1739, spp = "Tri", plot = 5, subplot = 14,
                    julian_date = 114, t(coef(tri_1739_114)))

# 4149
tri_4149_114 <- subset(li6800_merged, id == "4149" & julian_date == 114) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4149_114)

aci_coefs[92,] <- c(id = 4149, spp = "Tri", plot = 5, subplot = 13,
                    julian_date = 114, t(coef(tri_4149_114)))


#####################################################################
# 4/28/24: plot 3
#####################################################################

# 2329
tri_2329_119 <- subset(li6800_merged, id == "2329" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2329_119)

aci_coefs[93,] <- c(id = 2329, spp = "Tri", plot = 3, subplot = 20,
                    julian_date = 119, t(coef(tri_2329_119)))

# 6885
tri_6885_119 <- subset(li6800_merged, id == "6885" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6885_119)

aci_coefs[94,] <- c(id = 6885, spp = "Tri", plot = 3, subplot = 20,
                    julian_date = 119, t(coef(tri_6885_119)))

# 774
tri_774_119 <- subset(li6800_merged, id == "774" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_774_119)

aci_coefs[95,] <- c(id = 774, spp = "Tri", plot = 3, subplot = 20,
                    julian_date = 119, t(coef(tri_774_119)))

# 5479
tri_5479_119 <- subset(li6800_merged, id == "5479" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5479_119)

aci_coefs[96,] <- c(id = 5479, spp = "Tri", plot = 3, subplot = 21,
                    julian_date = 119, t(coef(tri_5479_119)))

# 4942
tri_4942_119 <- subset(li6800_merged, id == "4942" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4942_119)

aci_coefs[97,] <- c(id = 4942, spp = "Tri", plot = 3, subplot = 27,
                    julian_date = 119, t(coef(tri_4942_119)))

# 1666
tri_1666_119 <- subset(li6800_merged, id == "1666" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1666_119)

aci_coefs[98,] <- c(id = 1666, spp = "Tri", plot = 3, subplot = 27,
                    julian_date = 119, t(coef(tri_1666_119)))

# 6558
tri_6558_119 <- subset(li6800_merged, id == "6558" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6558_119)

aci_coefs[99,] <- c(id = 6558, spp = "Tri", plot = 3, subplot = 22,
                    julian_date = 119, t(coef(tri_6558_119)))

# 392
tri_392_119 <- subset(li6800_merged, id == "392" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_392_119)

aci_coefs[100,] <- c(id = 392, spp = "Tri", plot = 3, subplot = 22,
                    julian_date = 119, t(coef(tri_392_119)))

# 2337
mai_2337_119 <- subset(li6800_merged, id == "2337" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2337_119)

aci_coefs[101,] <- c(id = 2337, spp = "Mai", plot = 3, subplot = 22,
                     julian_date = 119, t(coef(mai_2337_119)))

# 583
tri_583_119 <- subset(li6800_merged, id == "583" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_583_119)

aci_coefs[102,] <- c(id = 583, spp = "Tri", plot = 3, subplot = 23,
                     julian_date = 119, t(coef(tri_583_119)))

# 4934
tri_4934_119 <- subset(li6800_merged, id == "4934" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4934_119)

aci_coefs[103,] <- c(id = 4934, spp = "Tri", plot = 3, subplot = 23,
                     julian_date = 119, t(coef(tri_4934_119)))

# 1157
mai_1157_119 <- subset(li6800_merged, id == "1157" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1157_119)

aci_coefs[104,] <- c(id = 1157, spp = "Mai", plot = 3, subplot = 30,
                     julian_date = 119, t(coef(mai_1157_119)))

# 179
mai_179_119 <- subset(li6800_merged, id == "174" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_179_119)

aci_coefs[105,] <- c(id = 179, spp = "Mai", plot = 3, subplot = 30,
                     julian_date = 119, t(coef(mai_179_119)))

# 4714
tri_4714_119 <- subset(li6800_merged, id == "4714" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_4714_119)

aci_coefs[106,] <- c(id = 4714, spp = "Tri", plot = 3, subplot = 7,
                     julian_date = 119, t(coef(tri_4714_119)))

# 902
tri_902_119 <- subset(li6800_merged, id == "902" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_902_119)

aci_coefs[107,] <- c(id = 902, spp = "Tri", plot = 3, subplot = 7,
                     julian_date = 119, t(coef(tri_902_119)))

# 5105
mai_5105_119 <- subset(li6800_merged, id == "5105" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5105_119)

aci_coefs[108,] <- c(id = 5105, spp = "Mai", plot = 3, subplot = 8,
                     julian_date = 119, t(coef(mai_5105_119)))

# 5184
mai_5184_119 <- subset(li6800_merged, id == "5184" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5184_119)

aci_coefs[109,] <- c(id = 5184, spp = "Mai", plot = 3, subplot = 13,
                     julian_date = 119, t(coef(mai_5184_119)))

# 141
mai_141_119 <- subset(li6800_merged, id == "141" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_141_119)

aci_coefs[110,] <- c(id = 141, spp = "Mai", plot = 3, subplot = 13,
                     julian_date = 119, t(coef(mai_141_119)))

# 552
tri_552_119 <- subset(li6800_merged, id == "552" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_552_119)

aci_coefs[111,] <- c(id = 552, spp = "Tri", plot = 3, subplot = 8,
                     julian_date = 119, t(coef(tri_552_119)))

# 4250
mai_4250_119 <- subset(li6800_merged, id == "4250" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4250_119)

aci_coefs[112,] <- c(id = 4250, spp = "Mai", plot = 3, subplot = 2,
                     julian_date = 119, t(coef(mai_4250_119)))

# 5488
tri_5488_119 <- subset(li6800_merged, id == "5488" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5488_119)

aci_coefs[113,] <- c(id = 5488, spp = "Tri", plot = 3, subplot = 1,
                     julian_date = 119, t(coef(tri_5488_119)))

# 9412
mai_9412_119 <- subset(li6800_merged, id == "9412" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_9412_119)

aci_coefs[114,] <- c(id = 9412, spp = "Mai", plot = 3, subplot = 3,
                     julian_date = 119, t(coef(mai_9412_119)))

# 5436
tri_5436_119 <- subset(li6800_merged, id == "5436" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5436_119)

aci_coefs[115,] <- c(id = 5436, spp = "Tri", plot = 3, subplot = 3,
                     julian_date = 119, t(coef(tri_5436_119)))

# 5495
tri_5495_119 <- subset(li6800_merged, id == "5495" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5495_119)

aci_coefs[116,] <- c(id = 5495, spp = "Tri", plot = 3, subplot = 3,
                     julian_date = 119, t(coef(tri_5495_119)))

# 3563
tri_3563_119 <- subset(li6800_merged, id == "3563" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_3563_119)

aci_coefs[117,] <- c(id = 3563, spp = "Tri", plot = 3, subplot = 4,
                     julian_date = 119, t(coef(tri_3563_119)))

# 2310
mai_2310_119 <- subset(li6800_merged, id == "2310" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2310_119)

aci_coefs[118,] <- c(id = 2310, spp = "Mai", plot = 3, subplot = 4,
                     julian_date = 119, t(coef(mai_2310_119)))

# 1926
tri_1926_119 <- subset(li6800_merged, id == "1926" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1926_119)

aci_coefs[119,] <- c(id = 1926, spp = "Tri", plot = 3, subplot = 5,
                     julian_date = 119, t(coef(tri_1926_119)))

# 5500
tri_5500_119 <- subset(li6800_merged, id == "5500" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5500_119)

aci_coefs[120,] <- c(id = 5500, spp = "Tri", plot = 3, subplot = 5,
                     julian_date = 119, t(coef(tri_5500_119)))

# 1686
tri_1686_119 <- subset(li6800_merged, id == "1686" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_1686_119)

aci_coefs[121,] <- c(id = 1686, spp = "Tri", plot = 3, subplot = 6,
                     julian_date = 119, t(coef(tri_1686_119)))

# 5797
tri_5797_119 <- subset(li6800_merged, id == "5797" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_5797_119)

aci_coefs[122,] <- c(id = 5797, spp = "Tri", plot = 3, subplot = 6,
                     julian_date = 119, t(coef(tri_5797_119)))

# 425
tri_425_119 <- subset(li6800_merged, id == "425" & julian_date == 119) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_425_119)

aci_coefs[123,] <- c(id = 425, spp = "Tri", plot = 3, subplot = 6,
                     julian_date = 119, t(coef(tri_425_119)))

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
         iwue = anet / gsw,
         id = ifelse(id == 1 & julian_date == 114, "flag1_tri", id),
         id = ifelse(id == 2 & julian_date == 114, "flag2_tri", id),
         id = ifelse(id == 3 & julian_date == 114, "flag3_tri", id))
  
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
  mutate(across(Tleaf:rd25, \(x) round(x, digits = 4))) %>%
  arrange(plot, julian_date)


write.csv(photo_cleaned_full,
          "../data/TT24_photo_traits_working.csv", 
          row.names = F)
