# TT24_fit_response_curves.R
# R script that merges LI-6800 files from "../cleaned_li6800/" and
# uses merged file to fit CO2 response curves
# 
# Note: all paths assume that the folder containing this R script
# is the working root directory

#####################################################################
# Libraries
#####################################################################
library(tidyverse)
library(plantecophys)

#####################################################################
# Merge cleaned LI-6800 files
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
#          row.names = F)

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
#plot(tri_583_104)

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
                    julian_date = 104, t(coef(tri_2924_105)))


# 6875
tri_6875_105 <- subset(li6800_merged, id == "6875" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_6875_105)

aci_coefs[19,] <- c(id = 6875, spp = "Tri", plot = 5, subplot = 29,
                    julian_date = 104, t(coef(tri_6875_105)))

# 2988
tri_2988_105 <- subset(li6800_merged, id == "2988" & julian_date == 105) %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(tri_2988_105)

aci_coefs[20,] <- c(id = 2988, spp = "Tri", plot = 5, subplot = 33,
                    julian_date = 104, t(coef(tri_2988_105)))

## STOPPING POINT 4/17/24
write.csv(aci_coefs, "../data/TT24_aci_coefs_working.csv", row.names = F)

