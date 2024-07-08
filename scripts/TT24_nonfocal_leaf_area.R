##############################################################################
## Prep  
##############################################################################
## Load libraries
library(tidyverse)
library(ggpubr)
library(LeafArea)

## Read .csv with leaf length and stem length for easy leaf area merge
allometry <- read.csv("../data/TT23_nonfocal_allometry.csv",
                      na.strings = c("NA", ""))


##############################################################################
## Calculate leaf disk area using ImageJ
##############################################################################

## Calculate leaf disk area
ij.path <- "/Applications/ImageJ.app"
imagepath.disk <- "../nonfocal_images/"

leaf_area <- run.ij(path.imagej = ij.path,
                    set.directory = imagepath.disk,
                    distance.pixel = 117.906,
                    known.distance = 1, low.size = 0.1,
                    set.memory = 10) %>%
  separate(sample, c("plot", "trt", "spp", "sample", "rep"))


allometry_full <- leaf_area %>%
  separate(sample, c("plot", "trt", "spp", "sample", "rep")) %>%
  mutate(rep = str_pad(rep, 2, pad = "0")) %>%
  unite("id", plot:rep, remove = FALSE) %>%
  full_join(allometry) %>%
  mutate(plot = gsub("plot", "", plot),
         focal = "n") %>%
  dplyr::select(id:spp, rep, focal, allometry_method,
                total_leaf_area = total.leaf.area,
                leaf_length, stem_length, notes)

write.csv(allometry_full, "../data/TT23_nonfocal_allometry_with_leaf_area.csv",
          row.names = FALSE)
