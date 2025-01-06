# TT24_create_C_budget_params.R
#
# R script that determines the effect of gm.trt across growing 
# season. Note: script uses two approaches: (1) by accounting 
# for size classes by discretizing individuals evenly across 
# size distributions, (2) by ignoring size distribution entirely
# 
# Note: all paths assume that the folder containing this R script
# is the working root directory

#####################################################################
# Libraries and file read in
#####################################################################
# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)
library(ggpubr)

# Read .csv file
photo_traits <- read.csv("../data/TT24_photo_traits.csv") %>% 
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")))

# Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

# Color palettes
gm.colors <- c("#00B2BE", "#F1B700")

#####################################################################
# Discretize size distribution into 5 size classes - Trillium
#####################################################################
# Visualize Trillium size distribution
ggplot(data = subset(distinct(photo_traits, id, .keep_all = T), spp == "Tri"),
       aes(x = init_leaf_area)) +
  geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))), 
                 bins = 5, fill = "gray", colour = "black") +
  labs(x = expression("Initial total leaf area (cm"^"2"*")"),
       y = "Relative Frequency") +
  theme_classic(base_size = 18)

# Create data frame that includes two options for bin ranges
tri_intervals <- data.frame(
  id = subset(
    distinct(photo_traits, id, .keep_all = T), spp == "Tri")$id,
  
  spp = "Tri",
  
  init_leaf_area = subset(
    distinct(photo_traits, id, .keep_all = T), spp == "Tri")$init_leaf_area,
  
  # bin ranges if the range of total leaf area is equally 
  # divided into 5 bins
  range_interval = cut_interval(
    x = subset(distinct(photo_traits, id, .keep_all = T), 
               spp == "Tri")$init_leaf_area, n = 5),
  
  # bin ranges if the number of individuals are equally 
  # divided into 5 bins
  number_interval = cut_number(
    x = subset(distinct(photo_traits, id, .keep_all = T), 
               spp == "Tri")$init_leaf_area, n = 5))

# Bin ranges if total leaf area is equally divided into 5 bins
unique(tri_intervals$range_interval)
# Levels: [53.4-167], [167-280], [280-394], [394-507], [507-621]

# Bin ranges such that each bin contains the same number of inds
unique(tri_intervals$number_interval)
# Levels: [53.4-99.7], [99.7-142], [142-203], [203-322], [322-621]

#####################################################################
# Discretize size distribution into 5 size classes - Maianthemum
#####################################################################
# Visualize Trillium size distribution
ggplot(data = subset(distinct(photo_traits, id, .keep_all = T), spp == "Mai"),
       aes(x = init_leaf_area)) +
  geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))), 
                 bins = 5, fill = "gray", colour = "black") +
  labs(x = expression("Initial total leaf area (cm"^"2"*")"),
       y = "Relative Frequency") +
  theme_classic(base_size = 18)

# Create data frame that includes two options for bin ranges
mai_intervals <- data.frame(
  id = subset(
    distinct(photo_traits, id, .keep_all = T), spp == "Mai")$id,
  
  spp = "Mai",
  
  init_leaf_area = subset(
    distinct(photo_traits, id, .keep_all = T), spp == "Mai")$init_leaf_area,
  
  # bin ranges if the range of total leaf area is equally 
  # divided into 5 bins
  range_interval = cut_interval(
    x = subset(distinct(photo_traits, id, .keep_all = T), 
               spp == "Mai")$init_leaf_area, n = 5),
  
  # bin ranges if the number of individuals are equally 
  # divided into 5 bins
  number_interval = cut_number(
    x = subset(distinct(photo_traits, id, .keep_all = T), 
               spp == "Mai")$init_leaf_area, n = 5))

# Bin ranges if total leaf area is equally divided into 5 bins
unique(mai_intervals$range_interval)
# Levels: [59.8-271], [271-483], [483-695], [695-906], [906-1120]

# Bin ranges such that each bin contains the same number of inds
unique(mai_intervals$number_interval)
# Levels: [59.8-201], [201-349], [349-554], [554-673], [673-1120]

#####################################################################
# Merge Tri and Mai intervals and join with photo_traits
#####################################################################

# Merge Tri and Mai intervals
intervals <- tri_intervals %>%
  full_join(mai_intervals)

# Join with photo_traits
photo_traits <- photo_traits %>%
  full_join(intervals)

# Take a peek at dataset
head(photo_traits)

# What is the min doy measurement for Tri and Mai?
min(subset(photo_traits, spp == "Tri")$doy)
min(subset(photo_traits, spp == "Mai")$doy)

# What is the min doy measurement for Tri and Mai?
max(subset(photo_traits, spp == "Tri")$doy)
max(subset(photo_traits, spp == "Mai")$doy)

#####################################################################
# Vcmax - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
vcmax_tri_size <- lmer(sqrt(vcmax) ~ gm.trt * doy * number_interval + (1|id), 
                  data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(vcmax_tri_size)
qqnorm(residuals(vcmax_tri_size))
qqline(residuals(vcmax_tri_size))
densityPlot(residuals(vcmax_tri_size))
shapiro.test(residuals(vcmax_tri_size))
outlierTest(vcmax_tri_size)

# Model output
summary(vcmax_tri_size)
Anova(vcmax_tri_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_vcmax_predicted <- data.frame(
  spp = "Tri",
  emmeans(vcmax_tri_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax_predicted = response,
                vcmax_SE = SE,
                vcmax_LCL = lower.CL,
                vcmax_UCL = upper.CL)
  
#####################################################################
# Vcmax25 - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
vcmax25_tri_size <- lmer(sqrt(vcmax25) ~ gm.trt * doy * number_interval + (1|id), 
                       data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(vcmax25_tri_size)
qqnorm(residuals(vcmax25_tri_size))
qqline(residuals(vcmax25_tri_size))
densityPlot(residuals(vcmax25_tri_size))
shapiro.test(residuals(vcmax25_tri_size))
outlierTest(vcmax25_tri_size)

# Model output
summary(vcmax25_tri_size)
Anova(vcmax25_tri_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_vcmax25_predicted <- data.frame(
  spp = "Tri",
  emmeans(vcmax25_tri_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax25_predicted = response,
                vcmax25_SE = SE,
                vcmax25_LCL = lower.CL,
                vcmax25_UCL = upper.CL)

#####################################################################
# Jmax - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
photo_traits$jmax[c(62, 296)] <- NA

jmax_tri_size <- lmer(sqrt(jmax) ~ gm.trt * doy * number_interval + (1|id), 
                       data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(jmax_tri_size)
qqnorm(residuals(jmax_tri_size))
qqline(residuals(jmax_tri_size))
densityPlot(residuals(jmax_tri_size))
shapiro.test(residuals(jmax_tri_size))
outlierTest(jmax_tri_size)

# Model output
summary(jmax_tri_size)
Anova(jmax_tri_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_jmax_predicted <- data.frame(
  spp = "Tri",
  emmeans(jmax_tri_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax_predicted = response,
                jmax_SE = SE,
                jmax_LCL = lower.CL,
                jmax_UCL = upper.CL)

#####################################################################
# Jmax25 - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
photo_traits$jmax25[c(62, 296)] <- NA

jmax25_tri_size <- lmer(sqrt(jmax25) ~ gm.trt * doy * number_interval + (1|id), 
                         data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(jmax25_tri_size)
qqnorm(residuals(jmax25_tri_size))
qqline(residuals(jmax25_tri_size))
densityPlot(residuals(jmax25_tri_size))
shapiro.test(residuals(jmax25_tri_size))
outlierTest(jmax25_tri_size)

# Model output
summary(jmax25_tri_size)
Anova(jmax25_tri_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_jmax25_predicted <- data.frame(
  spp = "Tri",
  emmeans(jmax25_tri_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax25_predicted = response,
                jmax25_SE = SE,
                jmax25_LCL = lower.CL,
                jmax25_UCL = upper.CL)

#####################################################################
# Ci - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
ci_tri_size <- lmer(ci ~ gm.trt * doy * number_interval + (1|id), 
                    data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(ci_tri_size)
qqnorm(residuals(ci_tri_size))
qqline(residuals(ci_tri_size))
densityPlot(residuals(ci_tri_size))
shapiro.test(residuals(ci_tri_size))
outlierTest(ci_tri_size)

# Model output
summary(ci_tri_size)
Anova(ci_tri_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_ci_predicted <- data.frame(
  spp = "Tri",
  emmeans(ci_tri_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(100, 196, 1)))) %>%
  dplyr::select(gm.trt:doy, spp,
                ci_predicted = emmean,
                ci_SE = SE,
                ci_LCL = lower.CL,
                ci_UCL = upper.CL)

#####################################################################
# Write data frame that compiles Trillium Ci, Vcmax, Jmax predictions
# across range in DOY accounting for gm.trt and size class intervals
#####################################################################
tri_predictions <- tri_ci_predicted %>%
  full_join(tri_vcmax_predicted) %>%
  full_join(tri_vcmax25_predicted) %>%
  full_join(tri_jmax_predicted) %>%
  full_join(tri_jmax25_predicted) %>%
  mutate(size_trt = factor(str_c(number_interval, "_", gm.trt),
                           levels = c("[53.4,99.7]_ambient",
                                      "[53.4,99.7]_weeded",
                                      "(99.7,142]_ambient",
                                      "(99.7,142]_weeded",
                                      "(142,203]_ambient",
                                      "(142,203]_weeded",
                                      "(203,322]_ambient",
                                      "(203,322]_weeded",
                                      "(322,621]_ambient",
                                      "(322,621]_weeded")))

# Visualize patterns
ggplot(data = tri_predictions, aes(x = doy, y = vcmax25_predicted)) +
  geom_smooth(method = "loess", aes(group = size_trt)) +
  geom_ribbon(aes(ymin = vcmax25_LCL,
                  ymax = vcmax25_UCL, 
                  group = size_trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Vcmax") +
  facet_wrap(~size_trt) +
  theme_bw(base_size = 18)

ggplot(data = tri_predictions, aes(x = doy, y = jmax25_predicted)) +
  geom_smooth(method = "loess", aes(group = size_trt)) +
  geom_ribbon(aes(ymin = jmax25_LCL,
                  ymax = jmax25_UCL, 
                  group = size_trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Jmax25") +
  facet_wrap(~size_trt) +
  theme_bw(base_size = 18)

ggplot(data = tri_predictions, aes(x = doy, y = ci_predicted)) +
  geom_smooth(method = "loess", aes(group = size_trt)) +
  geom_ribbon(aes(ymin = ci_LCL,
                  ymax = ci_UCL, 
                  group = size_trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(200, 400), breaks = seq(200, 400, 100)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Ci") +
  facet_wrap(~size_trt) +
  theme_bw(base_size = 18)

#####################################################################
# Vcmax - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
vcmax_mai_size <- lmer(log(vcmax) ~ gm.trt * doy * number_interval + (1|id), 
                       data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(vcmax_mai_size)
qqnorm(residuals(vcmax_mai_size))
qqline(residuals(vcmax_mai_size))
densityPlot(residuals(vcmax_mai_size))
shapiro.test(residuals(vcmax_mai_size))
outlierTest(vcmax_mai_size)

# Model output
summary(vcmax_mai_size)
Anova(vcmax_mai_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_vcmax_predicted <- data.frame(
  spp = "Mai",
  emmeans(vcmax_mai_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax_predicted = response,
                vcmax_SE = SE,
                vcmax_LCL = lower.CL,
                vcmax_UCL = upper.CL)

#####################################################################
# Vcmax25 - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
vcmax25_mai_size <- lmer(log(vcmax25) ~ gm.trt * doy * number_interval + (1|id), 
                         data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(vcmax25_mai_size)
qqnorm(residuals(vcmax25_mai_size))
qqline(residuals(vcmax25_mai_size))
densityPlot(residuals(vcmax25_mai_size))
shapiro.test(residuals(vcmax25_mai_size))
outlierTest(vcmax25_mai_size)

# Model output
summary(vcmax25_mai_size)
Anova(vcmax25_mai_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_vcmax25_predicted <- data.frame(
  spp = "Mai",
  emmeans(vcmax25_mai_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax25_predicted = response,
                vcmax25_SE = SE,
                vcmax25_LCL = lower.CL,
                vcmax25_UCL = upper.CL)


#####################################################################
# Jmax - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
jmax_mai_size <- lmer(log(jmax) ~ gm.trt * doy * number_interval + (1|id), 
                      data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(jmax_mai_size)
qqnorm(residuals(jmax_mai_size))
qqline(residuals(jmax_mai_size))
densityPlot(residuals(jmax_mai_size))
shapiro.test(residuals(jmax_mai_size))
outlierTest(jmax_mai_size)

# Model output
summary(jmax_mai_size)
Anova(jmax_mai_size)

# Write data.frame that predicts Jmax for a given day, plot, and size class
mai_jmax_predicted <- data.frame(
  spp = "Mai",
  emmeans(jmax_mai_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
    dplyr::select(gm.trt:doy, spp,
                  jmax_predicted = response,
                  jmax_SE = SE,
                  jmax_LCL = lower.CL,
                  jmax_UCL = upper.CL)

#####################################################################
# Jmax25 - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
jmax25_mai_size <- lmer(log(jmax25) ~ gm.trt * doy * number_interval + (1|id), 
                        data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(jmax25_mai_size)
qqnorm(residuals(jmax25_mai_size))
qqline(residuals(jmax25_mai_size))
densityPlot(residuals(jmax25_mai_size))
shapiro.test(residuals(jmax25_mai_size))
outlierTest(jmax25_mai_size)

# Model output
summary(jmax25_mai_size)
Anova(jmax25_mai_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_jmax25_predicted <- data.frame(
  spp = "Mai",
  emmeans(jmax25_mai_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax25_predicted = response,
                jmax25_SE = SE,
                jmax25_LCL = lower.CL,
                jmax25_UCL = upper.CL)

#####################################################################
# Ci - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
photo_traits$ci[c(367, 392, 450, 790, 791, 820)] <- NA

ci_mai_size <- lmer(ci ~ gm.trt * doy * number_interval + (1|id), 
                    data = subset(photo_traits, spp == "Mai" & ci > 0 & ci < 450))

# Check model assumptions
plot(ci_mai_size)
qqnorm(residuals(ci_mai_size))
qqline(residuals(ci_mai_size))
densityPlot(residuals(ci_mai_size))
shapiro.test(residuals(ci_mai_size))
outlierTest(ci_mai_size)

# Model output
summary(ci_tri_size)
Anova(ci_tri_size)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_ci_predicted <- data.frame(
  spp = "Mai",
  emmeans(ci_mai_size,
          ~gm.trt*number_interval,
          "doy",
          at = list(doy = seq(115, 250, 1)))) %>%
  dplyr::select(gm.trt:doy, spp,
                ci_predicted = emmean,
                ci_SE = SE,
                ci_LCL = lower.CL,
                ci_UCL = upper.CL)

#####################################################################
# Write data frame that compiles Maianthemum Ci, Vcmax, Jmax 
# predictions across range in DOY accounting for gm.trt and size 
# class intervals
#####################################################################

mai_predictions <- mai_ci_predicted %>%
  full_join(mai_vcmax_predicted) %>%
  full_join(mai_vcmax25_predicted) %>%
  full_join(mai_jmax_predicted) %>%
  full_join(mai_jmax25_predicted) %>%
  mutate(size_trt = factor(str_c(number_interval, "_", gm.trt),
                           levels = c("[59.8,201]_ambient",
                                      "[59.8,201]_weeded",
                                      "(201,349]_ambient",
                                      "(201,349]_weeded",
                                      "(349,554]_ambient",
                                      "(349,554]_weeded",
                                      "(554,673]_ambient",
                                      "(554,673]_weeded",
                                      "(673,1.12e+03]_ambient",
                                      "(673,1.12e+03]_weeded")))

# Visualize patterns
ggplot(data = mai_predictions, aes(x = doy, y = vcmax25_predicted)) +
  geom_smooth(method = "loess", aes(group = size_trt)) +
  geom_ribbon(aes(ymin = vcmax25_LCL,
                  ymax = vcmax25_UCL, 
                  group = size_trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Vcmax25") +
  facet_wrap(~size_trt) +
  theme_bw(base_size = 18)

ggplot(data = mai_predictions, aes(x = doy, y = jmax25_predicted)) +
  geom_smooth(method = "loess", aes(group = size_trt)) +
  geom_ribbon(aes(ymin = jmax25_LCL,
                  ymax = jmax25_UCL, 
                  group = size_trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Jmax25") +
  facet_wrap(~size_trt) +
  theme_bw(base_size = 18)

ggplot(data = mai_predictions, aes(x = doy, y = ci_predicted)) +
  geom_smooth(method = "loess", aes(group = size_trt)) +
  geom_ribbon(aes(ymin = ci_LCL,
                  ymax = ci_UCL, 
                  group = size_trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(200, 400), breaks = seq(200, 400, 100)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Ci") +
  facet_wrap(~size_trt) +
  theme_bw(base_size = 18)

#####################################################################
# Combine model predictions for Trillium and Maianthemum responses to
# gm treatment across growing season and size classes
#####################################################################

# Merge data frames and rearrange columns
combined_predictions <- tri_predictions %>%
  full_join(mai_predictions) %>%
  dplyr::select(gm_trt = gm.trt, number_interval, size_trt,
                doy:jmax25_UCL)

# Write combined_predictions to .csv
write.csv(combined_predictions, 
          "../data/c_budget/TT24_budget_params_with_size.csv",
          row.names = F)

#####################################################################
#####################################################################
# Now, create model predictions but ignore size classes
#####################################################################
#####################################################################

#####################################################################
# Vcmax - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
photo_traits$vcmax[452] <- NA

vcmax_tri_nosize <- lmer(sqrt(vcmax) ~ gm.trt * doy + (1|id), 
                       data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(vcmax_tri_nosize)
qqnorm(residuals(vcmax_tri_nosize))
qqline(residuals(vcmax_tri_nosize))
densityPlot(residuals(vcmax_tri_nosize))
shapiro.test(residuals(vcmax_tri_nosize))
outlierTest(vcmax_tri_nosize)

# Model output
summary(vcmax_tri_nosize)
Anova(vcmax_tri_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_vcmax_predicted_nosize <- data.frame(
  spp = "Tri",
  emmeans(vcmax_tri_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax_predicted = response,
                vcmax_SE = SE,
                vcmax_LCL = lower.CL,
                vcmax_UCL = upper.CL)

#####################################################################
# Vcmax25 - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
photo_traits$vcmax25[c(281, 296, 452)] <- NA

vcmax25_tri_nosize <- lmer(log(vcmax25) ~ gm.trt * doy + (1|id), 
                         data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(vcmax25_tri_nosize)
qqnorm(residuals(vcmax25_tri_nosize))
qqline(residuals(vcmax25_tri_nosize))
densityPlot(residuals(vcmax25_tri_nosize))
shapiro.test(residuals(vcmax25_tri_nosize))
outlierTest(vcmax25_tri_nosize)

# Model output
summary(vcmax25_tri_nosize)
Anova(vcmax25_tri_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_vcmax25_predicted_nosize <- data.frame(
  spp = "Tri",
  emmeans(vcmax25_tri_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax25_predicted = response,
                vcmax25_SE = SE,
                vcmax25_LCL = lower.CL,
                vcmax25_UCL = upper.CL)

#####################################################################
# Jmax - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
jmax_tri_nosize <- lmer(jmax ~ gm.trt * doy + (1|id), 
                      data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(jmax_tri_nosize)
qqnorm(residuals(jmax_tri_nosize))
qqline(residuals(jmax_tri_nosize))
densityPlot(residuals(jmax_tri_nosize))
shapiro.test(residuals(jmax_tri_nosize))
outlierTest(jmax_tri_nosize)

# Model output
summary(jmax_tri_nosize)
Anova(jmax_tri_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_jmax_predicted_nosize <- data.frame(
  spp = "Tri",
  emmeans(jmax_tri_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(100, 196, 1)))) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax_predicted = emmean,
                jmax_SE = SE,
                jmax_LCL = lower.CL,
                jmax_UCL = upper.CL)

#####################################################################
# Jmax25 - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
jmax25_tri_nosize <- lmer(sqrt(jmax25) ~ gm.trt * doy + (1|id), 
                        data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(jmax25_tri_nosize)
qqnorm(residuals(jmax25_tri_nosize))
qqline(residuals(jmax25_tri_nosize))
densityPlot(residuals(jmax25_tri_nosize))
shapiro.test(residuals(jmax25_tri_nosize))
outlierTest(jmax25_tri_nosize)

# Model output
summary(jmax25_tri_nosize)
Anova(jmax25_tri_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_jmax25_predicted_nosize <- data.frame(
  spp = "Tri",
  emmeans(jmax25_tri_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(100, 196, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax25_predicted = response,
                jmax25_SE = SE,
                jmax25_LCL = lower.CL,
                jmax25_UCL = upper.CL)

#####################################################################
# Ci - Tri (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
ci_tri_nosize <- lmer(ci ~ gm.trt * doy + (1|id), 
                    data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(ci_tri_nosize)
qqnorm(residuals(ci_tri_nosize))
qqline(residuals(ci_tri_nosize))
densityPlot(residuals(ci_tri_nosize))
shapiro.test(residuals(ci_tri_nosize))
outlierTest(ci_tri_nosize)

# Model output
summary(ci_tri_nosize)
Anova(ci_tri_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
tri_ci_predicted_nosize <- data.frame(
  spp = "Tri",
  emmeans(ci_tri_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(100, 196, 1)))) %>%
  dplyr::select(gm.trt:doy, spp,
                ci_predicted = emmean,
                ci_SE = SE,
                ci_LCL = lower.CL,
                ci_UCL = upper.CL)

#####################################################################
# Write data frame that compiles Trillium Ci, Vcmax, Jmax predictions
# across range in DOY accounting for gm.trt and size class intervals
#####################################################################
tri_predictions_nosize <- tri_ci_predicted_nosize %>%
  full_join(tri_vcmax_predicted_nosize) %>%
  full_join(tri_vcmax25_predicted_nosize) %>%
  full_join(tri_jmax_predicted_nosize) %>%
  full_join(tri_jmax25_predicted_nosize) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("ambient", "weeded")))

# Visualize patterns
ggplot(data = tri_predictions_nosize, aes(x = doy, y = vcmax25_predicted)) +
  geom_smooth(method = "loess", aes(color = gm.trt)) +
  geom_ribbon(aes(ymin = vcmax25_LCL,
                  ymax = vcmax25_UCL, 
                  fill = gm.trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Vcmax") +
  theme_bw(base_size = 18)

ggplot(data = tri_predictions_nosize, aes(x = doy, y = jmax25_predicted)) +
  geom_smooth(method = "loess", aes(color = gm.trt)) +
  geom_ribbon(aes(ymin = jmax25_LCL,
                  ymax = jmax25_UCL, 
                  fill = gm.trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Jmax25") +
  theme_bw(base_size = 18)

ggplot(data = tri_predictions_nosize, aes(x = doy, y = ci_predicted)) +
  geom_smooth(method = "loess", aes(color = gm.trt)) +
  geom_ribbon(aes(ymin = ci_LCL, ymax = ci_UCL, fill = gm.trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(200, 400), breaks = seq(200, 400, 100)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Ci") +
  theme_bw(base_size = 18)

#####################################################################
# Vcmax - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
vcmax_mai_nosize <- lmer(log(vcmax) ~ gm.trt * doy + (1|id), 
                       data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(vcmax_mai_nosize)
qqnorm(residuals(vcmax_mai_nosize))
qqline(residuals(vcmax_mai_nosize))
densityPlot(residuals(vcmax_mai_nosize))
shapiro.test(residuals(vcmax_mai_nosize))
outlierTest(vcmax_mai_nosize)

# Model output
summary(vcmax_mai_nosize)
Anova(vcmax_mai_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_vcmax_predicted_nosize <- data.frame(
  spp = "Mai",
  emmeans(vcmax_mai_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax_predicted = response,
                vcmax_SE = SE,
                vcmax_LCL = lower.CL,
                vcmax_UCL = upper.CL)

#####################################################################
# Vcmax25 - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
vcmax25_mai_nosize <- lmer(log(vcmax25) ~ gm.trt * doy + (1|id), 
                         data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(vcmax25_mai_nosize)
qqnorm(residuals(vcmax25_mai_nosize))
qqline(residuals(vcmax25_mai_nosize))
densityPlot(residuals(vcmax25_mai_nosize))
shapiro.test(residuals(vcmax25_mai_nosize))
outlierTest(vcmax25_mai_nosize)

# Model output
summary(vcmax25_mai_nosize)
Anova(vcmax25_mai_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_vcmax25_predicted_nosize <- data.frame(
  spp = "Mai",
  emmeans(vcmax25_mai_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                vcmax25_predicted = response,
                vcmax25_SE = SE,
                vcmax25_LCL = lower.CL,
                vcmax25_UCL = upper.CL)


#####################################################################
# Jmax - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
jmax_mai_nosize <- lmer(log(jmax) ~ gm.trt * doy + (1|id), 
                      data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(jmax_mai_nosize)
qqnorm(residuals(jmax_mai_nosize))
qqline(residuals(jmax_mai_nosize))
densityPlot(residuals(jmax_mai_nosize))
shapiro.test(residuals(jmax_mai_nosize))
outlierTest(jmax_mai_nosize)

# Model output
summary(jmax_mai_nosize)
Anova(jmax_mai_nosize)

# Write data.frame that predicts Jmax for a given day, plot, and size class
mai_jmax_predicted_nosize <- data.frame(
  spp = "Mai",
  emmeans(jmax_mai_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax_predicted = response,
                jmax_SE = SE,
                jmax_LCL = lower.CL,
                jmax_UCL = upper.CL)

#####################################################################
# Jmax25 - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
jmax25_mai_nosize <- lmer(log(jmax25) ~ gm.trt * doy + (1|id), 
                        data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(jmax25_mai_nosize)
qqnorm(residuals(jmax25_mai_nosize))
qqline(residuals(jmax25_mai_nosize))
densityPlot(residuals(jmax25_mai_nosize))
shapiro.test(residuals(jmax25_mai_nosize))
outlierTest(jmax25_mai_nosize)

# Model output
summary(jmax25_mai_nosize)
Anova(jmax25_mai_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_jmax25_predicted_nosize <- data.frame(
  spp = "Mai",
  emmeans(jmax25_mai_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(115, 250, 1)),
          type = "response")) %>%
  dplyr::select(gm.trt:doy, spp,
                jmax25_predicted = response,
                jmax25_SE = SE,
                jmax25_LCL = lower.CL,
                jmax25_UCL = upper.CL)

#####################################################################
# Ci - Mai (accounting for size classes where # of inds are equally
# divided across size distribution)
#####################################################################
photo_traits$ci[c(171, 393, 850)] <- NA

ci_mai_nosize <- lmer(ci ~ gm.trt * doy + (1|id), 
                    data = subset(photo_traits, spp == "Mai" & ci > 0 & ci < 450))

# Check model assumptions
plot(ci_mai_nosize)
qqnorm(residuals(ci_mai_nosize))
qqline(residuals(ci_mai_nosize))
densityPlot(residuals(ci_mai_nosize))
shapiro.test(residuals(ci_mai_nosize))
outlierTest(ci_mai_nosize)

# Model output
summary(ci_mai_nosize)
Anova(ci_mai_nosize)

# Write data.frame that predicts Vcmax for a given day, plot, and size class
mai_ci_predicted_nosize <- data.frame(
  spp = "Mai",
  emmeans(ci_mai_nosize,
          ~gm.trt,
          "doy",
          at = list(doy = seq(115, 250, 1)))) %>%
  dplyr::select(gm.trt:doy, spp,
                ci_predicted = emmean,
                ci_SE = SE,
                ci_LCL = lower.CL,
                ci_UCL = upper.CL)

#####################################################################
# Write data frame that compiles Maianthemum Ci, Vcmax, Jmax 
# predictions across range in DOY accounting for gm.trt and size 
# class intervals
#####################################################################

mai_predictions_nosize <- mai_ci_predicted_nosize %>%
  full_join(mai_vcmax_predicted_nosize) %>%
  full_join(mai_vcmax25_predicted_nosize) %>%
  full_join(mai_jmax_predicted_nosize) %>%
  full_join(mai_jmax25_predicted_nosize) %>%
  mutate(gm.trt = factor(gm.trt,
                           levels = c("ambient", "weeded")))

# Visualize patterns
ggplot(data = mai_predictions_nosize, aes(x = doy, y = vcmax25_predicted)) +
  geom_smooth(method = "loess", aes(color = gm.trt)) +
  geom_ribbon(aes(ymin = vcmax25_LCL,
                  ymax = vcmax25_UCL, 
                  fill = gm.trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Vcmax25") +
  theme_bw(base_size = 18)

ggplot(data = mai_predictions_nosize, aes(x = doy, y = jmax25_predicted)) +
  geom_smooth(method = "loess", aes(color = gm.trt)) +
  geom_ribbon(aes(ymin = jmax25_LCL,
                  ymax = jmax25_UCL, 
                  fill = gm.trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Jmax25") +
  theme_bw(base_size = 18)

ggplot(data = mai_predictions_nosize, aes(x = doy, y = ci_predicted)) +
  geom_smooth(method = "loess", aes(color = gm.trt)) +
  geom_ribbon(aes(ymin = ci_LCL,
                  ymax = ci_UCL, 
                  fill = gm.trt), 
              alpha = 0.1) +
  scale_y_continuous(limits = c(200, 400), breaks = seq(200, 400, 100)) +
  labs(x = "Day of Year",
       y = "Predicted Trillium Ci") +
  theme_bw(base_size = 18)

#####################################################################
# Combine model predictions for Trillium and Maianthemum responses to
# gm treatment across growing season and size classes
#####################################################################

# Merge data frames and rearrange columns
combined_predictions_nosize <- tri_predictions_nosize %>%
  full_join(mai_predictions_nosize) %>%
  dplyr::select(gm_trt = gm.trt, everything())

# Write combined_predictions to .csv
write.csv(combined_predictions_nosize, 
          "../data/c_budget/TT24_budget_params_without_size.csv",
          row.names = F)




