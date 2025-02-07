# TT24_crude_analysis.R
# R script that gives a first stab look at the effects of GM on photosynthetic
# processes in Trillium and Maianthemum
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
photo_traits <- read.csv("../data/TT24_photo_traits.csv")

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
  full_join(intervals) %>%
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")),
         id = factor(id))

#####################################################################
# Anet - Tri
#####################################################################
anet_tri <- gam(anet ~ gm.trt + s(doy, k = 4) +
                  s(doy, by = gm.trt, k = 4) +
                  s(id,  bs = "re"),
                data = subset(photo_traits, spp == "Tri" & anet > 0))

# Check model assumptions
plot(anet_tri)
qqnorm(residuals(anet_tri))
qqline(residuals(anet_tri))
densityPlot(residuals(anet_tri))
shapiro.test(residuals(anet_tri))
outlierTest(anet_tri)

# Model output
summary(anet_tri_gam)
Anova(anet_tri)
r.squaredGLMM(anet_tri)

# Pairwise comparisons
test(emtrends(anet_tri_gam, pairwise~gm.trt, "doy"))
# Weeded treatment exhibits stronger reduction in Trillium Anet as
# growing season progresses

# Plot prep
anet_tri_results <- data.frame(
  emmeans(anet_tri_gam, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1))))

# Plot
anet_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt)), 
                         aes(x = doy, y = anet, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(data = anet_tri_results,
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = anet_tri_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  labs(x = "Day of year",
       y = expression(bold("A"["net"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment")),
       color = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
anet_tri_plot

#####################################################################
# Anet - Mai
#####################################################################
photo_traits$anet[c(110)] <- NA

anet_mai <- gam(anet ~ gm.trt + s(doy, k = 4) +
                  s(doy, by = gm.trt, k = 4) +
                  s(id, plot, bs = "re"),
                data = subset(photo_traits, spp == "Mai" & anet > 0))

# Check model assumptions
plot(anet_mai)
qqnorm(residuals(anet_mai))
qqline(residuals(anet_mai))
densityPlot(residuals(anet_mai))
shapiro.test(residuals(anet_mai))
outlierTest(anet_mai)

# Model output
summary(anet_mai)
Anova(anet_mai)
r.squaredGLMM(anet_mai)

# Pairwise comparisons
test(emtrends(anet_mai, pairwise~gm.trt, "doy"))
# Ambient treatment exhibits stronger reduction in Trillium Anet as
# growing season progresses

# Are there any timepoints during the growing season where garlic mustard
# treatment significantly modifies Jmax25?
emmeans(anet_mai, pairwise~gm.trt, "doy", at = list(doy = seq(105, 240, 5)))

# Plot prep
anet_mai_results <- data.frame(
  emmeans(anet_mai, ~gm.trt, "doy",  at = list(doy = seq(105, 240, 5)),
          type = "response"))

# Plot
anet_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & 
                                        !is.na(gm.trt) & anet > 0), 
                        aes(x = doy, y = anet, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(data = anet_mai_results,
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = anet_mai_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4)) +
  labs(x = "Day of year",
       y = expression(bold("A"["net"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment")),
       color = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
anet_mai_plot

#####################################################################
#####################################################################
# TEMP STANDARDIZED RATES
#####################################################################
#####################################################################

#####################################################################
# Vcmax25 - Tri
#####################################################################
photo_traits$vcmax25[c(281, 296, 452)] <- NA

vcmax25_tri <- gam(vcmax25 ~ gm.trt + s(doy, k = 4) +
                     s(doy, by = gm.trt, k = 4) +
                     s(id, plot, bs = "re"),
                   data = subset(photo_traits, spp == "Tri" & anet > 0))

# Check model assumptions
plot(vcmax25_tri)
qqnorm(residuals(vcmax25_tri))
qqline(residuals(vcmax25_tri))
densityPlot(residuals(vcmax25_tri))
shapiro.test(residuals(vcmax25_tri))
outlierTest(vcmax25_tri)

# Model output
summary(vcmax25_tri)
Anova(vcmax25_tri)
r.squaredGLMM(vcmax25_tri)

# Pairwise comparisons
test(emtrends(vcmax25_tri, pairwise~gm.trt, "doy"))
# Weeded treatment exhibits stronger reduction in Trillium Vcmax25 as growing 
# season progresses

# Are there any timepoints during the growing season where garlic mustard
# treatment significantly modifies Vcmax25?
emmeans(vcmax25_tri, pairwise~gm.trt, "doy", at = list(doy = seq(105, 230, 5)))
# Ambient GM treatment exhibits *reduced* Vcmax25 in the early season
# (doy == 105), but this effect rapidly diminishes and ambient GM treatment 
# exhibits significantly greater Vcmax25 values at end of Tri sampling period 
# (i.e., doy == 165, marginal effect at doy == 155, 160)

# Plot prep
vcmax25_tri_results <- data.frame(
  emmeans(vcmax25_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1)),
          type = "response"))

# Plot
vcmax25_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt)), 
                         aes(x = doy, y = vcmax25, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = vcmax25_tri_results,
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = vcmax25_tri_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(100, 170), breaks = seq(100, 170, 20)) +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 60)) +
  labs(x = "Day of year",
       y = expression(bold("V"["cmax25"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
vcmax25_tri_plot

#####################################################################
# Vcmax25 - Tri including size classes (number_interval to keep number
# of reps within each class the same)
#####################################################################
vcmax25_tri_numbInt <- lmer(
  log(vcmax25) ~ gm.trt * doy * number_interval + (1|id), 
  data = subset(photo_traits, spp == "Tri"))

# Check model assumptions
plot(vcmax25_tri_numbInt)
qqnorm(residuals(vcmax25_tri_numbInt))
qqline(residuals(vcmax25_tri_numbInt))
densityPlot(residuals(vcmax25_tri_numbInt))
shapiro.test(residuals(vcmax25_tri_numbInt))
outlierTest(vcmax25_tri_numbInt)

# Model output
summary(vcmax25_tri_numbInt)
Anova(vcmax25_tri_numbInt)
r.squaredGLMM(vcmax25_tri_numbInt)

# Pairwise comparisons
test(emtrends(vcmax25_tri_numbInt, pairwise~gm.trt, "doy"))
## Slope that explains doy-vcmax25 relationship is significantly more 
## negative in weeded treatment than ambient treatment

emmeans(vcmax25_tri_numbInt, pairwise~gm.trt, "doy", at = list(doy = seq(105, 230, 5)))

cld(emmeans(vcmax25_tri_numbInt, pairwise~gm.trt|number_interval))
# Vcmax25 is greater in ambient plots, but only in lowest size class

#####################################################################
# Vcmax25 - Mai
#####################################################################
photo_traits$vcmax25[c(95, 99, 123)] <- NA

vcmax25_mai <- gam(vcmax25 ~ gm.trt + s(doy, k = 3) +
                     s(doy, by = gm.trt, k = 3) +
                     s(id, plot, bs = "re"),
                   data = subset(photo_traits, spp == "Mai" & anet > 0))

# Check model assumptions
plot(vcmax25_mai)
qqnorm(residuals(vcmax25_mai))
qqline(residuals(vcmax25_mai))
densityPlot(residuals(vcmax25_mai))
shapiro.test(residuals(vcmax25_mai))
outlierTest(vcmax25_mai)

# Model output
summary(vcmax25_mai)
Anova(vcmax25_mai)
r.squaredGLMM(vcmax25_mai)

# Pairwise comparisons
test(emtrends(vcmax25_mai, pairwise~gm.trt, "doy")) 
# GM presence decreases Vcmax more strongly across growing season

# Are there any timepoints during the growing season where garlic mustard
# treatment significantly modifies Vcmax25?
emmeans(vcmax25_mai, pairwise~gm.trt, "doy", at = list(doy = seq(120, 230, 5)))
# GM treatment decreases Vcmax25 starting at doy == 200 (marginal at doy == 195)
# and this effect continues throughout the remainder of the sampling period

# Plot prep
vcmax25_mai_results <- data.frame(
  emmeans(vcmax25_mai, ~gm.trt, "doy", at = list(doy = seq(120, 240, 1)),
          type = "response"))

# Plot
vcmax25_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt)), 
                         aes(x = doy, y = vcmax25, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = vcmax25_mai_results,
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = vcmax25_mai_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(120, 240), breaks = seq(120, 240, 30)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  labs(x = "Day of year",
       y = expression(bold("V"["cmax25"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
vcmax25_mai_plot

#####################################################################
# Vcmax25 - Mai including size classes (number_interval to keep number
# of reps within each class the same)
#####################################################################
photo_traits$vcmax25[c(99,122)] <- NA

vcmax25_mai_numbInt <- lmer(
  log(vcmax25) ~ gm.trt * doy * number_interval + (1|id), 
  data = subset(photo_traits, spp == "Mai"))

# Check model assumptions
plot(vcmax25_mai_numbInt)
qqnorm(residuals(vcmax25_mai_numbInt))
qqline(residuals(vcmax25_mai_numbInt))
densityPlot(residuals(vcmax25_mai_numbInt))
shapiro.test(residuals(vcmax25_mai_numbInt))
outlierTest(vcmax25_mai_numbInt)

# Model output
summary(vcmax25_mai_numbInt)
Anova(vcmax25_mai_numbInt)
r.squaredGLMM(vcmax25_mai_numbInt)

# Pairwise comparisons
test(emtrends(vcmax25_mai_numbInt, ~gm.trt, "doy"))
## Reductions in Vcmax25 across the growing season are significantly
## greater in the ambient treatment

cld(emtrends(vcmax25_mai_numbInt, ~gm.trt|number_interval, "doy"))
## Stronger Vcmax25 reductions across growing season in ambient
## treatment are driven by 554-673cm2 and 673-1120cm2 size classes
## (i.e., GM treatment more strongly influences larger Mai)

#####################################################################
# Jmax25 - Tri
#####################################################################
photo_traits$jmax25[c(281, 296, 452)] <- NA

jmax25_tri <- lmer(log(jmax25) ~ gm.trt * doy + (1 | id) + (1 | plot), 
                   data = subset(photo_traits, spp == "Tri"))

# Check model assumptions
plot(jmax25_tri)
qqnorm(residuals(jmax25_tri))
qqline(residuals(jmax25_tri))
densityPlot(residuals(jmax25_tri))
shapiro.test(residuals(jmax25_tri))
outlierTest(jmax25_tri)

# Model output
summary(jmax25_tri)
Anova(jmax25_tri)
r.squaredGLMM(jmax25_tri)

# Pairwise comparisons
test(emtrends(jmax25_tri, pairwise~gm.trt, "doy"))
# Weeded treatment exhibits stronger reduction in Jmax25 as growing season
# progresses

# Are there any timepoints during the growing season where garlic mustard
# treatment significantly modifies Jmax25?
emmeans(jmax25_tri, pairwise~gm.trt, "doy", at = list(doy = seq(105, 150, 5)))
## Ambient GM treatment has marginally reduced Jmax25 compared to weeded 
## treatment under doy == 105 and 110

# Plot prep
jmax25_tri_results <- data.frame(
  emmeans(jmax25_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1)),
          type = "response"))

# Plot
jmax25_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt)), 
                        aes(x = doy, y = jmax25, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = jmax25_tri_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = jmax25_tri_results,
              aes(x = doy, y = response, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) + 
  labs(x = "Day of year",
       y = expression(bold("J"["max25"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
jmax25_tri_plot

#####################################################################
# Jmax25 - Tri including size classes (number_interval to keep number
# of reps within each class the same)
#####################################################################
jmax25_tri_numbInt <- lmer(
  log(jmax25) ~ gm.trt * doy * number_interval + (1|id), 
  data = subset(photo_traits, spp == "Tri"))

# Check model assumptions
plot(jmax25_tri_numbInt)
qqnorm(residuals(jmax25_tri_numbInt))
qqline(residuals(jmax25_tri_numbInt))
densityPlot(residuals(jmax25_tri_numbInt))
shapiro.test(residuals(jmax25_tri_numbInt))
outlierTest(jmax25_tri_numbInt)

# Model output
summary(jmax25_tri_numbInt)
Anova(jmax25_tri_numbInt)
r.squaredGLMM(jmax25_tri_numbInt)

# Pairwise comparisons
test(emtrends(jmax25_tri_numbInt, pairwise~gm.trt, "doy"))
## Reductions in Jmax25 across the growing season are stronger
## in weeded treatment

#####################################################################
# Jmax25 - Mai
#####################################################################
jmax25_mai <- lmer(log(jmax25) ~ gm.trt * doy + (1|id) + (1 | plot), 
                   data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(jmax25_mai)
qqnorm(residuals(jmax25_mai))
qqline(residuals(jmax25_mai))
densityPlot(residuals(jmax25_mai))
shapiro.test(residuals(jmax25_mai))
outlierTest(jmax25_mai)

# Model output
summary(jmax25_mai)
Anova(jmax25_mai)
r.squaredGLMM(jmax25_mai)

# Pairwise comparisons
test(emtrends(jmax25_mai, pairwise~gm.trt, "doy"))
# Ambient treatment exhibits stronger reduction in Maianthemum Jmax25 across
# the growing season

# Are there any timepoints during the growing season where garlic mustard
# treatment significantly modifies Jmax25?
emmeans(jmax25_mai, pairwise~gm.trt, "doy", at = list(doy = seq(120, 230, 5)))
# Ambient GM treatment has marginally reduced Jmax25 compared to weeded 
# treatment later in growing season (at least marginally significant 
# at all timepoints after doy == 210)

# Plot prep
jmax25_mai_results <- data.frame(
  emmeans(jmax25_mai, ~gm.trt, "doy", at = list(doy = seq(120, 240, 1)),
          type = "response"))

# Plot
jmax25_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt)), 
                          aes(x = doy, y = jmax25, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = jmax25_mai_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = jmax25_mai_results,
              aes(x = doy, y = response, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(120, 240), breaks = seq(120, 240, 30)) +
  scale_y_continuous(limits = c(0, 210), breaks = seq(0, 200, 50)) +
  labs(x = "Day of year",
       y = expression(bold("J"["max25"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
jmax25_mai_plot

#####################################################################
# Jmax25 - Mai including size classes (number_interval to keep number
# of reps within each class the same)
#####################################################################
jmax25_mai_numbInt <- lmer(
  log(jmax25) ~ gm.trt * doy * number_interval + (1|id), 
  data = subset(photo_traits, spp == "Mai"))

# Check model assumptions
plot(jmax25_mai_numbInt)
qqnorm(residuals(jmax25_mai_numbInt))
qqline(residuals(jmax25_mai_numbInt))
densityPlot(residuals(jmax25_mai_numbInt))
shapiro.test(residuals(jmax25_mai_numbInt))
outlierTest(jmax25_mai_numbInt)

# Model output
summary(jmax25_mai_numbInt)
Anova(jmax25_mai_numbInt)
r.squaredGLMM(jmax25_mai_numbInt)

# Pairwise comparisons
test(emtrends(jmax25_mai_numbInt, pairwise~gm.trt, "doy"))
## Reductions in Jmax25 across the growing season are similar
## between gm treatments

cld(emtrends(jmax25_mai_numbInt, ~gm.trt|number_interval, "doy"))
## Jmax25 reductions across growing season are stronger in ambient
## treatment in larger size classes

#####################################################################
# Ci - Tri ignoring size classes
#####################################################################
ci_tri <- lmer(ci ~ gm.trt * doy + (1|id) + (1 | plot), 
               data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(ci_tri)
qqnorm(residuals(ci_tri))
qqline(residuals(ci_tri))
densityPlot(residuals(ci_tri))
shapiro.test(residuals(ci_tri))
outlierTest(ci_tri)

# Model output
summary(ci_tri)
Anova(ci_tri)
r.squaredGLMM(ci_tri)

# Pairwise comparisons
test(emtrends(ci_tri, pairwise~gm.trt, "doy"))
## Stronger positive effect of DOY in weeded treatment

# Plot prep
ci_tri_results <- data.frame(
  emmeans(ci_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1))))

# Plot
ci_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt) & ci > 0), 
                      aes(x = doy, y = ci, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = ci_tri_results,
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = ci_tri_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_y_continuous(limits = c(100, 400), breaks = seq(100, 400, 100)) +
  labs(x = "Day of year",
       y = expression(bold("C"["i"]*" ("*mu*"mol"*" mol"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
ci_tri_plot

#####################################################################
# Ci - Mai ignoring size classes
#####################################################################
photo_traits$ci[c(367, 392, 790, 791, 820)] <- NA
photo_traits$ci[c(171, 393, 450, 850)] <- NA

ci_mai <- lmer(ci ~ gm.trt * doy + (1 | id) + (1 | plot), 
               data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(ci_mai)
qqnorm(residuals(ci_mai))
qqline(residuals(ci_mai))
densityPlot(residuals(ci_mai))
shapiro.test(residuals(ci_mai))
outlierTest(ci_mai)

# Model output
summary(ci_mai)
Anova(ci_mai)
r.squaredGLMM(ci_mai)

# Plot prep
ci_mai_results <- data.frame(
  emmeans(ci_mai, ~gm.trt, "doy", at = list(doy = seq(120, 240, 1))))

# Plot
ci_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt) & ci > 0), 
                          aes(x = doy, y = ci, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = ci_mai_results,
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = ci_mai_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(120, 240), breaks = seq(120, 240, 30)) +
  scale_y_continuous(limits = c(100, 400), breaks = seq(100, 400, 100)) +
  labs(x = "Day of year",
       y = expression(bold("C"["i"]*" ("*mu*"mol"*" mol"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
ci_mai_plot


#####################################################################
# Compile plots
#####################################################################
png("../drafts/figs/TT24_temp_standardized_vcmax_jmax.png", 
    width = 16, height = 10, units = "in", res = 600)
ggarrange(vcmax25_tri_plot, jmax25_tri_plot, ci_tri_plot,
          vcmax25_mai_plot, jmax25_mai_plot, ci_mai_plot,
          nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom")
dev.off()


