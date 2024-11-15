# TT24_compile_datasheet.R
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
photo_traits <- read.csv("../data/TT24_photo_traits.csv") %>% 
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")))

## Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

## Color palettes
gm.colors <- c("#00B2BE", "#F1B700")
gm.colors.soil <- c("#000000", "#D95F0E")


#####################################################################
# Anet - Tri
#####################################################################
photo_traits$vcmax25[c(281, 296, 452)] <- NA

anet_tri <- lmer(log(anet) ~ gm.trt * doy + (1:doy|id), 
                 data = subset(photo_traits, spp == "Tri" & anet > 0))

# Check model assumptions
plot(anet_tri)
qqnorm(residuals(anet_tri))
qqline(residuals(anet_tri))
densityPlot(residuals(anet_tri))
shapiro.test(residuals(anet_tri))
outlierTest(anet_tri)

# Model output
summary(anet_tri)
Anova(anet_tri)
r.squaredGLMM(anet_tri)

# Pairwise comparisons
test(emtrends(anet_tri, pairwise~gm.trt, "doy"))
# Weeded treatment exhibits stronger reduction in Trillium Anet as
# growing season progresses

# Are there any timepoints during the growing season where garlic mustard
# treatment significantly modifies Jmax25?
emmeans(anet_tri, pairwise~gm.trt, "doy", at = list(doy = seq(105, 230, 5)))


# Plot prep
anet_tri_results <- data.frame(
  emmeans(anet_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1)),
          type = "response"))

# Plot
anet_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt)), 
                         aes(x = doy, y = anet, fill = gm.trt)) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(data = anet_tri_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = anet_tri_results,
              aes(x = doy, y = response, ymin = response - SE, 
                  ymax = response + SE, fill = gm.trt), alpha = 0.25) +
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
# Vcmax25 - Tri
#####################################################################
photo_traits$vcmax25[c(281, 296, 452)] <- NA

vcmax25_tri <- lmer(log(vcmax25) ~ gm.trt * doy + (1:doy|id), 
                    data = subset(photo_traits, spp == "Tri"))

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
# (e.g., doy == 105 and marginally at doy == 110), but this effect rapidly 
# diminishes and ambient GM treatment exhibits significantly greater Vcmax25 
# values at end of Tri sampling period (i.e., doy == 150, marginal effect 
# at doy == 145)

# Plot prep
vcmax25_tri_results <- data.frame(
  emmeans(vcmax25_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1)),
          type = "response"))

# Plot
vcmax_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt)), 
                         aes(x = doy, y = vcmax25, fill = gm.trt)) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = vcmax25_tri_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = vcmax25_tri_results,
              aes(x = doy, y = response, ymin = response - SE, 
                  ymax = response + SE, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(100, 170), breaks = seq(100, 170, 20)) +
  scale_y_continuous(limits = c(0, 170), breaks = seq(0, 160, 40)) +
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
vcmax_tri_plot

#####################################################################
# Vcmax25 - Mai
#####################################################################
photo_traits$vcmax25[95] <- NA

vcmax25_mai <- lmer(log(vcmax25) ~ gm.trt * doy + (1:doy|id), 
                    data = subset(photo_traits, spp == "Mai"))

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

# GM decreases Vcmax more strongly

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
vcmax_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt)), 
                         aes(x = doy, y = vcmax25, fill = gm.trt)) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = vcmax25_mai_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = vcmax25_mai_results,
              aes(x = doy, y = response, ymin = response - SE, 
                  ymax = response + SE, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(120, 240), breaks = seq(120, 240, 30)) +
  scale_y_continuous(limits = c(0, 170), breaks = seq(0, 160, 40)) +
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
vcmax_mai_plot

#####################################################################
# Jmax25 - Tri
#####################################################################
photo_traits$jmax25[c(281, 296, 452)] <- NA

jmax25_tri <- lmer(log(jmax25) ~ gm.trt * doy + (1:doy|id), 
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
jmax_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt)), 
                        aes(x = doy, y = jmax25, fill = gm.trt)) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = jmax25_tri_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = jmax25_tri_results,
              aes(x = doy, y = response, ymin = response - SE, 
                  ymax = response + SE, fill = gm.trt), alpha = 0.25) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
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

jmax_tri_plot

#####################################################################
# Jmax25 - Mai
#####################################################################
photo_traits$jmax25[9] <- NA

jmax25_mai <- lmer(log(jmax25) ~ gm.trt * doy + (1:doy|id), 
                   data = subset(photo_traits, spp == "Mai"))

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
jmax_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt)), 
                          aes(x = doy, y = jmax25, fill = gm.trt)) +
  geom_point(size = 2.5, shape = 21, alpha = 0.7) +
  geom_smooth(data = jmax25_mai_results,
              aes(x = doy, y = response, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = jmax25_mai_results,
              aes(x = doy, y = response, ymin = response - SE, 
                  ymax = response + SE, fill = gm.trt), alpha = 0.25) +
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
jmax_mai_plot

png("../data/plot_test.png", width = 10, height = 10,
    units = "in", res = 600)
ggarrange(vcmax_tri_plot, vcmax_mai_plot,
          jmax_tri_plot, jmax_mai_plot,
          nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
dev.off()





