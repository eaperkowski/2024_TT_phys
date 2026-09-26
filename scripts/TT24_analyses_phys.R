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
library(performance)

#####################################################################
# When is physiology acclimated to canopy closure?
#####################################################################
## Using a rolling 10-day mean approach to determine when 10-day
## mean maximum (PARmax10), average PAR (PARavg10), and average
## daytime PAR (PARavg10_day) values cease to significantly differ. 
##
## Physiology is acclimated to canopy closure after at least 5
## consecutive days where there is no difference in 10-day mean and 
## 11-20 day mean light availability. This is true for all light
## availability metrics, but is mostly relevant for average daytime
## PAR.

# Load daily weather station dataset
wx_daily <- read.csv("../data/TT24_daily_weather_summary.csv") %>% filter(doy > 100)
head(wx_daily)
n <- nrow(wx_daily)

# Create empty data frame used in for-loop. Data frame includes
# current day's 10-day max PAR and mean PAR value, the 
# previous day's 10-day max PAR and mean PAR value, and the 
# test summary statistics comparing current and previous day 
# values
df_par_comps <- data.frame(doy = wx_daily$doy,
                           parmax_this10 = wx_daily$parmax10,
                           parmax_last10 = NA,
                           parmax_pval = NA,
                           paravg_this10 = wx_daily$paravg10,
                           paravg_last10 = NA,
                           paravg_pval = NA,
                           paravg_day_this10 = wx_daily$paravg10_day,
                           paravg_day_last10 = NA,
                           paravg_day_pval = NA)

# Previous 10-day window means = rolling means from 10 days earlier
df_par_comps$parmax_last10 <- c(rep(NA, 10), wx_daily$parmax10[seq_len(n - 10)])
df_par_comps$paravg_last10 <- c(rep(NA, 10), wx_daily$paravg10[seq_len(n - 10)])
df_par_comps$paravg_day_last10 <- c(rep(NA, 10), wx_daily$paravg10_day[seq_len(n - 10)])

# For-loop to determine whether rolling windows are different from each other
for (i in seq_len(n)) {
  k <- wx_daily$doy[i]
  
  this_max <- wx_daily$par_max[wx_daily$doy >= (k - 9)  & wx_daily$doy <= k]
  last_max <- wx_daily$par_max[wx_daily$doy >= (k - 19) & wx_daily$doy <= (k - 10)]
  this_avg <- wx_daily$par_mean[wx_daily$doy >= (k - 9)  & wx_daily$doy <= k]
  last_avg <- wx_daily$par_mean[wx_daily$doy >= (k - 19) & wx_daily$doy <= (k - 10)]
  this_avg_day <- wx_daily$par_mean_day[wx_daily$doy >= (k - 9)  & wx_daily$doy <= k]
  last_avg_day <- wx_daily$par_mean_day[wx_daily$doy >= (k - 19) & wx_daily$doy <= (k - 10)]
  
  if (length(this_max) >= 10 && length(last_max) >= 10) {
    df_par_comps$parmax_pval[i] <- t.test(this_max, last_max)$p.value
  }
  if (length(this_avg) >= 10 && length(last_avg) >= 10) {
    df_par_comps$paravg_pval[i] <- t.test(this_avg, last_avg)$p.value
  }
  
  if (length(this_avg_day) >= 10 && length(last_avg_day) >= 10) {
    df_par_comps$paravg_day_pval[i] <- t.test(this_avg_day, last_avg_day)$p.value
  }}

# Boolean dataframe to note when PAR values are greater than 0.05
df_par_comps_bool <- df_par_comps %>%
  mutate(parmax_pval = round(parmax_pval, 3),
         paravg_pval = round(paravg_pval, 3),
         paravg_day_pval = round(paravg_day_pval, 3)) %>%
  dplyr::select(doy, parmax_pval, paravg_pval, paravg_day_pval) %>%
  mutate(parmax_differ = ifelse(parmax_pval < 0.05, TRUE, FALSE),
         paravg_differ = ifelse(paravg_pval < 0.05, TRUE, FALSE),
         paravg_day_differ = ifelse(paravg_day_pval < 0.05, TRUE, FALSE))

# Make a plot 
png("../drafts/figs/TT24_canopy_closure_plot.png", 
    width = 8, height = 6, units = "in", res = 600)
ggplot() +
  geom_rect(aes(xmin = 124.5, xmax = 200, ymin = 0, ymax = Inf),
            fill = "lightgray", alpha = 0.3) +
  geom_vline(xintercept = 140, linewidth = 1, linetype = "dashed") +
  geom_bar(data = subset(wx_daily, doy > 110 & doy < 200), 
           aes(x = doy, y = par_mean_day), stat = "identity", alpha = 0.3) +
  geom_line(data = subset(df_par_comps, doy < 200),
            aes(x = doy, y = paravg_day_this10,
                color = "current"),
            linewidth = 2) +
  geom_line(data = subset(df_par_comps, doy < 200),
            aes(x = doy, y = paravg_day_last10,
                color = "previous"),
            linewidth = 2) +
  scale_x_continuous(limits = c(110, 200), breaks = seq(110, 200, 30)) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  scale_color_manual(values = c("darkblue", "darkred"),
                     labels = c("Current 10-day mean",
                                "Prior 10-day mean")) +
  labs(x = "Day of year",
       y = expression(bold("Mean daytime PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       color = "") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.position = c(0.95, 1.05),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
dev.off()

#####################################################################
# Read in datasets
#####################################################################

# Knot at 140 days (threshold where 11-20 day PAR does not
# differ from 1-10 day PAR)
knot <- 140

# Read .csv file for seasonal carbon budget
cbudget <- read.csv("../data/c_budget/TT24_seasonal_c_budget.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "ambient")),
         plot = factor(plot, levels = c("3", "5", "6")))

# Write dataframe with ID, species, and n_meas to merge with photo_traits
# object
nmeas <- cbudget %>% distinct(id, spp, .keep_all = TRUE) %>%
  dplyr::select(id, spp, n_meas)

# Read .csv file containing gas exchange data
photo_traits <- read.csv("../data/TT24_photo_traits.csv") %>%
  full_join(nmeas, by = c("id", "spp")) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "ambient")),
         plot = factor(plot, levels = c("3", "5", "6")),
         id = factor(id),
         after = pmax(doy - knot, 0))

# Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

# Color palettes
gm.colors <- c("#00B2BE", "#F1B700")

#####################################################################
# Anet - Tri
#####################################################################
anet_tri <- lmer(anet ~ gm.trt * (doy + after) + plot + (1 | id),
                 data = subset(photo_traits, 
                               spp == "Tri" & anet > 0 & n_meas >= 3))

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
performance(anet_tri)

# Pairwise comparisons
test(emtrends(anet_tri, pairwise~gm.trt, "doy"))
test(emtrends(anet_tri, ~1, "doy"))
test(emtrends(anet_tri, ~1, "after"))
emmeans(anet_tri, pairwise~plot)


# Plot prep
# paired doy + after grid (do not let emmeans cross them)
doy_grid_tri <- 100:170

anet_tri_results <- lapply(doy_grid_tri, function(d) {
  out <- as.data.frame(
    emmeans(anet_tri, ~ gm.trt,
            at = list(doy = d, after = pmax(d - knot, 0))))
  out$doy <- d
  out}) %>% 
  bind_rows()

# Plot
anet_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & anet > 0 & n_meas >= 3), 
                        aes(x = doy, y = anet, fill = gm.trt)) +
  geom_rect(aes(xmin = 140, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E5E5E5") +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
  geom_ribbon(data = anet_tri_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt),
              alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = anet_tri_results,
            aes(x = doy, y = emmean, color = gm.trt),
            linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Day of year",
       y = expression(bold("A"["net"]*" ("*mu*"mol"*" m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment")),
       color = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
anet_tri_plot

#####################################################################
# Anet - Mai
#####################################################################
photo_traits$anet[c(139, 579)] <- NA

anet_mai <- lmer(anet ~ gm.trt * (doy + after) + plot + (1 | id),
                 data = subset(photo_traits, spp == "Mai" & anet > 0 & n_meas >= 4 & doy < 210))

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
performance(anet_mai)

# Pairwise comparisons
test(emtrends(anet_mai, ~gm.trt, "doy"))
test(emtrends(anet_mai, ~gm.trt, "after"))
test(emtrends(anet_mai, ~1, "doy"))
test(emtrends(anet_mai, ~1, "after"))
emmeans(anet_mai, pairwise~gm.trt)
emmeans(anet_mai, pairwise~plot)

# Plot prep
doy_grid_mai <- 119:200

anet_mai_results <- lapply(doy_grid_mai, function(d) {
  out <- as.data.frame(
    emmeans(anet_mai, ~ gm.trt,
            at = list(doy = d, after = pmax(d - knot, 0))))
  out$doy <- d
  out}) %>% 
  bind_rows()

# Plot
anet_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & 
                                        !is.na(gm.trt) & anet > 0 & n_meas >= 4), 
                        aes(x = doy, y = anet, fill = gm.trt)) +
  geom_rect(aes(xmin = 140, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E5E5E5") +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
  geom_ribbon(data = anet_mai_results,
              aes(x = doy, y = emmean, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt),
              alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = anet_mai_results,
            aes(x = doy, y = emmean, color = gm.trt),
            linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(118, 202), breaks = seq(120, 200, 20)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
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
photo_traits$vcmax25[c(127, 302, 472, 545)] <- NA

vcmax25_tri <- lmer(sqrt(vcmax25) ~ gm.trt * (doy + after) + plot + (1 | id),
                 data = subset(photo_traits, spp == "Tri" & anet > 0 & 
                                 n_meas >= 3 & id != "5495"))

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

# Plot prep
# paired doy + after grid (do not let emmeans cross them)
doy_grid_tri <- 100:170

vcmax25_tri_results <- lapply(doy_grid_tri, function(d) {
  out <- as.data.frame(
    emmeans(vcmax25_tri, ~ gm.trt,
            at = list(doy = d, after = pmax(d - knot, 0)),
            type = "response"))
  out$doy <- d
  out}) %>% 
  bind_rows()

# Plot
vcmax25_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & 
                                           !is.na(gm.trt) & anet > 0 &
                                           n_meas >= 3), 
                           aes(x = doy, y = vcmax25, fill = gm.trt)) +
  geom_rect(aes(xmin = 140, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E5E5E5") +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
  geom_ribbon(data = vcmax25_tri_results,
              aes(x = doy, y = response, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt),
              alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = vcmax25_tri_results,
            aes(x = doy, y = response, color = gm.trt),
            linewidth = 1.5, inherit.aes = FALSE) +
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
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
vcmax25_tri_plot

#####################################################################
# Vcmax25 - Mai
#####################################################################
vcmax25_mai <- lmer(log(vcmax25) ~ gm.trt * (doy + after) + plot + (1 | id),
                    data = subset(photo_traits, spp == "Mai" & n_meas >= 4 &
                                    anet > 0 & vcmax25 < 180))

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
test(emtrends(vcmax25_mai, ~1, "doy"))
test(emtrends(vcmax25_mai, ~1, "after"))
emmeans(vcmax25_mai, pairwise~plot)

# Plot prep
vcmax25_mai_results <- lapply(doy_grid_mai, function(d) {
  out <- as.data.frame(
    emmeans(vcmax25_mai, ~ gm.trt,
            at = list(doy = d, after = pmax(d - knot, 0)),
            type = "response"))
  out$doy <- d
  out}) %>% 
  bind_rows()

# Plot
vcmax25_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & 
                                           !is.na(gm.trt) & n_meas >= 4 &
                                           vcmax25 < 180), 
                         aes(x = doy, y = vcmax25, fill = gm.trt)) +
  geom_rect(aes(xmin = 140, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E5E5E5") +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
  geom_ribbon(data = vcmax25_mai_results,
              aes(x = doy, y = response, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt),
              alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = vcmax25_mai_results,
            aes(x = doy, y = response, color = gm.trt),
            linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(118, 202), breaks = seq(120, 200, 20)) +
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
# Jmax25 - Tri
#####################################################################
photo_traits$jmax25[c(127, 302, 472, 545)] <- NA

jmax25_tri <- lmer(sqrt(jmax25) ~ gm.trt * (doy + after) + plot + (1 | id), 
                   data = subset(photo_traits, spp == "Tri" & 
                                   !is.na(gm.trt) & n_meas >= 3))

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

# Plot prep
jmax25_tri_results <- lapply(doy_grid_tri, function(d) {
  out <- as.data.frame(
    emmeans(jmax25_tri, ~ gm.trt,
            at = list(doy = d, after = pmax(d - knot, 0)),
            type = "response"))
  out$doy <- d
  out}) %>% 
  bind_rows()

# Plot
jmax25_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & 
                                          !is.na(gm.trt) & n_meas >= 3), 
                        aes(x = doy, y = jmax25, fill = gm.trt)) +
  geom_rect(aes(xmin = 140, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E5E5E5") +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
  geom_ribbon(data =jmax25_tri_results,
              aes(x = doy, y = response, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt),
              alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = jmax25_tri_results,
            aes(x = doy, y = response, color = gm.trt),
            linewidth = 1.5, inherit.aes = FALSE) +
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
# Jmax25 - Mai
#####################################################################
jmax25_mai <- lmer(log(jmax25) ~ gm.trt * (doy + after) + plot + (1|id), 
                   data = subset(photo_traits, spp == "Mai" & n_meas >= 4 & doy < 220))

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
performance(jmax25_mai)

# Pairwise comparisons
test(emtrends(jmax25_mai, ~1, "doy"))
test(emtrends(jmax25_mai, ~1, "after"))
emmeans(jmax25_mai, pairwise~plot)

# Plot prep
jmax25_mai_results <- lapply(doy_grid_mai, function(d) {
  out <- as.data.frame(
    emmeans(jmax25_mai, ~ gm.trt,
            at = list(doy = d, after = pmax(d - knot, 0)),
            type = "response"))
  out$doy <- d
  out}) %>% 
  bind_rows()

# Plot
jmax25_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt)), 
                          aes(x = doy, y = jmax25, fill = gm.trt)) +
  geom_rect(aes(xmin = 140, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#E5E5E5") +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
  geom_ribbon(data = jmax25_mai_results,
              aes(x = doy, y = response, ymin = lower.CL, 
                  ymax = upper.CL, fill = gm.trt),
              alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = jmax25_mai_results,
            aes(x = doy, y = response, color = gm.trt),
            linewidth = 1.5, inherit.aes = FALSE) +
  scale_fill_manual(values = gm.colors) +
  scale_color_manual(values = gm.colors) +
  scale_x_continuous(limits = c(118, 202), breaks = seq(120, 200, 20)) +
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
# Ci - Tri ignoring size classes
#####################################################################
photo_traits$ci[c(77, 261)] <- NA

gsw_tri <- lmer(log(gsw) ~ gm.trt * (doy + after) + plot + (1 | id), 
               data = subset(photo_traits, spp == "Tri" & n_meas >= 3))

# Check model assumptions
plot(gsw_tri)
qqnorm(residuals(gsw_tri))
qqline(residuals(gsw_tri))
densityPlot(residuals(gsw_tri))
shapiro.test(residuals(gsw_tri))
outlierTest(gsw_tri)

# Model output
summary(gsw_tri)
Anova(gsw_tri)
performance(gsw_tri)

# Pairwise comparisons
test(emtrends(gsw_tri, pairwise~gm.trt, "doy"))
## Stronger positive effect of DOY in weeded treatment



# Plot prep
ci_tri_results <- data.frame(
  emmeans(gsw_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1))))

# Plot
ci_tri_plot <- ggplot(data = subset(photo_traits, spp == "Tri" & !is.na(gm.trt) & ci > 0), 
                      aes(x = doy, y = ci, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.05) +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
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
photo_traits$ci[c(13, 82, 183, 324, 455, 722)] <- NA

ci_mai <- lmer(sqrt(gsw) ~ gm.trt * (doy + after) + plot + (1 | id), 
               data = subset(photo_traits, spp == "Mai" & n_meas >= 4 & ci > 200))

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
performance(ci_mai)

# Plot prep
ci_mai_results <- data.frame(
  emmeans(ci_mai, ~gm.trt, "doy", at = list(doy = seq(120, 240, 1))))

# Plot
ci_mai_plot <- ggplot(data = subset(photo_traits, spp == "Mai" & !is.na(gm.trt) & ci > 0), 
                          aes(x = doy, y = ci, fill = gm.trt)) +
  geom_line(aes(group = id), alpha = 0.05) +
  geom_point(size = 2.5, shape = 21, alpha = 0.2) +
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
# Total net C gain (gC/m2/yr) - Trillium
#####################################################################
cbudget$total_netc_assim[c(24, 71)] <- NA

netCgain_tri <- lmer(total_netc_assim ~ gm.trt + (1 | plot), 
                     data = subset(cbudget, spp == "Tri" & n_meas > 2))


# Check model assumptions
plot(netCgain_tri)
qqnorm(residuals(netCgain_tri))
qqline(residuals(netCgain_tri))
densityPlot(residuals(netCgain_tri))
shapiro.test(residuals(netCgain_tri))
outlierTest(netCgain_tri)

# Model output
summary(netCgain_tri)
Anova(netCgain_tri)
r.squaredGLMM(netCgain_tri)

# Post-hoc comparisons
emmeans(netCgain_tri, pairwise~gm.trt)

# Plot prep
netCgain_tri_prep <- cld(emmeans(netCgain_tri, pairwise~gm.trt), Letters = LETTERS,
                         reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
netCgain_tri_plot <- ggplot(data = subset(cbudget, spp == "Tri" & n_meas > 1),
                               aes(x = gm.trt, y = total_netc_assim, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = netCgain_tri_prep, 
            aes(label = .group, y = 300),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(150, 300), breaks = seq(150, 300, 50)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold("Total C gain (gC m"^"-2"*" yr"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

####################################################################
# Total net C gain (gC/m2/yr) - Maianthemum
#####################################################################
cbudget$total_netc_assim[c(32)] <- NA

netCgain_mai <- lmer(log(total_netc_assim) ~ gm.trt + (1 | plot), 
                     data = subset(cbudget, spp == "Mai" & n_meas > 2))


# Check model assumptions
plot(netCgain_mai)
qqnorm(residuals(netCgain_mai))
qqline(residuals(netCgain_mai))
densityPlot(residuals(netCgain_mai))
shapiro.test(residuals(netCgain_mai))
outlierTest(netCgain_mai)

# Model output
summary(netCgain_mai)
Anova(netCgain_mai)
r.squaredGLMM(netCgain_mai)

# Plot prep
netCgain_mai_prep <- cld(emmeans(netCgain_mai, pairwise~gm.trt), Letters = LETTERS,
                         reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
netCgain_mai_plot <- ggplot(data = subset(cbudget, spp == "Mai" & n_meas > 1),
                            aes(x = gm.trt, y = total_netc_assim, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = netCgain_mai_prep, 
            aes(label = .group, y = 300),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(150, 300), breaks = seq(150, 300, 50)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold("Total C gain (gC m"^"-2"*" yr"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

#####################################################################
# Total net C gain (gC/yr) - Trillium
#####################################################################
cbudget$total_netc_assim_tla[c(46)] <- NA

netCgainTLA_tri <- lmer(total_netc_assim_tla ~ gm.trt + (1 | plot), 
                        data = subset(cbudget, spp == "Tri" & n_meas > 2 & plot != 3))


# Check model assumptions
plot(netCgainTLA_tri)
qqnorm(residuals(netCgainTLA_tri))
qqline(residuals(netCgainTLA_tri))
densityPlot(residuals(netCgainTLA_tri))
shapiro.test(residuals(netCgainTLA_tri))
outlierTest(netCgainTLA_tri)

# Model output
summary(netCgainTLA_tri)
Anova(netCgainTLA_tri)
r.squaredGLMM(netCgainTLA_tri)

# Plot prep
netCgainTLA_tri_prep <- cld(emmeans(netCgainTLA_tri, pairwise~gm.trt), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
netCgainTLA_tri_plot <- ggplot(data = subset(cbudget, spp == "Tri" & n_meas > 2 & plot != 3),
                               aes(x = gm.trt, y = total_netc_assim_tla, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = netCgainTLA_tri_prep, 
            aes(label = .group, y = 30),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold("Total C gain (gC yr"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

#####################################################################
# Total net C gain (gC/yr) - Trillium
#####################################################################
cbudget$total_netc_assim_tla[c(46, 173)] <- NA

netCgainTLA_mai <- lmer(total_netc_assim_tla ~ gm.trt + (1 | plot), 
                        data = subset(cbudget, spp == "Mai" & n_meas > 3))


# Check model assumptions
plot(netCgainTLA_mai)
qqnorm(residuals(netCgainTLA_mai))
qqline(residuals(netCgainTLA_mai))
densityPlot(residuals(netCgainTLA_mai))
shapiro.test(residuals(netCgainTLA_mai))
outlierTest(netCgainTLA_mai)

# Model output
summary(netCgainTLA_mai)
Anova(netCgainTLA_mai)
r.squaredGLMM(netCgainTLA_mai)

# Pairwise comparisons
emmeans(netCgainTLA_mai, pairwise~gm.trt)

# Plot prep
netCgainTLA_mai_prep <- cld(emmeans(netCgainTLA_mai, pairwise~gm.trt), 
                            Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
netCgainTLA_mai_plot <- ggplot(data = subset(cbudget, spp == "Mai" & n_meas > 1),
                        aes(x = gm.trt, y = total_netc_assim_tla, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = netCgainTLA_mai_prep, 
            aes(label = .group, y = 30),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold("Total C gain (gC yr"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

#####################################################################
# Compile plots
#####################################################################
png("../drafts/figs/TT24_temp_standardized_vcmax_jmax.png", 
    width = 16, height = 10, units = "in", res = 600)
ggarrange(vcmax25_tri_plot, jmax25_tri_plot, cica_tri_plot,
          vcmax25_mai_plot, jmax25_mai_plot, cica_mai_plot,
          nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom")
dev.off()

png("../drafts/figs/TT24_cbudget.png", 
    width = 8, height = 10, units = "in", res = 600)
ggarrange(netCgain_tri_plot, netCgain_mai_plot,
          netCgainTLA_tri_plot, netCgainTLA_mai_plot,
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()
