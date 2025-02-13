# Libraries
library(tidyverse)
library(lme4)
library(emmeans)
library(car)
library(ggpubr)
library(mgcv)

# Read in C budget data
c_budget <- read.csv("../data/c_budget/TT24_cBudget_data.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c("ambient", "weeded")),
         id = factor(id))

# Double check C budget data are formulated correctly
head(c_budget)
## Looks good

# Some plot aesthetics
gm.colors <- c("#F1B700", "#00B2BE")
facet_lab_tri <- "Trillium spp."
names(facet_lab_tri) <- "Tri"
facet_lab_mai <- "M. racemosum"
names(facet_lab_mai) <- "Mai"

# Separate species into individual data frames
tri_data <- filter(c_budget, spp == "Tri")
mai_data <- filter(c_budget, spp == "Mai")

#####################################################################
# Trillium GAM for net daily C assimilation (g C / m2 / day)
#####################################################################
daily_netC_tri <- gam(daily_netC_assim ~ gm.trt +
                        s(doy, k = 3) +
                        s(doy, by = gm.trt, k = 3) +
                        s(id, bs = "re"),
                      data = subset(c_budget, spp == "Tri" & doy < 170),
                      model = "REML")

# Model results
summary(daily_netC_tri)

# Some patterns
test(emtrends(daily_netC_tri, ~1, "doy"))
test(emtrends(daily_netC_tri, pairwise~gm.trt, "doy"))

# Prepare regression line for plot
daily_netC_tri_regline <- data.frame(
  emmeans(
    daily_netC_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1))))

# Plot
daily_netC_tri_plot <- ggplot(
  data = subset(c_budget, spp == "Tri" & doy < 170 & !is.na(gm.trt)),
  aes(x = doy, y = daily_netC_assim)) +
  geom_point(aes(fill = gm.trt), shape = 21, alpha = 0.7, size = 2) +
  geom_ribbon(data = daily_netC_tri_regline,
              aes(x = doy, y = emmean, ymin = lower.CL, ymax = upper.CL,
                  fill = gm.trt), alpha = 0.5) +
  geom_line(data = daily_netC_tri_regline,
            aes(x = doy, y = emmean, color = gm.trt), linewidth = 1) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(-1, 12), breaks = seq(0, 12, 4)) +
  labs(x = "Day of year",
       y = expression(bold("Net daily C assimilation (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())
  
#####################################################################
# Trillium GAM for net daily C assimilation accounting for leaf/plant
# size (g C / day)
#####################################################################
daily_netC_tla_tri <- gam(daily_netC_assim_tla ~ gm.trt + 
                            s(doy, k = 3) +
                            s(doy, by = gm.trt, k = 3) + 
                            s(id, bs = "re"), 
                          data = subset(c_budget, spp == "Tri" & doy < 170),
                        model = "REML")

# Model results
summary(daily_netC_tla_tri)
anova(daily_netC_tla_tri)

# Let's look at some trends
test(emtrends(daily_netC_tla_tri, ~1, "doy"))
test(emtrends(daily_netC_tla_tri, pairwise~gm.trt, "doy"))

# Prepare regression line for plot
daily_netC_tla_tri_regline <- data.frame(
  emmeans(
    daily_netC_tla_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1))))

# Plot
daily_netC_tla_tri_plot <- ggplot(
  data = subset(c_budget, spp == "Tri" & doy < 170 & !is.na(gm.trt)),
  aes(x = doy, y = daily_netC_assim_tla)) +
  geom_point(aes(fill = gm.trt), shape = 21, alpha = 0.7, size = 2) +
  geom_ribbon(data = daily_netC_tla_tri_regline,
              aes(x = doy, y = emmean, ymin = lower.CL, ymax = upper.CL,
                  fill = gm.trt), alpha = 0.5) +
  geom_line(data = daily_netC_tla_tri_regline,
            aes(x = doy, y = emmean, color = gm.trt), linewidth = 1) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(-0.02, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Day of year",
       y = expression(bold("Total net daily C assimilation (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())

#####################################################################
# Maianthemum GAM for net daily C assimilation (g C / m2 / day)
#####################################################################
daily_netC_mai <- lmer(log(daily_netC_assim) ~ gm.trt * doy + (1 | id),
                       data = subset(
                         c_budget, spp == "Mai" & daily_netC_assim > 0))

# Model results
summary(daily_netC_mai)

# Some patterns
test(emtrends(daily_netC_mai, ~1, "doy"))
test(emtrends(daily_netC_mai, pairwise~gm.trt, "doy"))

# Prepare regression line for plot
daily_netC_mai_regline <- data.frame(
  emmeans(
    daily_netC_mai, ~gm.trt, "doy", at = list(doy = seq(120, 240, 1)),
    type = "response"))


# Plot
daily_netC_mai_plot <- ggplot(data = subset(c_budget, spp == "Mai" & !is.na(gm.trt) & 
                       daily_netC_assim > 0 & doy > 120 & doy < 240),
       aes(x = doy, y = daily_netC_assim)) +
  geom_point(aes(fill = gm.trt), shape = 21, alpha = 0.7, size = 2) +
  geom_ribbon(data = daily_netC_mai_regline,
              aes(x = doy, y = response, ymin = asymp.LCL, ymax = asymp.UCL,
                  fill = gm.trt), alpha = 0.5) +
  geom_line(data = daily_netC_mai_regline,
            aes(x = doy, y = response, color = gm.trt), linewidth = 1) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  labs(x = "Day of year",
       y = expression(bold("Net daily C assimilation (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())

#####################################################################
# Maianthemum GAM for net daily C assimilation accounting for leaf/plant
# size (g C / day)
#####################################################################
daily_netC_tla_mai <- gam(daily_netC_assim_tla ~ gm.trt + 
                            s(doy, k = 3) +
                            s(doy, by = gm.trt, k = 3) + 
                            s(id, bs = "re"), 
                          data = subset(c_budget, spp == "Mai" & daily_netC_assim_tla > 0),
                          model = "REML")

# Model results
summary(daily_netC_tla_mai)
anova(daily_netC_tla_mai)

# Let's look at some trends
test(emtrends(daily_netC_tla_mai, ~1, "doy"))
test(emtrends(daily_netC_tla_mai, pairwise~gm.trt, "doy"))

# Prepare regression line for plot
daily_netC_tla_mai_regline <- data.frame(
  emmeans(
    daily_netC_tla_mai, ~gm.trt, "doy", at = list(doy = seq(120, 240, 1))))

# Plot
daily_netC_tla_mai_plot <- ggplot(
  data = subset(c_budget, spp == "Mai" & !is.na(gm.trt) & 
                  daily_netC_assim > 0 & doy > 120 & doy < 240),
       aes(x = doy, y = daily_netC_assim_tla)) +
  geom_point(aes(fill = gm.trt), shape = 21, alpha = 0.7, size = 2) +
  geom_ribbon(data = daily_netC_tla_mai_regline,
              aes(x = doy, y = emmean, ymin = lower.CL, ymax = upper.CL,
                  fill = gm.trt), alpha = 0.5) +
  geom_line(data = daily_netC_tla_mai_regline,
            aes(x = doy, y = emmean, color = gm.trt), linewidth = 1) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Day of year",
       y = expression(bold("Total net daily C assimilation (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid.minor.y = element_blank())


#####################################################################
# Make figure!
#####################################################################

png(filename = "../drafts/figs/TT24_cBudget_bothSpecies.png", height = 12,
    width = 12, units = "in", res = 600)
ggarrange(daily_netC_tri_plot, daily_netC_tla_tri_plot, 
          daily_netC_mai_plot, daily_netC_tla_mai_plot, 
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)"), font.label = list(font.size = 18))
dev.off()
