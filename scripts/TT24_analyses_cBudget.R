# Libraries
library(tidyverse)
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
ggplot(data = subset(c_budget, spp == "Tri" & doy < 170 & !is.na(gm.trt)),
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
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))
  
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
ggplot(data = subset(c_budget, spp == "Tri" & doy < 170 & !is.na(gm.trt)),
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
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))


#####################################################################
# Trillium GAM for cumulative daily C assimilation (g C / m2 / day)
#####################################################################

cumul_netC_tri <- gam(cumul_netC_assim ~ gm.trt +
                        s(doy, k = 3) +
                        s(doy, by = gm.trt, k = 3) +
                        s(id, bs = "re"),
                      data = subset(c_budget, spp == "Tri" & doy < 170),
                      model = "REML")

# Model results
summary(cumul_netC_tri)

# Some patterns
test(emtrends(cumul_netC_tri, ~1, "doy"))
test(emtrends(cumul_netC_tri, pairwise~gm.trt, "doy"))

# Prepare regression line for plot
daily_netC_tri_regline <- data.frame(
  emmeans(
    cumul_netC_tri, ~gm.trt, "doy", at = list(doy = seq(100, 170, 1))))

# Plot
ggplot(data = subset(c_budget, spp == "Tri" & doy < 170 & !is.na(gm.trt)),
       aes(x = doy, y = cumul_netC_assim_fract)) +
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
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"))












#####################################################################
# Trillium plots
#####################################################################
# Net C assimilation without plant size (gC/m2/d)
tri_netC_assim_plot <- ggplot(data = subset(daily_model_results, 
                                            spp == "Tri" & daily_C_assim > 0 & n_measurements > 3 &
                                              doy < 170),
                              aes(x = doy, y = daily_netC_assim)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 3)) +
  labs(x = "Day of year",
       y = expression(bold("Daily net C assimilation (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# C assimilation without plant size (gC/m2/d)
tri_C_assim_plot <- ggplot(data = subset(daily_model_results, 
                                         spp == "Tri" & daily_C_assim > 0 & n_measurements > 3 &
                                           doy < 170),
                           aes(x = doy, y = daily_C_assim)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 3)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C assimilation (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# Respiration without plant size (gC/m2/d)
tri_C_resp_plot <- ggplot(data = subset(daily_model_results, 
                                        spp == "Tri" & daily_C_assim > 0 & n_measurements > 3 &
                                          doy < 170),
                          aes(x = doy, y = daily_C_resp)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C respiration (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# Net C assimilation accounting for plant size
tri_netC_assim_tla_plot <- ggplot(data = subset(daily_model_results, 
                                                spp == "Tri" & daily_C_assim > 0 & n_measurements > 3 &
                                                  doy < 170),
                                  aes(x = doy, y = daily_netC_assim_tla)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Day of year",
       y = expression(bold("Daily net C assimilation (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# C assimilation accounting for plant size
tri_C_assim_tla_plot <- ggplot(data = subset(daily_model_results, 
                                             spp == "Tri" & daily_C_assim > 0 & n_measurements > 3 &
                                               doy < 170),
                               aes(x = doy, y = daily_C_assim_tla)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C assimilation (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# Respiration accounting for plant size
tri_C_resp_tla_plot <- ggplot(data = subset(daily_model_results, 
                                            spp == "Tri" & daily_C_assim > 0 & n_measurements > 3 &
                                              doy < 170),
                              aes(x = doy, y = daily_C_resp_tla)) +
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C respiration (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

#####################################################################
# Mai plots
#####################################################################
# Net C assimilation without plant size (gC/m2/d)
mai_netC_assim_plot <- ggplot(data = subset(daily_model_results, 
                                            spp == "Mai" & daily_C_assim > 0 & n_measurements > 3 &
                                              doy > 110 & doy < 250),
                              aes(x = doy, y = daily_netC_assim)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(-1, 6), breaks = seq(0, 6, 2)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily net C assimilation (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# C assimilation without plant size (gC/m2/d)
mai_C_assim_plot <- ggplot(data = subset(daily_model_results, 
                                         spp == "Mai" & daily_C_assim > 0 & n_measurements > 3 &
                                           doy > 110 & doy < 250),
                           aes(x = doy, y = daily_C_assim)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C assimilation (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))
mai_C_assim_plot

# Respiration without plant size (gC/m2/d)
mai_C_resp_plot <- ggplot(data = subset(daily_model_results, 
                                        spp == "Mai" & daily_C_assim > 0 & n_measurements > 3 &
                                          doy > 110 & doy < 250),
                          aes(x = doy, y = daily_C_resp)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C respiration (g C m"^"-2"*" d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# Net C assimilation accounting for plant size
mai_netC_assim_tla_plot <- ggplot(data = subset(daily_model_results, 
                                                spp == "Mai" & daily_C_assim > 0 & n_measurements > 3 &
                                                  doy > 110 & doy < 250),
                                  aes(x = doy, y = daily_netC_assim_tla)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily net C assimilation (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# C assimilation accounting for plant size
mai_C_assim_tla_plot <- ggplot(data = subset(daily_model_results, 
                                             spp == "Mai" & daily_C_assim > 0 & n_measurements > 3 &
                                               doy > 110 & doy < 250),
                               aes(x = doy, y = daily_netC_assim_tla)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C assimilation (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))

# Respiration accounting for plant size
mai_C_resp_tla_plot <- ggplot(data = subset(daily_model_results, 
                                            spp == "Mai" & daily_C_assim > 0 & n_measurements > 3 &
                                              doy > 110 & doy < 250),
                              aes(x = doy, y = daily_netC_assim_tla)) +
  geom_rect(aes(xmin = 125, xmax = 250, ymin = -Inf, ymax = Inf),
            alpha = 0.8, fill = "#ECECEC") +
  geom_point(aes(fill = gm.trt), shape = 21, size = 2, alpha = 0.7) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(110, 250), breaks = seq(110, 250, 20)) +
  labs(x = "Day of year",
       y = expression(bold("Daily C respiration (g C d"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment")),
       color = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet_lab_mai)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18))


png("../drafts/figs/TT24_cBudget_mai.png", width = 16, height = 12, units = "in", res = 600)
ggarrange(mai_C_assim_plot, mai_C_resp_plot, mai_netC_assim_plot,
          mai_C_assim_tla_plot, mai_C_resp_tla_plot, mai_netC_assim_tla_plot,
          nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))
dev.off()


