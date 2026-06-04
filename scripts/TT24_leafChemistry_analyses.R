## TT24_leafChemistry_analyses.R

# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)

# Read function for calculating chi and beta
source("../functions/calc_chi.R")
source("../functions/calc_beta.R")

# Read leaf chemistry data sheet
leaf_chemistry <- read.csv("../data/TT24_leaf_chemistry_data.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "ambient")))

## Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

## Color palettes
gm.colors <- c("#00B2BE", "#F1B700")

# ------------------------------------------------------------
# Chi analyses - Tri
# ------------------------------------------------------------

# Model
chi_tri <- lmer(chi ~ gm.trt + (1 | plot), 
                 data = subset(leaf_chemistry, spp == "Tri" & !is.na(subplot)))

# Check normality assumptions
plot(chi_tri)
qqnorm(residuals(chi_tri))
qqline(residuals(chi_tri))
hist(residuals(chi_tri))
shapiro.test(residuals(chi_tri))
outlierTest(chi_tri)

# Model output
summary(chi_tri)
Anova(chi_tri) # No GM treatment effect

# Plot prep
chi_tri_prep <- cld(emmeans(chi_tri, pairwise~gm.trt),
                    reversed = TRUE, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
chi_tri_plot <- ggplot(data = subset(leaf_chemistry, spp == "Tri" & !is.na(subplot)),
                       aes(x = gm.trt, y = chi, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = chi_tri_prep, 
            aes(label = .group, y = 1),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0.8, 1), breaks = seq(0.8, 1, 0.05)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold(chi*" (unitless)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

# ------------------------------------------------------------
# Chi analyses - Mai
# ------------------------------------------------------------

# Model
chi_mai <- lmer(chi ~ gm.trt + (1 | plot), 
                data = subset(leaf_chemistry, spp == "Mai" & !is.na(subplot)))

# Check normality assumptions
plot(chi_mai)
qqnorm(residuals(chi_mai))
qqline(residuals(chi_mai))
hist(residuals(chi_mai))
shapiro.test(residuals(chi_mai))
outlierTest(chi_mai)

# Model output
summary(chi_mai)
Anova(chi_mai) # Marginal GM treatment effect

# Post-hoc comparisons
emmeans(chi_mai, pairwise~gm.trt) 
## Chi is marginally greater in weeded treatment 
## (i.e., stomatal conductance is greater)

# Plot prep
chi_mai_prep <- cld(emmeans(chi_mai, pairwise~gm.trt),
                     reversed = TRUE, Letters = LETTERS, alpha = 0.1) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
chi_mai_plot <- ggplot(data = subset(leaf_chemistry, spp == "Mai" & !is.na(subplot)),
                        aes(x = gm.trt, y = chi, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = chi_mai_prep, 
            aes(label = .group, y = 1),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0.8, 1), breaks = seq(0.8, 1, 0.05)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold(chi*" (unitless)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

# ------------------------------------------------------------
# Beta analyses - Tri
# ------------------------------------------------------------

# Remove outliers
leaf_chemistry$beta[c(17, 29, 75)] <- NA

# Model
beta_tri <- lmer(sqrt(beta) ~ gm.trt + (1 | plot), 
                data = subset(leaf_chemistry, spp == "Tri" & !is.na(subplot)))

# Check normality assumptions
plot(beta_tri)
qqnorm(residuals(beta_tri))
qqline(residuals(beta_tri))
hist(residuals(beta_tri))
shapiro.test(residuals(beta_tri))
outlierTest(beta_tri)

# Model output
summary(beta_tri)
Anova(beta_tri) # No GM treatment effect

# Plot prep
beta_tri_prep <- cld(emmeans(beta_tri, pairwise~gm.trt),
                     reversed = TRUE, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

beta_tri_plot <-ggplot(data = subset(leaf_chemistry, spp == "Tri" & !is.na(subplot)),
       aes(x = gm.trt, y = beta, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = beta_tri_prep, 
            aes(label = .group, y = 3),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(beta)) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

# ------------------------------------------------------------
# Beta analyses - Mai
# ------------------------------------------------------------

# Remove outliers
leaf_chemistry$beta[c(9, 118, 135, 238)] <- NA

# Model
beta_mai <- lmer(log(beta) ~ gm.trt + (1 | plot), 
                 data = subset(leaf_chemistry, spp == "Mai" & 
                                 !is.na(subplot) & beta < 5000))

# Check normality assumptions
plot(beta_mai)
qqnorm(residuals(beta_mai))
qqline(residuals(beta_mai))
hist(residuals(beta_mai))
shapiro.test(residuals(beta_mai))
outlierTest(beta_mai)

# Model output
summary(beta_mai)
Anova(beta_mai) # GM treatment effect

# Post-hoc comparisons
emmeans(beta_mai, pairwise~gm.trt, type = "response") 
## Beta is marginally greater in weeded treatment 
## (likely driven by increase in cost to acquire water in ambient
## treatment)

# Plot prep
beta_mai_prep <- cld(emmeans(beta_mai, pairwise~gm.trt),
    reversed = TRUE, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
beta_mai_plot <- ggplot(data = subset(leaf_chemistry, spp == "Mai" & !is.na(subplot) & beta < 5000),
       aes(x = gm.trt, y = beta, fill = gm.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3, shape = 21) +
  geom_text(data = beta_mai_prep, 
            aes(label = .group, y = 6),
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  scale_fill_manual(values = gm.colors) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold(beta*" (unitless)"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))

# ------------------------------------------------------------
# Make plot
# ------------------------------------------------------------
png("../drafts/figs/TT24_leaf_chemistry.png", height = 10, width = 8,
    units = "in", res = 600)
ggarrange(chi_tri_plot, beta_tri_plot,
          chi_mai_plot, beta_mai_plot,
          nrow = 2, ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)"),
          font.label = list(size = 18))
dev.off()


