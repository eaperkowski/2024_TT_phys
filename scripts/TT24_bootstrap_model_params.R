#####################################################################
# Libraries and file read in
#####################################################################
# Libraries
library(tidyverse)
library(lme4)
library(ggpubr)

# Read .csv file
photo_traits <- read.csv("../data/TT24_photo_traits.csv") %>% 
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")))

# Define a function for bootstrapping model predictions
predict_function <- function(fit) {
  # Define the new data for prediction
  trt_df <- expand.grid(
    doy = seq(100, 250, by = 1),
    gm.trt = levels(photo_traits$gm.trt)
  )
  
  # Return raw predictions for each bootstrap sample
  predict(fit, newdata = trt_df, re.form = NA)
}

# Make placeholder to merge bootstrapped estimates to
trt_df_placeholder <- expand.grid(
  doy = seq(100, 250, by = 1),
  gm.trt = levels(photo_traits$gm.trt))

# Set seed for reproducibility
set.seed(10)

#####################################################################
# Vcmax25 - Tri (does not account for size... yet)
#####################################################################
photo_traits$vcmax25[c(281, 296, 452)] <- NA

vcmax25_tri <- lmer(log(vcmax25) ~ gm.trt * doy + (1|id), 
                    data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(vcmax25_tri)
qqnorm(residuals(vcmax25_tri))
qqline(residuals(vcmax25_tri))
densityPlot(residuals(vcmax25_tri))
shapiro.test(residuals(vcmax25_tri))
outlierTest(vcmax25_tri)

# Bootstrap model predictions to get 1000 Vcmax values per doy for
# each of gm.trt == "ambient" and gm.trt == "weeded"
vcmax25_tri_boot <- bootMer(
  vcmax25_tri,
  predict_function,
  nsim = 1000,
  type = "parametric")

# Extract bootstrapped predictions
vcmax25_tri_boot_preds <- vcmax25_tri_boot$t  
# Each column contains bootstrapped values (n = 1000) for each row 
# in `trt_df_placeholder`

# Combine bootstrapped predictions with doy and gm.trt descriptorys
full_vcmax25_tri_boot_data <- data.frame(t(vcmax25_tri_boot_preds)) %>% 
  cbind(trt_df_placeholder) %>%
  dplyr::select(doy, gm.trt, X1:X1000) %>%
  mutate(spp = "Tri",
         across(X1:X1000, \(x) exp(x))) # Backtransform estimates to response scale
head(full_vcmax25_tri_boot_data)

# Move bootstrap data to long-format for easy load into climate dataset
vcmax25_tri_boot_data_long <- full_vcmax25_tri_boot_data %>%
  pivot_longer(cols = X1:X1000, 
               names_to = "boot_iter", 
               values_to = "Vcmax25_est")

# View full bootstrap prediction data sheet for Trillium Vcmax25
head(vcmax25_tri_boot_data_long)

# Quick plot - does this look similar to original plot?
ggplot(data = vcmax25_tri_boot_data_long, 
       aes(x = doy, y = Vcmax25_est, fill = gm.trt)) +
  geom_point(shape = 21, alpha = 0.2) + 
  labs(x = "Day of Year", y = expression("V"["cmax25_pred"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")"),
       fill = expression(italic("Alliaria")* " treatment")) +
  theme_bw(base_size = 18)
# Yes!

#####################################################################
# Jmax25 - Tri (does not account for size... yet)
#####################################################################
photo_traits$jmax25[c(296)] <- NA

jmax25_tri <- lmer(log(jmax25) ~ gm.trt * doy + (1|id), 
                          data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(jmax25_tri)
qqnorm(residuals(jmax25_tri))
qqline(residuals(jmax25_tri))
densityPlot(residuals(jmax25_tri))
shapiro.test(residuals(jmax25_tri))
outlierTest(jmax25_tri)

# Bootstrap model predictions to get 1000 Jmax values per doy for
# each of gm.trt == "ambient" and gm.trt == "weeded"
jmax25_tri_boot <- bootMer(
  jmax25_tri,
  predict_function,
  nsim = 1000,
  type = "parametric")

# Extract bootstrapped predictions
jmax25_tri_boot_preds <- jmax25_tri_boot$t  
# Each column contains bootstrapped values (n = 1000) for each row 
# in `trt_df_placeholder`

# Combine bootstrapped predictions with doy and gm.trt descriptorys
full_jmax25_tri_boot_data <- data.frame(t(jmax25_tri_boot_preds)) %>% 
  cbind(trt_df_placeholder) %>%
  dplyr::select(doy, gm.trt, X1:X1000) %>%
  mutate(spp = "Tri",
         across(X1:X1000, \(x) exp(x))) # Backtransform estimates to response scale
head(full_jmax25_tri_boot_data)

# Move bootstrap data to long-format for easy load into climate dataset
jmax25_tri_boot_data_long <- full_jmax25_tri_boot_data %>%
  pivot_longer(cols = X1:X1000, 
               names_to = "boot_iter", 
               values_to = "Jmax25_est")

# View full bootstrap prediction data sheet for Trillium Vcmax25
head(jmax25_tri_boot_data_long)

# Quick plot - does this look similar to original plot?
ggplot(data = jmax25_tri_boot_data_long, 
       aes(x = doy, y = Jmax25_est, fill = gm.trt)) +
  geom_point(shape = 21, alpha = 0.2) + 
  labs(x = "Day of Year", y = expression("J"["max25_pred"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")"),
       fill = expression(italic("Alliaria")* " treatment")) +
  theme_bw(base_size = 18)
# Yes!

#####################################################################
# Ci - Tri (does not account for size... yet)
#####################################################################
ci_tri <- lmer(ci ~ gm.trt * doy + (1|id), 
                    data = subset(photo_traits, spp == "Tri" & ci > 0))

# Check model assumptions
plot(ci_tri)
qqnorm(residuals(ci_tri))
qqline(residuals(ci_tri))
densityPlot(residuals(ci_tri))
shapiro.test(residuals(ci_tri))
outlierTest(ci_tri)

# Bootstrap model predictions to get 1000 Vcmax values per doy for
# each of gm.trt == "ambient" and gm.trt == "weeded"
ci_tri_boot <- bootMer(
  ci_tri,
  predict_function,
  nsim = 1000,
  type = "parametric"
)

# Extract bootstrapped predictions
ci_tri_boot_preds <- ci_tri_boot$t  
# Each column contains bootstrapped values (n = 1000) for each row 
# in `trt_df_placeholder`

# Combine bootstrapped predictions with doy and gm.trt descriptors
full_ci_tri_boot_data <- data.frame(t(ci_tri_boot_preds)) %>% 
  cbind(trt_df_placeholder) %>%
  dplyr::select(doy, gm.trt, X1:X1000) %>%
  mutate(spp = "Tri")
head(full_ci_tri_boot_data)

# Move bootstrap data to long-format for easy load into climate dataset
ci_tri_boot_data_long <- full_ci_tri_boot_data %>%
  pivot_longer(cols = X1:X1000, 
               names_to = "boot_iter", 
               values_to = "Ci_est")

# View full bootstrap prediction data sheet for Trillium Vcmax25
head(ci_tri_boot_data_long)

# Quick plot - does this look similar to original plot?
ggplot(data = ci_tri_boot_data_long, 
       aes(x = doy, y = Ci_est, fill = gm.trt)) +
  geom_point(shape = 21, alpha = 0.2) + 
  labs(x = "Day of Year", y = expression("C"["i_pred"]*" ("*mu*"mol mol"^"-1"*")"),
       fill = expression(italic("Alliaria")* " treatment")) +
  theme_bw(base_size = 18)
# Yes!

#####################################################################
# Vcmax25 - Mai (does not account for size... yet)
#####################################################################
vcmax25_mai <- lmer(log(vcmax25) ~ gm.trt * doy + (1 | id), 
                    data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(vcmax25_mai)
qqnorm(residuals(vcmax25_mai))
qqline(residuals(vcmax25_mai))
densityPlot(residuals(vcmax25_mai))
shapiro.test(residuals(vcmax25_mai))
outlierTest(vcmax25_mai)

# Bootstrap model predictions to get 1000 Vcmax values per doy for
# each of gm.trt == "ambient" and gm.trt == "weeded"
vcmax25_mai_boot <- bootMer(
  vcmax25_mai,
  predict_function,
  nsim = 1000,
  type = "parametric"
)

# Extract bootstrapped predictions
vcmax25_mai_boot_preds <- vcmax25_mai_boot$t  
# Each column contains bootstrapped values (n = 1000) for each row 
# in `trt_df_placeholder`

# Combine bootstrapped predictions with doy and gm.trt descriptorys
full_vcmax25_mai_boot_data <- data.frame(t(vcmax25_mai_boot_preds)) %>% 
  cbind(trt_df_placeholder) %>%
  dplyr::select(doy, gm.trt, X1:X1000) %>%
  mutate(spp = "Mai",
         across(X1:X1000, \(x) exp(x))) # Backtransform estimates to response scale
head(full_vcmax25_mai_boot_data)

# Move bootstrap data to long-format for easy load into climate dataset
vcmax25_mai_boot_data_long <- full_vcmax25_mai_boot_data %>%
  pivot_longer(cols = X1:X1000, 
               names_to = "boot_iter", 
               values_to = "Vcmax25_est")

# View full bootstrap prediction data sheet for Trillium Vcmax25
head(vcmax25_mai_boot_data_long)

# Quick plot - does this look similar to original plot?
ggplot(data = vcmax25_mai_boot_data_long, 
       aes(x = doy, y = Vcmax25_est, fill = gm.trt)) +
  geom_point(shape = 21, alpha = 0.2) + 
  labs(x = "Day of Year", y = expression("V"["cmax25_pred"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")"),
       fill = expression(italic("Alliaria")* " treatment")) +
  theme_bw(base_size = 18)
# Yes!

#####################################################################
# Jmax25 - Mai (does not account for size... yet)
#####################################################################
jmax25_mai <- lmer(log(jmax25) ~ gm.trt * doy + (1|id), 
                   data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(jmax25_mai)
qqnorm(residuals(jmax25_mai))
qqline(residuals(jmax25_mai))
densityPlot(residuals(jmax25_mai))
shapiro.test(residuals(jmax25_mai))
outlierTest(jmax25_mai)

# Bootstrap model predictions to get 1000 Jmax values per doy for
# each of gm.trt == "ambient" and gm.trt == "weeded"
jmax25_mai_boot <- bootMer(
  jmax25_mai,
  predict_function,
  nsim = 1000,
  type = "parametric")

# Extract bootstrapped predictions
jmax25_mai_boot_preds <- jmax25_mai_boot$t  
# Each column contains bootstrapped values (n = 1000) for each row 
# in `trt_df_placeholder`

# Combine bootstrapped predictions with doy and gm.trt descriptors
full_jmax25_mai_boot_data <- data.frame(t(jmax25_mai_boot_preds)) %>% 
  cbind(trt_df_placeholder) %>%
  dplyr::select(doy, gm.trt, X1:X1000) %>%
  mutate(spp = "Mai",
         across(X1:X1000, \(x) exp(x))) # Backtransform estimates to response scale
head(full_jmax25_mai_boot_data)

# Move bootstrap data to long-format for easy load into climate dataset
jmax25_mai_boot_data_long <- full_jmax25_mai_boot_data %>%
  pivot_longer(cols = X1:X1000, 
               names_to = "boot_iter", 
               values_to = "Jmax25_est")

# View full bootstrap prediction data sheet for Maianthemum Jmax25
head(jmax25_mai_boot_data_long)

# Quick plot - does this look similar to original plot?
ggplot(data = jmax25_mai_boot_data_long, 
       aes(x = doy, y = Jmax25_est, fill = gm.trt)) +
  geom_point(shape = 21, alpha = 0.2) + 
  labs(x = "Day of Year", y = expression("J"["max25_pred"]*" ("*mu*"mol m"^"-2"*"s"^"-1"*")"),
       fill = expression(italic("Alliaria")* " treatment")) +
  theme_bw(base_size = 18)
# Yes!

#####################################################################
# Ci - Mai (does not account for size... yet)
#####################################################################
photo_traits$ci[c(367, 392, 790, 791, 820)] <- NA
photo_traits$ci[c(171, 393, 450)] <- NA

ci_mai <- lmer(ci ~ gm.trt * doy + (1|id), 
               data = subset(photo_traits, spp == "Mai" & ci > 0))

# Check model assumptions
plot(ci_mai)
qqnorm(residuals(ci_mai))
qqline(residuals(ci_mai))
densityPlot(residuals(ci_mai))
shapiro.test(residuals(ci_mai))
outlierTest(ci_mai)

# Bootstrap model predictions to get 1000 Vcmax values per doy for
# each of gm.trt == "ambient" and gm.trt == "weeded"
ci_mai_boot <- bootMer(
  ci_mai,
  predict_function,
  nsim = 1000,
  type = "parametric")

# Extract bootstrapped predictions
ci_mai_boot_preds <- ci_mai_boot$t  
# Each column contains bootstrapped values (n = 1000) for each row 
# in `trt_df_placeholder`

# Combine bootstrapped predictions with doy and gm.trt descriptors
full_ci_mai_boot_data <- data.frame(t(ci_mai_boot_preds)) %>% 
  cbind(trt_df_placeholder) %>%
  dplyr::select(doy, gm.trt, X1:X1000) %>%
  mutate(spp = "Mai")
head(full_ci_mai_boot_data)

# Move bootstrap data to long-format for easy load into climate dataset
ci_mai_boot_data_long <- full_ci_mai_boot_data %>%
  pivot_longer(cols = X1:X1000, 
               names_to = "boot_iter", 
               values_to = "Ci_est")

# View full bootstrap prediction data sheet for Trillium Vcmax25
head(ci_mai_boot_data_long)

# Quick plot - does this look similar to original plot?
ggplot(data = ci_mai_boot_data_long, 
       aes(x = doy, y = Ci_est, fill = gm.trt)) +
  geom_point(shape = 21, alpha = 0.2) + 
  labs(x = "Day of Year", y = expression("C"["i_pred"]*" ("*mu*"mol mol"^"-1"*")"),
       fill = expression(italic("Alliaria")* " treatment")) +
  theme_bw(base_size = 18)
# Yes!

#####################################################################
# Merge bootstrapped dataset for each
#####################################################################
tri_boot_model_params <- vcmax25_tri_boot_data_long %>%
  full_join(jmax25_tri_boot_data_long) %>%
  full_join(ci_tri_boot_data_long)

mai_boot_model_params <- vcmax25_mai_boot_data_long %>%
  full_join(jmax25_mai_boot_data_long) %>%
  full_join(ci_mai_boot_data_long)


boot_model_params <- tri_boot_model_params %>%
  full_join(mai_boot_model_params)

write.csv(boot_model_params, 
          "../data/c_budget/TT24_budget_params_without_size.csv",
          row.names = F)



