#####################################################################
# Load libraries, photosynthesis model, plot aesthetics, and datasets
#####################################################################

# Libraries
library(tidyverse)
library(lubridate)
library(lme4)
library(emmeans)
library(car)
library(ggpubr)

# Load photosynthesis model
source("../functions/photosynthesis_model.R")

# Some plot aesthetics
gm.colors <- c("#00B2BE", "#F1B700")
facet_lab_tri <- "Trillium spp."
names(facet_lab_tri) <- "Tri"
facet_lab_mai <- "M. racemosum"
names(facet_lab_mai) <- "Mai"

# Read in weather datasets
## Hourly climate data
weather_hourly <- read.csv("../data/TT24_hourly_weather_data.csv") %>%
  filter(doy > 99 & doy < 250) %>%
  mutate(atm_pressure_pa = atm_pressure_kpa * 1000,
         hour = hour(date)) %>%
  dplyr::select(date_time = date, 
                doy,
                hour,
                par = par_umol_m2_s,
                temp_c = air_temperature_c,
                patm = atm_pressure_pa) %>%
  separate(date_time, c("date", "time"), sep = " ", remove = FALSE)

## Daily climate data to get rolling means for temps
weather_daily <- read.csv("../data/TT24_daily_weather_summary.csv") %>%
  filter(doy > 99 & doy < 250) %>%
  dplyr::select(doy, tavg10)

# Merge rolling mean with hourly climate
weather <- weather_hourly %>%
  left_join(weather_daily, by = "doy")
head(weather)

# Read in compiled trait dataset
photo_traits <- read.csv("../data/TT24_photo_traits.csv") %>%
  mutate(doy_photo = doy) %>%
  group_by(id) %>%
  mutate(doy_max = max(doy)) %>%
  select(-doy)
head(photo_traits)

# How many measurements per ID?
n_measurements <- photo_traits %>%
  group_by(id, spp, plot, subplot, gm.trt) %>%
  summarize(n_measurements = length(id)) %>%
  arrange(n_measurements)

# 
id <- unique(photo_traits$id)

photo_traits_doy <- photo_traits %>%
  select(id, doy_photo, doy_max)

#####################################################################
# Create data frame that adds id of all individuals in hourly weather 
# dataset, then adds doy of each photo measurement and the final 
# photo measurement. Then, add photo_traits for each individual, 
# filter dataset to include only the most recent photo measurement 
# per doy of weather data. Final prep for spitting data into 
# photosynthesis model
#####################################################################
weather_photo_merged <- weather %>%
  crossing(id = id) %>%
  left_join(photo_traits_doy, by = "id", relationship = "many-to-many") %>%
  left_join(photo_traits, by = c("id", "doy_photo", "doy_max")) %>%
  mutate(time_diff = abs(doy - doy_photo)) %>%
    group_by(doy, id) %>%
    filter(time_diff == min(time_diff)) %>% # Keep the closest measurement
    ungroup() %>%
    dplyr::select(-time_diff, date_time, weather_date = date.x, time:patm, tavg10_weather = tavg10.x, 
           machine, photo_date = date.y, doy_photo, doy_max, id,spp:Tleaf, tavg10_vcmax25 = tavg10.y,
           init_leaf_length_cm:init_leaf_area, -(vp10:vpd10), -anet, -vcmax, -jmax, -rd) 
head(weather_photo_merged)

#####################################################################
# Run photosynthesis model!
#####################################################################
model_results <- photosynthesis(temp_c = weather_photo_merged$temp_c,
                                temp_c_rolling = weather_photo_merged$tavg10_weather,
                                ci = weather_photo_merged$ci,
                                par = weather_photo_merged$par,
                                patm = weather_photo_merged$patm,
                                vcmax25 = weather_photo_merged$vcmax25,
                                jmax25 = weather_photo_merged$jmax25)
head(model_results)

# Merge photosynthesis model results with `weather_photo_merged`
hourly_model_results <- cbind(weather_photo_merged, model_results)
head(hourly_model_results)

# Remove repeat cols (Vcmax, Jmax, Ci)
hourly_model_results <- hourly_model_results[, -c(35:37, 39)]

#####################################################################
# Convert A, Rd, and Anet from umol/m2/s to g C/ms/s, then to
# g C / day
#####################################################################
hourly_model_results <- hourly_model_results %>%
  mutate(hourly_C_assim = a * 3600 / 1e6 * 12,
         hourly_C_resp = rd * 3600 / 1e6 * 12,
         hourly_netC_assim = anet * 3600 / 1e6 * 12) %>%
  mutate(hourly_C_assim_tla = hourly_C_assim * (total_leaf_area_cm2 / 10000),
         hourly_C_resp_tla = hourly_C_resp * (total_leaf_area_cm2 / 10000),
         hourly_netC_assim_tla = hourly_netC_assim * (total_leaf_area_cm2 / 10000))

#####################################################################
# Calculate total daily C estimates
#####################################################################
daily_model_results <- hourly_model_results %>%
  group_by(doy, id, spp, plot, subplot, gm.trt) %>%
  summarize(daily_C_assim = sum(hourly_C_assim, na.rm = TRUE),
            daily_C_resp = sum(hourly_C_resp, na.rm = TRUE),
            daily_netC_assim = sum(hourly_netC_assim, na.rm = TRUE),
            
            daily_C_assim_tla = sum(hourly_C_assim_tla, na.rm = TRUE),
            daily_C_resp_tla = sum(hourly_C_resp_tla, na.rm = TRUE),
            daily_netC_assim_tla = sum(hourly_netC_assim_tla, na.rm = TRUE)) %>%
  ungroup(doy) %>%
  mutate(cumul_C_assim = cumsum(daily_C_assim),
         cumul_C_resp = cumsum(daily_C_resp),
         cumul_netC_assim = cumsum(daily_netC_assim),
         
         cumul_C_assim_tla = cumsum(daily_C_assim_tla),
         cumul_C_resp_tla = cumsum(daily_C_resp_tla),
         cumul_netC_assim_tla = cumsum(daily_netC_assim_tla),
         
         
         cumul_C_assim_fract = cumul_C_assim / max(cumul_C_assim),
         cumul_C_resp_fract = cumul_C_resp / max(cumul_C_assim),
         cumul_netC_assim_fract = cumul_netC_assim / max(cumul_netC_assim)) %>%
  left_join(n_measurements)
head(daily_model_results)

# Check that cumsum() worked
filter(daily_model_results, id == 9412)

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
  geom_rect(aes(xmin = 125, xmax = 170, ymin = -Inf, ymax = Inf),
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
  facet_grid(~spp, labeller = labeller(spp = facet_lab_tri)) +
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


