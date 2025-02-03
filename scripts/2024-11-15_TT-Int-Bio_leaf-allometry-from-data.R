## Trillium Trail Integrative Bio
## Leaf allometry equations
## Code written and annotated by JM 11-2024

# Set working directory
setwd("~/")

# Read in data of plant size and leaf area measurements for Trillium spp. and
#  Maianthemum racemosum.
allo_dat <- read.csv("../data/TT24_nonfocal_data_with_fluorescence.csv")

# I'm gonna plot some stuff using ggplot.
require(ggplot2)

## =============================================================================
## (1) Total leaf area of stem from longest leaf length for Trillium spp.
## =============================================================================

# A previous estimate of Trillium total leaf area (developed by Robert M) was 
#  total_leaf_area = leaf_length * leaf_length * (1 + cos(pi/3) + sin(pi/3)),
#  which estimates the area of 3 equilateral triangles where leaf_length is the 
#  height of each triangle.

# Start by visualizing relationship between longest leaf length and total leaf area.
tri_plot <- ggplot() +
  geom_point(data = subset(allo_dat, spp == "tri"),
             aes(x = leaf_length_cm, y = total_leaf_area_cm2)) +
  theme_bw()

tri_plot # Probably looking at a log-log relationship here

# Fit log-linear model with intercept set to zero.
tri_allo <- lm(log(total_leaf_area_cm2) ~ log(leaf_length_cm) - 1,
               data = subset(allo_dat, spp == "tri"))
summary(tri_allo)

# From this model, we can use total_leaf_area = exp(log(leaf_length)*2.2).

# Plot this to the data to visually check fit, and add in the previous allometric
#  equation (based on leaves as equilateral triangles) for comparison
tri_plot + 
  #log-linear with intercept set to zero (solid line)
  geom_line(aes(x = seq(7,20), y = exp(2.2*log(seq(7,20))))) + 
  #previous allometric equation based on area of 3 equilateral triangles (dashed line)
  geom_line(aes(x = seq(7,20), y = seq(7,20)*seq(7,20)*(1 + cos(pi/3) + sin(pi/3))),
            linetype = "dashed")

# The assumption of leaves as 3 equilateral triangles overestimates leaf area,
#  a bit more so for larger leaves, but the shape of the relationship is decent.


## =============================================================================
## (2) Total leaf area of stem from EITHER (a) longest leaf length OR (b) stem
##     height in Maianthemum.
## =============================================================================

# (a) Estimating from longest leaf length. -----
# A previous estimate of total leaf area of a Maianthemum stem from longest leaf 
#  length is total_leaf_area = exp(2.4*log(leaf_length)), where leaf areas were
#  estimated as rectangles from leaf length and width.

# Start by visualizing the relationship between longest leaf length and total 
#  leaf area.
mai_plot1 <- ggplot() +
  geom_point(data = subset(allo_dat, spp == "mai"),
             aes(x = leaf_length_cm, y = total_leaf_area_cm2)) +
  theme_bw()

mai_plot1 # Looking like a log-log (or possibly quadratic?) relationship

# Fit log linear model with intercept set to zero.
mai_allo1 <- lm(log(total_leaf_area_cm2) ~ log(leaf_length_cm) - 1,
                data = subset(allo_dat, spp == "mai"))
summary(mai_allo1)

# This model suggests that total_leaf_area = exp(log(leaf_length)*2.2). (Same
#  as with Trillium -- interesting.)

# Plot this to the data to visually check fit.
mai_plot1 +
  geom_line(aes(x = seq(0,25), y = exp(log(seq(0,25))*2.2))) 

# This underestimates total leaf area for larger leaves.

# Try quadratic instead?
mai_allo2 <- lm(total_leaf_area_cm2 ~ leaf_length_cm + I(leaf_length_cm^2) - 1,
                data = subset(allo_dat, spp == "mai"))
summary(mai_allo2)

# This model suggests total_leaf_area = -20.2*leaf_length + 3.2*leaf_length^2. 
#  Because of the negative linear coefficient, the predictions go negative in 
#  the range of approx. 0-7.5 (which is not ideal):
mai_plot1 +
  geom_line(aes(x = seq(0,25), y = -20.2*seq(0,25) + 3.2*seq(0,25)^2))

# Trying a log-linear model again, but without fixing the intercept to zero,
#  which will mean that there can be a non-zero total leaf area when longest leaf
#  length = 0, but it will be negligibly small, and positive (in terms of biological
#  reality, this seems preferable to a rather large range where leaf area is
#  predicted to be negative)
mai_allo3 <- lm(log(total_leaf_area_cm2) ~ log(leaf_length_cm),
                data = subset(allo_dat, spp == "mai"))
summary(mai_allo3)

# This model suggests that total_leaf_area = exp(log(leaf_length)*2.7 - 1.4)

# Plot predictions from all three models, with data, for comparison:
mai_plot1 +
  #log-linear with intercept set to zero (solid line)
  geom_line(aes(x = seq(0,25), y = exp(log(seq(0,25))*2.2))) +
  #quadratic with intercept set to zero (dashed line)
  geom_line(aes(x = seq(0,25), y = -20.2*seq(0,25) + 3.2*seq(0,25)^2),
            linetype = "dashed") +
  #log-linear with non-zero intercept (dotted line)
  geom_line(aes(x = seq(0,25), y = exp(log(seq(0,25))*2.7 - 1.4)),
            linetype = "dotted")
  

# (b) Estimating from stem height. -----
# A previous estimate of total leaf area (developed by Robert M) of a Maianthemum 
#  stem from stem height is total_leaf_area = 1.456*leaf_length + 0.149*leaf_length^2.

# Start by visualizing the relationship between stem height and total leaf area.
mai_plot2 <- ggplot() +
  geom_point(data = subset(allo_dat, spp == "mai"),
             aes(x = stem_length_cm, y = total_leaf_area_cm2)) +
  theme_bw()
mai_plot2

# First fitting log-linear model with intercept set to zero.
mai_allo1 <- lm(log(total_leaf_area_cm2) ~ log(stem_length_cm) - 1,
                data = subset(allo_dat, spp == "mai"))
summary(mai_allo1)

# This model suggests that total_leaf_area = exp(log(stem_height)*1.6)

# Plot this to the data to visually check fit:
mai_plot2 +
  geom_line(aes(x = seq(10,80), y = exp(log(seq(10,80))*1.6)))

# Looks pretty good.

# Fit quadratic model as well, to align with the previous allometric formula.
mai_allo2 <- lm(total_leaf_area_cm2 ~ stem_length_cm + I(stem_length_cm^2) - 1,
                data = subset(allo_dat, spp == "mai"))
summary(mai_allo2)

# This model suggests that total_leaf_area = 5.0*stem_height + 0.11*stem_height^2.

# Plot predictions from the above two models and from the previous allometric 
#  formula, with the data, for comparison.
mai_plot2 +
  #log-linear with intercept set to zero (solid line)
  geom_line(aes(x = seq(10,80), y = exp(log(seq(10,80))*1.6))) +
  #quadratic with intercept set to zero (dashed line)
  geom_line(aes(x = seq(10,80), y = seq(10,80)*5 + 0.11*seq(10,80)^2),
              linetype = "dashed") +
  #previous allometric equation, assuming quadratic relationship (dotted line)
  geom_line(aes(x = seq(10,80), y = seq(10,80)*1.456 + 0.149*seq(10,80)^2),
            linetype = "dotted")

# The fitted log-linear and quadratic equations perform similarly.

