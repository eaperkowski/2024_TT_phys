calc_chi_c3 <- function(leaf.d13c = NA, air = NA, year = NA) {
  
  ## Derive d13C in air as a function of calendar year, 
  ## as done in Feng et al. (1999), with option to directly modify
  ## air fractionation value
  if(is.na(year) & is.na(air)) {air =-8}
  if(is.na(year) & !is.na(air)) {air = air}
  if(!is.na(year) & is.na(air)) {air = -6.429 - 0.006*exp(0.0217*(year-1740))}
  
  # Global constants
  a = 4.4
  b = 27
  
  # Leaf carbon discrimination relative to air
  delta <- (air - leaf.d13c) / (1 + leaf.d13c * 0.001)
  chi <- (delta - a) / (b - a)
  
  return(data.frame(delta = delta, chi = chi))
  
}

calc_chi_c4 <- function(leaf.d13c = NA, air = NA, year = NA) {

  ## Derive d13C in air as a function of calendar year, 
  ## as done in Feng et al. (1999), with option to directly modify
  ## air fractionation value
  if(is.na(year) & is.na(air)) {air =-8}
  if(is.na(year) & !is.na(air)) {air = air}
  if(!is.na(year) & is.na(air)) {air = -6.429 - 0.006*exp(0.0217*(year-1740))}
  
  # Global constants
  a = 4.4
  c = -5.7
  d = 30
  phi = 0.40
  
  b.c4 = -5.7 + (d * phi)
  
  # Leaf carbon discrimination relative to air
  delta <- (air - leaf.d13c) / (1 + leaf.d13c * 0.001)
  chi <- (delta - a) / (b.c4 - a)
  
  return(data.frame(delta = delta, chi = chi))
  
}



# Test fxn when year is not included (air.d13C should = -8)
# calc_chi_c3(leaf.d13c = -30)
# calc_chi_c4(leaf.d13c = -13)
# 
# # Test fxn when year is included (air.d13C should be different than -8)
#calc_chi_c3(leaf.d13c = -31, year = 2020)
#calc_chi_c4(leaf.d13c = -13, year = 2020)

