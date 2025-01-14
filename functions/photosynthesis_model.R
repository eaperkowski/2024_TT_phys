# Simple photosynthesis model derived from Farquhar et al. (1980)
# with some additional helper functions for standardizing
# Vcmax/Jmax estimates to a given air temperature, calculating
# Michaelis-Menton andd gammaStar coefficients. Equations based
# on those explained in Smith et al. (2019) and Stocker et al. (2020).
# Original helper functions can be found at following URL:
# https://github.com/SmithEcophysLab/optimal_vcmax_R

###############################################################################
# Helper functions (should be loaded into working environment before 
# photosynthesis model is run)
###############################################################################

# Temp correction for Vcmax, Jmax
temp_standardize_vcmax = function(tleaf, tmean, tref) {
  
  temp <- tleaf + 273.15
  Ha <- 71513
  Hd <- 200000
  adelS <- 668.39
  bdelS <- -1.07
  tmeanK <- tmean + 273.15
  trefK <- tref + 273.15
  R <- 8.314
  kbeg <- exp(Ha*(temp - trefK) / (trefK * R * temp))
  kend <- ((1 + exp((trefK * (adelS + bdelS * tmean) - Hd) / 
                      (trefK * R))) / (1 + exp((temp * (adelS + bdelS * tmean) - Hd) / (temp * R))))
  kbeg * kend
  
}

# tLeaf = leaf temperature, tmean = growing season temp, tref = 
temp_standardize_jmax = function(tleaf, tmean, tref) {

  temp <- tleaf + 273.15
  Ha <- 49884
  Hd <- 200000
  adelS <- 659.7
  bdelS <- -0.75
  tmeanK <- tmean + 273.15
  trefK <- tref + 273.15
  R <- 8.314
  kbeg <- exp(Ha * (temp - trefK) / (trefK * R * temp))
  kend <- ((1 + exp((trefK * (adelS + bdelS * tmean) - Hd) / (trefK * R))) / 
             (1 + exp((temp * (adelS + bdelS * tmean) - Hd) / (temp * R))))
  kbeg * kend
    
  }


# Calculate Km (needed for mc to calc Ac)
calc_km_pa = function(temp) {

  R = 8.314        
  O2 = 2.09476e5      
  Kc25 = 41.03 
  Ko25 = 28210 
  Hkc = 79430  
  Hko = 36380 
  
  temp_k = 273.15 + temp
  
  Kc_pa = Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  Ko_pa = Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  O2_pa = O2 * (1e-6)
  
  Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
  
  Km_pa 
}

# Calculate gammastar (needed for mc and m to calc Ac and Aj)
calc_gammastar_pa = function(temp) {
  
  #gammastar25 = 42.75  # ppm
  gammastar25 = 4.332 # Pa
  Hgm = 37830 # J mol-1
  R = 8.314 # J K-1 mol-1
  O2 = 2.09476e5 # ppm
  O2_0 = O2 * 1e-6
  
  temp_k = 273.15 + temp
  
  gStar_pa = gammastar25 * exp((Hgm / R) * ((1 / 298.15) - (1 / temp_k)))
  
  gStar_pa
  
}

## Full photosynthesis model
photosynthesis <- function(temp_c = 25,
                           temp_c_rolling = 25, # rolling 10 day mean temp, for Vcmax/Jmax temp correction
                           ci = 400,
                           par = 800,
                           patm = 101325,
                           q0 = 0.257, # Quantum yield of electron transport, from Smith et al. (2019)
                           theta = 0.85, # Distribution of light intensity relative to distribution of photosynthetic capacity,
                                         # from Smith et al. (2019)
                           vcmax25 = 100, 
                           jmax25 = 200) {
  
  #  Michaelis-Menton coefs and gammastar from temp eqs. from Bernacchi et al. (2001)
  km <- calc_km_pa(temp_c)
  gammastar <- calc_gammastar_pa(temp_c)
  
  # Calculate Vcmax and Jmax at temp_c using temp corrections from Kattge & Knorr (2007)
  vcmax = vcmax25 * temp_standardize_vcmax(tleaf = temp_c, tmean = temp_c_rolling, tref = 25)
  jmax = jmax25 * temp_standardize_jmax(tleaf = temp_c, tmean = temp_c_rolling, tref = 25)
  
  # Convert Ci from umol/mol to Pa
  ci_pa <- ci * 1e-6 * patm
  
  # Ac
  mc <- ((ci_pa - gammastar) / (ci_pa + km))
  ac <- vcmax * mc
  
  # Aj
  m <- ((ci_pa - gammastar) / ( ci_pa + (2 * gammastar)))
  aj <- (m / 4) * (q0 * par + jmax - sqrt((q0 * par + jmax)^2 - (4 * theta * q0 * par * jmax))) / (2 * theta)

  # Anet is the minimum of Ac and Aj
  a <- pmin(ac, aj)
  rd <- vcmax * 0.015
  anet <- a - (vcmax * 0.015)
  
  results <- data.frame(temp_c = temp_c, temp_c_rolling = temp_c_rolling, 
                        ci = ci, ci_pa = ci_pa, par = par, q0 = q0, theta = theta, 
                        vcmax = vcmax, jmax = jmax, km = km, gammastar = gammastar, 
                        mc = mc, m = m, ac = ac, aj = aj, a = a, rd = rd, 
                        anet = anet)
  
  return(results)
  
}
