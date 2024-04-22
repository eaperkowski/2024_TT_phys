## stomatal_limitation(A_net, Vcmax, elevation, temp):
# Created by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Last edit: April 22, 2024
# Last edit by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Calculates Γ* and the Michaelis-Menton for CO2 in Rubisco following 
# leaf temperature equations explained in Bernacchi et al. (2001). Then,
# calculates stomatal limitation from equations explained in Sharkey &
# Farquhar (1982)
#
#
# Arguments:
#    - A_net        =  Net photosynthesis (A) at 400ppm CO2
#    - Vcmax        =  Unstandardized Vcmax measurement (because A_net is 
#                      not standardized) 
#    - temp         =  leaf temperature by which A_net and Vcmax were measured.
#                      Required to calculate Γ* and Km
#    - Rd.meas      =  Default FALSE. If false, calculates as function of Vcmax
#    - Rd           =  measured dark respiration. Rd.meas should be set to TRUE
#                      if Rd data are included
#    - K            =  Boolean noting whether temperature is in K or not. Default
#                      is FALSE (noting temperature is not in K)
#
# Returns:
#    - Data frame containing Michaelis-Menton constants (Kc, Ko, Km), gammastar,
#      A_mod, and stomatal limitation value (stomLim).
# 
# Note: Stomatal limitation should be a value less than or equal to one
stomatal_limitation <- function(A_net, Vcmax, leaf.temp, 
                                Rd.meas = FALSE, Rd, temp = c("C", "K")) {
  
  ## Global constants
  R = 8.314 # universal gas constant (J mol^-1 K^-1)
  Oi = 210 # leaf intercellular O2 concentration (mmol mol^-1)
  
  ## Leaf temperature correction
  leaf.temp <- ifelse(temp == "C", 
                      leaf.temp + 273.15, 
                      ifelse(temp == "K",
                             leaf.temp,
                             NA))

  ## Calculate Kc
  Kc <- 404.9*exp((79430*(leaf.temp - 298))/(298*R*leaf.temp))

  ## Calculate Ko
  Ko <- 278.4*exp((36380*(leaf.temp - 298))/(298*R*leaf.temp))
  
  ## Calculate Km
  Km <- Kc * (1 + (Oi / Ko))
  
  ## Calculate gamma.star
  gammastar <- 42.75*exp((37830*(leaf.temp - 298))/(298*R*leaf.temp))
   
   # Add contingency over whether Rd was measured or not
   if(Rd.meas == FALSE) {
     Rd <- 0.015 * Vcmax
   }
   
   if(Rd.meas == TRUE) {
     Rd <- Rd
   }

   # Calculate A_mod
   A_mod = Vcmax * ((400 - gammastar) / (400 + Km)) - Rd
   
   # Determine stomatal limitation
   l <- 1 - (A_net / A_mod)
   
   return(data.frame(Kc, Ko, Km, gammastar, l))
}


