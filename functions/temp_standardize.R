# temp_standardize (estimate, estimate.type, standard.to, tLeaf, tGrow):
# Standardizes Vcmax, Jmax, and Rd estimates to single temperature. Vcmax and
# Jmax are temperature standardized using a modified Arrhenius function as 
# explained in Kattge & Knorr (2007), while Rd is temperature standardized using
# a log polynomial from O'Sullivan (2013) and parameters set using Heskel et al.
# (2016). Function is capable of standardizing to any temperature.
# 
# Created by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Last edit: April 22, 2024
# Last edit by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Arguments:
#    - estimate        =  Vcmax, Jmax, or Rd estimate
#    - estimate.type   =  Notes whether estimate type is a Vcmax, Jmax, or Rd
#                         estimate
#    - pft             =  plant functional type for log polynomial parameters
#                         from Heskel et al. (2016). Use "GM" if unsure. ONLY
#                         NEEDED FOR DARK RESPIRATION
#                           - "BlDcTmp" -> broad leaf deciduous temperate
#                           - "BlDcTrp" -> broad leaf deciduous tropical
#                           - "BlEvTmp" -> broad leaf evergreen temperate
#                           - "BlEvTrp" -> broad leaf evergreen tropical
#                           - "C3H"     -> C3 herbaceous
#                           - "NlEv"    -> needle-leaf evergreen
#                           - "SEv"     -> needle-leaf shrub
#                           - "GM"      -> Global mean
#    - standard.to     =  temperature to standardize rate estimate to. In deg C
#    - tLeaf           =  average leaf temp. during estimate. Must be in degC
#    - tGrow           =  avg. growing season temp. or avg. temp across entire 
#                         experiment
#
# Returns:
# Vector of each Vcmax, Jmax, or Rd estimate standardized to the designated
# temperature
temp_standardize <- function(estimate,
                            estimate.type = c("Vcmax", "Jmax", "Rd"),
                            standard.to,
                            pft,
                            tLeaf,
                            tGrow) {

  ## delta S constants as per Kattge & Knorr (2007)
  a_vcmax <- 668.39
  b_vcmax <- -1.07
  a_jmax <- 659.70
  b_jmax <- -0.75
  
  ## Calculate delta S for Vcmax
  S_vcmax <- a_vcmax + b_vcmax * tGrow
  
  ## Calculate delta S for Jmax
  S_jmax <- a_jmax + b_jmax * tGrow
  
  ## Constants to standardize Vcmax and Jmax Rd estimates to 25 degC
  tK <- tLeaf + 273.15
  tO <- standard.to + 273.15
  Ha_vcmax <- 71513
  Ha_jmax <- 49884
  Hd <- 200000
  R <- 8.314
  
  if (estimate.type == "Vcmax") {
  multOneVcmax <- exp((Ha_vcmax * (tK - tO)) / (R * tK * tO))
  multTwoVcmax <- (1 + exp((tO * S_vcmax - Hd)/(R * tO))) / (1 + exp((tK * S_vcmax - Hd)/(R * tK)))
  
  multipliersVcmax <- multOneVcmax * multTwoVcmax
  VcmaxStandard <- estimate / multipliersVcmax
  
  return(VcmaxStandard)
  
  }
  
  if(estimate.type == "Jmax") {
  
  multOneJmax <- exp((Ha_jmax * (tK - tO)) /(R * tK * tO))
  multTwoJmax <- (1 + exp((tO * S_jmax - Hd)/(R * tO))) / (1 + exp((tK * S_jmax - Hd)/(R * tK)))
  
  multipliersJmax <- multOneJmax * multTwoJmax
  JmaxStandard <- estimate / multipliersJmax
  
  return(JmaxStandard)
  
  }
  
  if(estimate.type == "Rd") {
    
    if(pft == "BlDcTmp") {
      a_rd = -2.2264
      b_rd = 0.0993
      c_rd = -0.0005
    }
    else if(pft == "BlDcTrp") {
      a_rd = -2.727
      b_rd = 0.1125
      c_rd = 0.00058
    }
    
    else if(pft == "BlEvTmp") {
      a_rd = -1.8106
      b_rd = 0.0896
      c_rd = -0.00021
    }
    
    else if(pft == "BlEvTrp") {
      a_rd = -2.6105
      b_rd = 0.1022
      c_rd = -0.00052
    }
    
    else if(pft == "C3H") {
      a_rd = -1.7507
      b_rd = 0.1271
      c_rd = -0.00110
    }
    
    else if(pft == "NlEv") {
      a_rd = -2.0464
      b_rd = 0.1125
      c_rd = -0.00063
    }
    
    else if(pft == "SEv") {
      a_rd = -1.815
      b_rd = 0.0971
      c_rd = -0.00047
    }
    
    else if(pft == "GM") {
      a_rd = -2.2276
      b_rd = 0.01012
      c_rd = -0.0005
    }
    else{
      warning("No recognized plant functional type. Value calculated using global mean parameters")
      
      a_rd = -2.2276
      b_rd = 0.01012
      c_rd = -0.0005
    }

    ## Equation from Heskel et al. (2016) and O'Sullivan et al. (2013)
    rd.25 = estimate * exp(b_rd * 
                             (standard.to - tLeaf) + c_rd * (
                               (standard.to)^2 - (tLeaf)^2))
    
    return(rd.25)
    
  }
}


