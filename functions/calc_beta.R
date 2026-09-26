# calc_beta(chi = NA, temp = NA, vpd = NA, ca = 420, z = NA):
#
# Calculates beta following Eq. 3 in Lavergne et al. (2020) and originally
# derived in Prentice et al. (2014). Beta is a parameter in photosynthetic 
# least-cost theory frameworks that predicts resource use demand. It is
# a unitless ratio that can also be derived as the ratio of the unit cost of 
# Rubisco carboxylation maintenance (b) to the unit cost of transpirational 
# maintenance (a). This function includes functions that can be accessed at the
# Smith Plantecophys Lab GitHub repository titled "optimal_vcmax_R". Link to the
# repository is included below.
#
# Function inputs:
#   - chi       = Carbon-13 isotope derived estimate of Ci:Ca. Must be between 
#                 0 and 1. Code will otherwise not run.
#   - temp      = acclimated temperature (degrees Celsius)
#   - vpd       = acclimated vapor pressure deficit (Pa)
#   - ca        = atmospheric CO2 (ppm)
#   - z         = elevation (m)
#
# Function outputs:
#   - beta      = resource demand term, calculated per Lavergne et al. (2020)
#
# Repository link: https://github.com/SmithEcophysLab/optimal_vcmax_R
#
# References:
# Lavergne A., Sandoval D., Hare V.J., Graven H., Prentice I.C. 2020. Impacts of 
# soil water stress on the acclimated stomatal limitation of photosynthesis: 
# insights from stable carbon isotope data. Global Change Biology 26: 7158-7172.
#
# Prentice I.C., Dong N., Gleason S.M., Maire V. & Wright I.J. (2014) Balancing 
# the costs of carbon gain and water transport: testing a new theoretical 
# framework for plant functional ecology. Ecology Letters 17, 82–91.
calc_beta <- function(chi = NA, temp = NA, vpd = NA, ca = 420, z = 0) {

  # Calculate atmospheric pressure (Pa) from elevation (m)
  calc_patm = function(z) {
    
    kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
    kTo = 298.15   # base temperature, K (Prentice, unpublished)
    kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
    kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    
    patm = kPo*(1.0 - kL * z / kTo)**(kG * kMa / (kR * kL))
    
    patm
  }
  
  # Calculate viscosity of water
  calc_viscosity_h2o = function(temp, z){ #temp in °C and z in m
    
    tk_ast  = 647.096    # Kelvin
    rho_ast = 322.0      # kg/m**3
    mu_ast  = 1e-6       # Pa s
    
    patm = calc_patm(z)
    
    rho = calc_density_h2o(temp, z) # density of water (kg m-3)
    
    # Calculate dimensionless parameters:
    tbar = (temp + 273.15)/tk_ast
    tbarx = tbar**(0.5)
    tbar2 = tbar**2
    tbar3 = tbar**3
    rbar = rho/rho_ast
    
    # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
    mu0 = 1e2*tbarx/mu0
    
    # Create Table 3, Huber et al. (2009):
    
    h_array1 = c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)
    h_array2 = c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573)
    h_array3 = c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0)
    h_array4 = c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0)
    h_array5 = c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0)
    h_array6 = c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0)
    h_array7 = c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264)
    
    h_array = rbind(h_array1, h_array2, h_array3, h_array4, h_array5, h_array6, h_array7)
    
    # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    
    mu1 = 0.0
    ctbar = (1.0/tbar) - 1.0
    for (i in 1:6){
      
      coef1 = ctbar**(i-1)
      coef2 = 0.0
      
      for (j in 1:7){
        
        coef2 = coef2 + h_array[j,i] * (rbar - 1.0)**(j-1)
        
      }
      
      mu1 = mu1 + coef1 * coef2 
      
    }
    
    mu1 = exp( rbar * mu1 )
    
    # Calculate mu_bar (Eq. 2, Huber et al., 2009)
    #   assumes mu2 = 1
    mu_bar = mu0 * mu1
    
    # Calculate mu (Eq. 1, Huber et al., 2009)
    viscosity_h2o = mu_bar * mu_ast    # Pa s
    viscosity_h2o
    
  }
  
  ## Calculate the density of h2o (kg m-3) given tmeperature and pressure
  calc_density_h2o = function(temp, z){ # temp in °C and z in m
    
    patm = calc_patm(z)
    
    # Calculate lambda, (bar cm**3)/g:
    my_lambda = 1788.316 + 1.55053*temp + -0.4695911*temp*temp + (3.096363e-3)*temp*temp*temp + -(7.341182e-6)*temp*temp*temp*temp
    
    # Calculate po, bar
    po = 5918.499 + 58.05267*temp + -1.1253317*temp*temp + (6.6123869e-3)*temp*temp*temp + -(1.4661625e-5)*temp*temp*temp*temp
    
    # Calculate vinf, cm**3/g
    vinf = 0.6980547 + -(7.435626e-4)*temp + (3.704258e-5)*temp*temp + -(6.315724e-7)*temp*temp*temp + (9.829576e-9)*temp*temp*temp*temp + -(1.197269e-10)*temp*temp*temp*temp*temp + (1.005461e-12)*temp*temp*temp*temp*temp*temp + -(5.437898e-15)*temp*temp*temp*temp*temp*temp*temp + (1.69946e-17)*temp*temp*temp*temp*temp*temp*temp*temp + -(2.295063e-20)*temp*temp*temp*temp*temp*temp*temp*temp*temp
    
    # Convert pressure to bars (1 bar = 100000 Pa)
    pbar = (1e-5)*patm
    
    # Calculate the specific volume (cm**3 g**-1):
    vau = vinf + my_lambda/(po + pbar)
    
    # Convert to density (g cm**-3) -> 1000 g/kg; 1000000 cm**3/m**3 -> kg/m**3:
    density_h2o = (1e3/vau)
    
    density_h2o
    
  }
  
  ## Calculate nstar (unitless relative viscosity of h2o at temperature relative to 25°C)
  calc_nstar = function(temp, z){ # temp in °C and z in m
    
    patm = calc_patm(z)
    
    # viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
    ns      = calc_viscosity_h2o( temp, z )  # Pa s 
    ns25    = calc_viscosity_h2o( 25, z )  # Pa s 
    nstar = ns / ns25                       # (unitless)
    
    nstar
    
  }

  # Calculate gammastar (Pa)
  calc_gammastar_pa = function(temp, z) {
    
    patm = calc_patm(z)
    rat = calc_patm(z) / calc_patm(0)
    
    #gammastar25 = 42.75  # ppm
    gammastar25 = 4.332 * rat  # Pa
    Hgm=37830 # J mol-1
    R = 8.314        # J K-1 mol-1
    O2 = 2.09476e5 # ppm
    O2_0 = O2 * 1e-6 * calc_patm(0)
    O2_z = O2 * 1e-6 * calc_patm(z)
    
    temp_k = 273.15+ temp
    
    gStar_pa = gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
    
    gStar_pa
    
  }

  # Calculate the Michaelis-Menton coefficient (Pa) for Rubisco from temperature
  calc_km_pa = function(temp, z) {
    
    patm = calc_patm(z) 
    rat = patm / calc_patm(0)
    
    R = 8.314        
    O2 = 2.09476e5      
    Kc25 = 41.03 * rat 
    Ko25 = 28210 * rat 
    Hkc = 79430  
    Hko = 36380 
    
    temp_k = 273.15 + temp
    
    Kc_pa =Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
    Ko_pa =Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
    
    O2_pa = O2 * (1e-6) * patm 
    
    Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
    
    Km_pa 
    
  }
  
  # Determine CO2 concentration, nstar, gammaStar, and Km given fxns above
  patm = calc_patm(z)
  
  Ca = ca * 1e-6 * patm
  
  nstar = calc_nstar(temp, z)
  gammaStar = calc_gammastar_pa(temp, z)
  K = calc_km_pa(temp, z)
  
  beta = 1.6 * nstar * vpd * ((chi - (gammaStar / ca))^2 / ((1 - chi)^2 * (K + gammaStar)))
  
  return(list(eta_star = nstar, gamma_star = gammaStar, K = K, beta = beta))
}


# Tests
calc_beta(chi = c(0.9), temp = 27, vpd = 1000, z = 1000)


