import matplotlib.pyplot as plt
import numpy as np
from gas import Gas, mix_func
from unit_conversions import *

# STEP 1 - Obtain data required for calcultions
#-----------------------------------------------------------------------
m_dot_c = 0.01          # kg/s      Mass flow rate of coolant fluid
m_dot_g = 0.121         # kg/s      Mass flow rate of hot combustion gas
P_cc = psia2Pa(420)     # Pa        Combustion chamber pressure
D_cc = mm2m(42.11)      # m         Diameter of combustion chamber
of_ratio = 1.5          # N/A       Weight ratio of fuel to oxidizer

# STEP 2 - Calculate free stream gas properties
#-----------------------------------------------------------------------
T_cc = 3222.64          # K         Combustion gas temperature

# Gases fractions from combustion from CEA
g_CO = Gas("CO",0.3094, g2kg(28.0101))
g_CO2 = Gas("CO2",0.27809, g2kg(44.0095))
g_H2O = Gas("H2O",0.36675, g2kg(18.01528))
g_OH = Gas("OH",0.02538, g2kg(17.00734))
g_O2 = Gas("O2",0.00805, g2kg(31.99880))
g_O = Gas("O",0.00230,g2kg(15.99940))
g_H2 = Gas("H2",0.00916, g2kg(2.01588))
g_H = Gas("H",0.0008,g2kg(1.007940))
gases = [g_CO, g_CO2, g_H2O, g_OH, g_O2, g_O, g_H2, g_H]

# Scale mass fractions so the total is 1.0
mf_sum = sum(gas.mass_fraction for gas in gases)
print(f"Mass fraction before scaling: {mf_sum}")
for gas in gases:
    gas.mass_fraction /= mf_sum

# Calculate specific heats of gases using table A-9 from:
# Norton, Appendices, in: F. Kreith, D.Y. Goswami (Eds.), The CRC Handbook of Mechanical Engineering, CRC Press, 2005.
g_CO.calc_cp(0.297,2.654,2.226,-1.146,2.851,-.2762, T_cc)
g_CO2.calc_cp(0.189,3.247,5.847,-3.412,9.469,-1.009, T_cc)
g_H2O.calc_cp(0.462,2.798,2.693,-0.5392,-0.01783,0.09027, T_cc)
g_OH.calc_cp(0.489,3.229,0.2014,0.4357,-2.043,0.2696, T_cc)
g_O2.calc_cp(0.26,3.156,1.809,-1.052,3.190,-0.3629, T_cc)
g_O.calc_cp(0.520,2.662,-.3051,.2250,-.7447,0.09383, T_cc)
g_H2.calc_cp(4.12,3.717,-0.9220,1.221,-4.328, 0.5202, T_cc)
g_H.calc_cp(8.25,2.567,-0.1509,0.1219,-0.4184,0.05182, T_cc)

# Converting from mass fraction to mol fractions
M_tot = sum(gas.mass_fraction/gas.M for gas in gases)
for gas in gases:
    gas.mol_fraction = gas.mass_fraction/gas.M/M_tot
# Take mass averaged specific heat for the mixture
c_pg = sum(gas.mol_fraction*gas.cp for gas in gases)
print(f"Gas mixture specific heat: {c_pg}") # NASA CEA value: c_pg = 2.2341 pretty close to calculation

# Calculating viscosities
g_CO.calc_visc(23.811,5.3944e-1, -1.5411e-4,T_cc)
g_CO2.calc_visc(11.811, 4.9838e-1, -1.0851e-4, T_cc)
g_H2O.calc_visc(-36.826, 4.2900e-1, -1.6200e-5,T_cc)
g_O2.calc_visc(44.224,5.6200e-1,-1.1300e-4,T_cc)
g_O.visc = g_O2.visc/2 # NO DATA, ASSUMING HALF AS H2
g_H2.calc_visc(27.758,2.1200e-1,-3.2800e-5,T_cc)
g_H.visc = g_H2.visc/2 # NO DATA, ASSUMING HALF AS H2
g_OH.visc = (g_H2.visc+g_O2.visc)/2 # NO DATA, ASSUMING AVERAGE OF O2 and H2
# The viscosities that were assumed only make up ~3% of the mass of the gas


# Wikes method for finding viscosity of mixture
visc_g = sum(
    gas_i.mol_fraction*gas_i.visc/
    sum(
        gas_j.mol_fraction*mix_func(gas_i,gas_j)
        for gas_j in gases
    )
    for gas_i in gases
)
print(f"Gas mixture viscosity: {visc_g} kg/m-s")

# Calculating thermal conductivity of gas


# c_pg = 2.2341              # kJ/kg*K    Specific heat of combustion gas at constant pressure
u_g = ...
Pr_g = ...
rho_g = 2.4041

# STEP 3 - Obtain liquid coolant properties
#-----------------------------------------------------------------------

# Enthalpy of vaporization, it isn't published at 420 psia
# Instead, it will calculate it from vapor pressure curve and the Clausius-Clapeyron equation

# STEP 4 - Calculate Reynolds number
#-----------------------------------------------------------------------

# STEP 5 - Calculate the friction factor
#-----------------------------------------------------------------------

# STEP 6 - Calculate h_o
#-----------------------------------------------------------------------

# STEP 7 - Calculate h_fgs
#-----------------------------------------------------------------------

# STEP 8 - Calculate Q_rad
#-----------------------------------------------------------------------

# STEP 9 - Correct for transpiration
#-----------------------------------------------------------------------

# STEP 10 - Calculate Q_tot and m_dot_v
#-----------------------------------------------------------------------

# STEP 11 - Calculate the entrainment fraction
#-----------------------------------------------------------------------

# STEP 12 - Calculate the liquid film cooled length
#-----------------------------------------------------------------------
