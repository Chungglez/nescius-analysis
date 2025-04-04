from nescius_heat_transfer.unit_conversions import *
import numpy as np
from math import log

# Inputs
burn_time = 30
fo_ratio = 1.1
m_dot = 0.08687942639039442
fuel_density = 789
ox_density = 1141
tank_diameter = in2m(3.26)

# Proppellant masses
m_dot_ox = m_dot*fo_ratio/(1+fo_ratio)
m_dot_fuel = m_dot*1/(1+fo_ratio)
m_fuel = burn_time*m_dot_fuel
m_ox = burn_time*m_dot_ox

# Tank volumes and sizing
Vol_fuel = m_fuel/fuel_density
Vol_ox = m_ox/ox_density
L_cyl_fuel = Vol_fuel/d2area(tank_diameter)
L_cyl_ox = Vol_ox/d2area(tank_diameter)
print(f"LOX tank length: {m2in(L_cyl_ox)} (in)")
print(f"Ethanol tank length: {m2in(L_cyl_fuel)} (in)")

# Capacitive level sensor sizes
R_inner = in2m(.25)
R_outer = in2m(.5)
e_0 = 1/(4e-7*pi*299792458**2)
e_r_fuel = 24.3
e_r_ox = 1.5
print(e_0)

C_fuel = 2*pi*e_0*e_r_fuel*L_cyl_fuel/log(R_outer/R_inner)
C_empty = 2*pi*e_0*L_cyl_fuel/log(R_outer/R_inner)

print(f"Filled (pF): {C_fuel*1e12}")
print(f"Empty (pF): {C_empty*1e12}")


C_ox = 2*pi*e_0*e_r_ox*L_cyl_fuel/log(R_outer/R_inner)

