from nescius_heat_transfer.unit_conversions import *
import numpy as np
from math import log

# Inputs
burn_time = 30
of_ratio = 1.6
m_dot = 0.12543*2
fuel_density = 789
ox_density = 1141
tank_diameter = in2m(4.26)

# Proppellant masses
m_dot_ox = m_dot*of_ratio/(1+of_ratio)
m_dot_fuel = m_dot*1/(1+of_ratio)
m_fuel = burn_time*m_dot_fuel
m_ox = burn_time*m_dot_ox
fo_ratio = m_fuel/m_ox
print(fo_ratio)

# Tank volumes and sizing
Vol_fuel = m_fuel/fuel_density
Vol_ox = m_ox/ox_density
L_cyl_fuel = Vol_fuel/d2area(tank_diameter)
L_cyl_ox = Vol_ox/d2area(tank_diameter)
print(f"LOX tank minimum length: {m2in(L_cyl_ox)} (in)")
print(f"Ethanol tank minimum length: {m2in(L_cyl_fuel)} (in)")

L_tank = in2m(18)

# Capacitive level sensor sizes
R_inner = in2m(.3125/2)
R_outer = in2m(.493/2)
e_0 = 1/(4e-7*pi*299792458**2)
e_r_fuel = 24.3
e_r_ox = 1.5
e_r_water = 80

# Coaxial capacitance calc
C_fuel = 2*pi*e_0*e_r_fuel*L_tank/log(R_outer/R_inner)
C_ox = 2*pi*e_0*e_r_ox*L_tank/log(R_outer/R_inner)
C_water = 2*pi*e_0*e_r_water*L_tank/log(R_outer/R_inner)
C_empty = 2*pi*e_0*L_tank/log(R_outer/R_inner)

# Stray capacitance from wires
L_wires = 0.25 # m
r_wires = mm2m(.812/2) # 20 GA
wire_seperation = mm2m(5)
C_wire = pi * e_0/log(wire_seperation/r_wires)
print(C_wire*1e12)

# Printing results
print(f"Filled fuel (pF): {C_fuel*1e12}")
print(f"Filled water (pF): {C_water*1e12}")
print(f"Filled lox (pF): {C_ox*1e12}")
print(f"Empty (pF): {C_empty*1e12}")
print(f"Wires (pF): {C_wire*1e12}\n")

# Fuel bring into range...
C_series = 47e-12
C_new_full = 1/(1/(C_fuel+C_wire) + 1/C_series)
C_new_empty = 1/(1/(C_empty+C_wire) + 1/C_series)
print(f"Fuel input range (pF): {(C_new_full-C_new_empty)*1e12}")
print(f"Fuel required offseet (pF): {C_new_full*1e12-15}")

# Water bringing into range...
C_series = 33e-12
C_new_full = 1/(1/(C_water+C_wire) + 1/C_series)
C_new_empty = 1/(1/(C_empty+C_wire) + 1/C_series)

print(f"Water input range (pF): {(C_new_full-C_new_empty)*1e12}")
print(f"Water required offseet (pF): {C_new_full*1e12-15}")

# LOX bringing into range...
C_series = 470e-12 # No capacitor in series
C_new_full = 1/(1/(C_ox+C_wire) + 1/C_series)
C_new_empty = 1/(1/(C_empty+C_wire) + 1/C_series)

print(f"LOX input range (pF): {(C_new_full-C_new_empty)*1e12}")
print(f"LOX required offseet (pF): {C_new_full*1e12-15}")