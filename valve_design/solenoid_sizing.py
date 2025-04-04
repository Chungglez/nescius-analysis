import math
from math import floor

from nescius_heat_transfer.unit_conversions import *
import itertools

# FIXED VARIABLES
#-----------------------------------------------------------------------------------------------------------------------
# Value                       | Units       | Description
#-----------------------------+-------------+------------------
mu_0 = 4e-7*pi                # H/m         | Permeability of free space
mu_r = 700                    # N/A         | Relative permeability of 410 stainless steel
resistivity = 1.724e-8        # Ω-m         | Resistivity of copper
fluid_pressure = psia2Pa(500) # Pa          | Maximum absolute pump pressure
stroke = in2m(.25)            # m           | Travel distance of the solenoid plunger
housing_height = in2m(0.875)  # m           | Dimension of the plunger housing cap that contains the fluid and plunger
bobbin_thickness = mm2m(1)    # m           | Thickness of the bobbin that the coild is wound around
bot_thickness = in2m(.1875)   # m           | Thickness of the circular steel disc at the bottom of the solenoid
top_thickness = in2m(.125)    # m           | Thickness of the circular steel disc at the top of the solenoid
tube_thickness = in2m(.125)   # m           | Wall thickness of the steel tube surrounding the solenoid
solenoid_ID = in2m(.625)      # m           | Inner diameter of the solenoid coils
solenoid_ID += 2*bobbin_thickness

# VARIABLE PARAMETERS
#-----------------------------------------------------------------------------------------------------------------------
# Each combination of every option will be anaylzed to find the best performance within the range of possibilities
voltages = [5, 12, 24] # V
solenoid_diamneters = [in2m(1.375), in2m(1.5), in2m(1.625),in2m(1.75),in2m(2)] # m
solenoid_heights = [in2m(1.25),in2m(1.5)] # m
wire_sizes = {
 # AWG: (Copper Diameter (m), Wire Diameter (m), Max Current (A))
    20: (in2m(.0320), in2m(.0346), 1.5),
    22: (in2m(.0254), in2m(.0276), .92),
    24: (in2m(.0201),in2m(.0223), .577),
    26: (in2m(.0159),in2m(.0178), .361),
    28: (in2m(.0126),in2m(.0144), .226),
    30: (in2m(.01),in2m(.0117), .142),
    32: (in2m(.008),in2m(.0094), .091)
}

# Generate every combination of the parameters
options = list(itertools.product(voltages, solenoid_diamneters, solenoid_heights, wire_sizes))
best_option = None
best_score = float('-inf')
counter = 1

# Loop through each option, and find the best one
for voltage, solenoid_diamneter, solenoid_height, wire_size in options:
    # Unpack and do simple dimension calculations
    copper_diam, wire_diam, max_current = wire_sizes[wire_size]
    coil_height = solenoid_height - top_thickness - bot_thickness - 2*bobbin_thickness
    coil_width = solenoid_diamneter - solenoid_ID - 2*bobbin_thickness - 2*tube_thickness

    # Calculate the number of coils
    n_c = floor(coil_height/wire_diam - 1)
    m_c = floor(coil_width/wire_diam)
    N = n_c*m_c

    # Calculate electrical quantities
    l_avg = pi*(solenoid_ID+coil_width)
    resistance = resistivity * l_avg * N / d2area(copper_diam) * (234.5+80)/(234.5+20)
    current = voltage/resistance
    power_loss = current**2 * resistance

    # Magnetic circuit
    F_mag = N * current
    R_top = math.log(solenoid_diamneter/(solenoid_ID-2*bobbin_thickness))/(mu_r*mu_0*2*pi*top_thickness)
    R_tube = (solenoid_height-top_thickness-bot_thickness)/(mu_0*mu_r*pi/4*(solenoid_diamneter**2 - (solenoid_diamneter-2*tube_thickness)**2))
    R_core = (solenoid_height-housing_height)/(mu_0*mu_r*pi/4*(solenoid_ID-2*bobbin_thickness)**2)
    R_gap_radial = math.log((solenoid_ID-2*bobbin_thickness)/in2m(0.5))/(mu_0*2*pi*bot_thickness)
    R_gap_axial = in2m(.30503357)/(mu_0**pi/4*(solenoid_ID-2*bobbin_thickness)**2)
    R_plunger = in2m(.56996643)/(mu_0*mu_r*pi/4*(in2m(0.5)**2-in2m(0.25)**2))
    R_bot = math.log(solenoid_diamneter/solenoid_ID)/(mu_r*mu_0*2*pi*bot_thickness)
    R_m = R_top+R_tube+R_core+R_gap_radial+R_gap_axial+R_plunger+R_bot
    flux = F_mag/R_m

    # Force on solenoid
    B_plunger = flux/(pi/4*(in2m(0.5)**2-in2m(0.25)**2))
    F_plunger = B_plunger**2*(pi/4*(in2m(0.5)**2-in2m(0.25)**2))/(2*mu_0)

    # Evaluate (and save) result
    score = F_plunger - 10*power_loss
    if current*1.3 > max_current:
        continue    # Don't allow results that are close to the limit of the wire
    if score > best_score:
        best_score = score
        best_option = (voltage,solenoid_diamneter,solenoid_height, wire_size,F_plunger, current, N, solenoid_diamneter, solenoid_height, B_plunger)
    counter+=1

# Unpack and display best result
(voltage,solenoid_diamneter,solenoid_height, wire_size,F_plunger, current, N, solenoid_diamneter, solenoid_height, B_plunger) = best_option
print(f"Turns: {N} \n"
      f"Current: {current} \n"
      f"Voltage: {voltage} \n"
      f"Force (lbs): {N2lb(F_plunger)} \n"
      f"Wire size: {wire_size} \n"
      f"Height (in): {m2in(solenoid_height)} \n"
      f"Diameter (in): {m2in(solenoid_diamneter)} \n"
      f"B: {B_plunger}")
