import math
from math import floor

from nescius_heat_transfer.unit_conversions import *
import itertools

from valve_design.solenoid_sizing import voltage

# FIXED VARIABLES
#-----------------------------------------------------------------------------------------------------------------------
# Value                       | Units       | Description
#-----------------------------+-------------+------------------
mu_0 = 4e-7*pi                # H/m         | Permeability of free space
mu_r = 100                    # N/A         | Relative permeability of 410 stainless steel
resistivity = 1.724e-8        # Ω-m         | Resistivity of copper
fluid_pressure = psia2Pa(500) # Pa          | Maximum absolute pump pressure
stroke = in2m(1/16)            # m           | Trave`l distance of the solenoid plunger
bobbin_thickness = mm2m(1)    # m           | Thickness of the bobbin that the coild is wound around
bot_thickness = in2m(3/16)   # m           | Thickness of the circular steel disc at the bottom of the solenoid
top_thickness = in2m(3/16)    # m           | Thickness of the circular steel disc at the top of the solenoid
outer_tube_thickness = in2m(.065)   # m           | Wall thickness of the steel tube surrounding the solenoid
core_casing_thickness = in2m(3/64)
oring_flange_width = in2m(.1875)
    # m           | Inner diameter of the solenoid coils

# VARIABLE PARAMETERS
#-----------------------------------------------------------------------------------------------------------------------
# Each combination of every option will be anaylzed to find the best performance within the range of possibilities
voltages = [5, 12, 24] # V
solenoid_diameters = [in2m(2)] # m
solenoid_heights = [in2m(2)] # m
core_IDs = [in2m(.875), in2m(1), in2m(1.0623), in2m(1.125), in2m(1.1875), in2m(1.25), in2m(1.3125), in2m(1.375)]
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
options = list(itertools.product(voltages, solenoid_diameters, solenoid_heights, wire_sizes,core_IDs))
best_option = None
best_score = float('-inf')
counter = 1

# Loop through each option, and find the best one
for voltage, solenoid_diameter, solenoid_height, wire_size, core_ID in options:
    # Unpack and do simple dimension calculations
    solenoid_ID = core_ID + 2 * oring_flange_width + 2 * bobbin_thickness
    copper_diam, wire_diam, max_current = wire_sizes[wire_size]
    coil_height = solenoid_height - top_thickness - bot_thickness - 2*bobbin_thickness
    coil_width = (solenoid_diameter - solenoid_ID - 2*outer_tube_thickness)/2
    print(m2in(coil_width))

    # Calculate the number of coils
    n_c = floor(coil_height/wire_diam - 1)
    m_c = floor(coil_width/wire_diam)
    N = n_c*m_c
    if N <= 0:
        continue

    # Calculate electrical quantities
    l_avg = pi*(solenoid_ID+coil_width)
    resistance = resistivity * l_avg * N / d2area(copper_diam) * (234.5+80)/(234.5+20)
    current = voltage/resistance
    power_loss = current**2 * resistance

    # Magnetic circuit
    F_mag = N * current
    R_top = math.log(solenoid_diameter/core_ID)/(mu_r*mu_0*2*pi*top_thickness)
    R_tube = (solenoid_height-top_thickness-bot_thickness)/(mu_0*mu_r*(d2area(solenoid_diameter) - d2area(solenoid_diameter-2*outer_tube_thickness)))
    R_core = (solenoid_height-stroke-in2m(.375))/(mu_0*mu_r*d2area(solenoid_ID-2*bobbin_thickness))
    R_gap_radial = math.log((core_ID+2*core_casing_thickness)/core_ID)/(mu_0*2*pi*in2m(.375))
    R_gap_axial = stroke/(mu_0*(d2area(core_ID) - d2area(in2m(0.25))))
    R_plunger = in2m(.5)/(mu_0*mu_r*(d2area(core_ID)-d2area(in2m(0.375))))
    R_bot = math.log(solenoid_diameter/solenoid_ID)/(mu_r*mu_0*2*pi*bot_thickness)

    R_m = R_top+R_tube+R_core+R_gap_radial+R_gap_axial+R_plunger+R_bot
    R_gap = R_gap_radial+R_gap_axial
    R_steel = R_m-R_gap
    #print(
    #    f"Top={R_top / R_m * 100:.1f}%, Tube={R_tube / R_m * 100:.1f}%, Core={R_core / R_m * 100:.1f}%, GapRad={R_gap_radial / R_m * 100:.1f}%, GapAx={R_gap_axial / R_m * 100:.1f}%, Plung={R_plunger / R_m * 100:.1f}%, Bot={R_bot / R_m * 100:.1f}%")

    flux = F_mag/(R_m)

    # Force on solenoid
    A_plunger = d2area(core_ID)-d2area(in2m(.25))
    B_plunger = flux/A_plunger
    F_plunger = B_plunger**2*A_plunger/(2*mu_0)
    F_compare=1/2*flux**2 * (1/(mu_0*(d2area(core_ID) - d2area(in2m(0.25))))  - math.log((core_ID+2*core_casing_thickness)/core_ID)/(mu_0*2*pi*in2m(.375)**2))

    print(f"{F_plunger} {F_compare}")


    L = N**2/R_m
    i1 = voltage/L*.001
    # Evaluate (and save) result
    score = F_compare
    if current > 8*max_current:
        continue    # Don't allow results that are close to the limit of the wire
    if score > best_score:
        best_score = score
        best_option = (voltage,solenoid_diameter,solenoid_height, wire_size,F_compare, current, N, solenoid_diameter, solenoid_height, B_plunger, core_ID, flux)
    counter+=1

# Unpack and display best result
(voltage,solenoid_diameter,solenoid_height, wire_size,F_plunger, current, N, solenoid_diameter, solenoid_height, B_plunger, core_ID, flux) = best_option
print(f"Turns: {N} \n"
      f"Current: {current} \n"
      f"Voltage: {voltage} \n"
      f"Force (lbs): {N2lb(F_plunger)} \n"
      f"Wire size: {wire_size} \n"
      f"Height (in): {m2in(solenoid_height)} \n"
      f"Diameter (in): {m2in(solenoid_diameter)} \n"
      f"B: {B_plunger} \n"
      f"Core ID (in): {m2in(core_ID)} \n"
      f"Flux: {flux} \n")
