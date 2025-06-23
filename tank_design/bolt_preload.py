from nescius_heat_transfer.unit_conversions import *
from math import ceil

# 10-32 screw clamp force
major_diam = .190                        # in
minor_diam = .156                        # in
pitch_diam = (major_diam + minor_diam)/2 # in
screw_area = d2area(pitch_diam)          # in^2

yield_strength = 30000                   # psi
proof_strength = 0.85 * yield_strength   # psi
proof_load = proof_strength * screw_area # lb
preload = 0.8 * proof_load               # lb
print(f"#10-32 Stainless Steel Screw Preload: {preload:.2f} (lb)")

# Thermal contraction
ptfe_length = 1/16
delta_T = -297.33 - 70  # Lox temp - assembly temp
cte = 70e-6
delta_L = cte*ptfe_length*delta_T
print(f"PTFE shrinkage: {delta_L:.4f} (in)")
screw_length = .588
cte = 8.89e-6
delta_L = cte*screw_length*delta_T
print(f"Screw shrinkage: {delta_L:.4f} (in)")

# Belleville Washer
# https://www.mcmaster.com/9713K12/
washer_load = 125
washer_deflection = .003
n_parallel = ceil(preload/washer_load)
n_series = ceil(abs(delta_L)/washer_deflection)
print(f"Belleville waashers in series: {n_series}")
print(f"Belleville waashers in parallel: {n_parallel}")