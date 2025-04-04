from math import exp
from scipy.optimize import fsolve
import numpy as np
from unit_conversions import *

L = .12 # m
f = .03
d_hyd = mm2m(3.7) # m
rho = 789 # kg/m^3
m_dot = 0.07240876657074044 # kg/s
n_channels = 20
V = m_dot/n_channels/(rho*d2area(d_hyd)) # m/s
visc = 0.0011
Re = rho*V*d_hyd/visc

deltaP = .03*L/d_hyd * rho * V**2 /2 # N/m^2


E_316 = 1.8e11
cte_316 = 16e-6
deltaT = 430-360
k_316 = 15
poisson_316 = 0.25
S = E_316*cte_316*deltaT/(2*(1-poisson_316))
print(Pa2MPa(S))