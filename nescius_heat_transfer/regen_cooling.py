import math
import matplotlib.pyplot as plt
from dask.array import indices
from scipy.optimize import fsolve
from math import exp

from unit_conversions import *
from util import *
import numpy as np
from scipy.integrate import quad

# ENGINE PARAMETERS
#-----------------------------------------------------------------------------------------------------------------------

# Value                    | Units       | Description
#--------------------------+-------------+------------------
g = 9.80665                # m/s^2       | Acceleration due to gravity
gamma = 1.1598             # N/A         | Specific heat ratio
R = 373.8                  # J/kg-K      | Gas constant (From RPA)
visc_g = 1.0405e-4         # Pa-s        | Gas viscosity
Cp_g = gamma*R/(gamma-1)   # J/kg-K      | Gas specific heat
Tc_ns = 1410            # K           | Combustion temperature (From CEA)
Pr_g = 0.52                # N/A         | Gas Prandtl's number (From RPA, NASA equation didn't match well)
pc_ns = psia2Pa(420)       # N/m^2       | Nozzle stagnation pressure
F = lb2N(152)               # N           | Engine thrust
pe = psia2Pa(p_atm)        # N/m^2       | Nozzle exist pressure, atmospheric presure
fo_ratio = 1.1             # N/A         | Weight ratio of fuel to oxidizer
contraction_ratio = 22.5   # N/A         | Ratio of combustion chamber area to throaat area
Lc = 1.667                 # m           | Characteristic length, aka time for combustion
theta_conv = 40            # deg         | Nozzle converging half cone angle
r_2ratio = 0.65

# ENGINE SIZING
#-----------------------------------------------------------------------------------------------------------------------
# Calculate exit velo   city (m/s) from combustion temperature and pressure
v_e = sqrt((2 * gamma) / (gamma - 1) * R * Tc_ns * (1 - (pe / pc_ns) ** ((gamma - 1) / gamma)))
Isp = v_e / g # s
# Calculate total mass flow rate
m_dot = F / v_e # kg/s
m_dot_ox = m_dot*fo_ratio/(1+fo_ratio) # kg/s
m_dot_fuel = m_dot*1/(1+fo_ratio) # kg/s
print(f"m_dot_fuel {m_dot_fuel}")

# Calculate throat and exit areas    m^2
A_t = m_dot * sqrt(Tc_ns) / (
        pc_ns * sqrt(gamma / R) * ((gamma + 1) / 2) ** (
        -(gamma + 1) / (2 * (gamma - 1))
    )
)
A_e = A_t * (2 / (gamma + 1)) ** (1 / (gamma - 1)) * (pc_ns / pe) ** (1 / gamma) / sqrt(
    (gamma + 1) / (gamma - 1) * (1 - (pe / pc_ns) ** ((gamma - 1) / gamma))
)

# ENGINE GEOMETREY FUNCTION
#-----------------------------------------------------------------------------------------------------------------------
# Parabolic nozzle:
exp_ratio = A_e / A_t
# Parabola starting and exit angles
divg_angle_func = lambda e: 0.0413 * e ** 2 + 0.475 * e + 20.888
exit_anlge_funct = lambda e: 0.0188 * e ** 2 - 0.9092 * e + 15.663
theta_e = 12.9  # exit_anlge_funct(exp_ratio)
theta_n = 22.9  # divg_angle_func(exp_ratio)
R_t = area2r(A_t)
D_t = 2*R_t
R_e = area2r(A_e)
R_c = area2r(A_t * contraction_ratio)

R_5 = lambda x: -np.sqrt((0.382 * R_t) ** 2 - x ** 2) + 1.382 * R_t
x_5 = 0.382 * R_t * sind(theta_n)
x_4 = 0
R_4 = lambda x: -np.sqrt((1.5 * R_t) ** 2 - x ** 2) + 2.5 * R_t
x_3 = 1.5 * R_t * sind(-theta_conv)
r_2max = (R_c-R_4(x_3))/(1-cosd(theta_conv))
r_2 = r_2max*r_2ratio
R_3 = lambda x: -tand(theta_conv) * (x - x_3) + R_4(x_3)
x_2 = x_3 - (R_c - R_4(x_3) - r_2*(1 - cosd(theta_conv))) / tand(theta_conv)
x_1 = x_2 - r_2 * sind(theta_conv)
R_2 = lambda x: np.sqrt(r_2 ** 2 - (x - x_1) ** 2) + (R_c-r_2)

def R_conv(z):
    return np.piecewise(
        z,
        [z < x_2, (z >= x_2) & (z < x_3), z >= x_3],
        [R_2, R_3, R_4]
    )

A_conv = lambda x: pi * R_conv(x)**2
Vol_conv = quad(A_conv,x_1,x_4)[0]
Vol_cyl = Lc*A_t - Vol_conv
L_cyl = Vol_cyl/r2area(R_c)
print(f"L_cyl {m2mm(L_cyl)}")
x_0 = x_1-L_cyl

M1 = tand(theta_n)
M2 = tand(theta_e)
x1 = x_5
y1 = R_5(x1)
y2 = area2r(A_e)

c_func = lambda c: np.sqrt(2*M2*(y2-c))*np.sqrt(M2/(2*M1)*(y2-c))+c-y1
c = fsolve(c_func,0)[0]
b = M2/(2*M1)*(y2-c)-x1
a = sqrt(2*M2*(y2-c))
x_6 =((y2-c)/a)**2-b
R_6 = lambda x: a*np.sqrt(x+b)+c

def R_engine(z):
    return np.piecewise(
        z,
        [z < x_1, (z >= x_1)& (z < x_2), (z >= x_2) & (z < x_3), (z >= x_3) & (z < x_4),(z >= x_4) & (z < x_5), z>=x_5],
        [R_c, R_2, R_3, R_4, R_5,R_6]
    )
print(f"{m2mm(D_t)} {m2mm(R_c*2)} {m2mm(R_e)} {m2mm(L_cyl)}")

x_list = np.linspace(x_6, x_0, 1000)
plt.plot(x_list, R_engine(x_list))
plt.plot(x_list, -R_engine(x_list))
plt.axis('equal')
plt.xlabel("Nozzle Position (m)")
plt.ylabel("Nozzle Radius (m)")
plt.title(f"Contraction Ratio: {contraction_ratio}")
plt.show()

# THERMODYNAMIC PROPERTIES
#-----------------------------------------------------------------------------------------------------------------------
A_list = r2area(R_engine(x_list))
area_mach_eq = lambda M,A: A / A_t - pow(2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M ** 2),
                                                            (gamma + 1) / (2 * (gamma - 1))) / M
M_list = np.array([
    fsolve(
        area_mach_eq,
        0.1 if x < x_4 else 2,
        args=(A)
    )[0]
    for A, x in zip(A_list, x_list)
])

# Find gas properties using isentropic relations
recovery_factor = Pr_g**0.33
T_list = Tc_ns*pow(1+(gamma-1)/2*M_list**2,-1)
Taw_list = Tc_ns*(1+recovery_factor*(gamma-1)/2*M_list**2)/(1+(gamma-1)/2*M_list**2)

T_comp = Taw_list/T_list
p_list = pc_ns*pow(1+(gamma-1)/2*M_list**2,-gamma/(gamma-1))
rho_c = pc_ns/(R*Tc_ns)
rho_list = rho_c*pow(1+(gamma-1)/2*M_list**2,-gamma/(gamma-1))
Vel_list = np.sqrt(gamma*R*T_list)*M_list

# CHAMBER ACOUSTIC MODES
#-----------------------------------------------------------------------------------------------------------------------
L_c = x_4 - x_0             # Length of combustion chamber
A_c = sqrt(gamma*R*Tc_ns)   # Speed of sound m/s
N_L = A_c/(2*L_c)           # Longitudinal frequency, Hz
N_T = 0.59*A_c/(2*R_c)      # Transverse frequency, Hz
N_R = 1.22*A_c/(2*R_c)      # Radial frequency

# Display results
print(f"Longitudinal frequency: {N_L} Hz")
print(f"Longitudinal frequency: {N_T} Hz")
print(f"Longitudinal frequency: {N_R} Hz")

# REGEN COOLING HEAT TRANSFER
#-----------------------------------------------------------------------------------------------------------------------

# Value                      | Units       | Description
#----------------------------+-------------+------------------
k_316 = 15                   # W/m-K       | 316L wall thermal conductivity
E_316 = 1.8e11               # Pa          | Modulus of elasticity for 316L
poisson_316 = 0.25           # N/A         | Poisson's ratio for 316L
cte_316 = 16e-6              # m/m-K         | Coeffecient of thermal expansion for 316L
k_co = 0.167                 # W/m-k       | Ethanol coolant thermal conductivity
shell_thick_i = mm2m(1)      # m           | Thickness of inner chamber wall
shell_thick_o = mm2m(1)      # m           | Outer chamber wall thickness
channel_w = mm2m(1)          # m           | Thickness of channel divider
channel_h = mm2m(1.5)          # m           | Channel hieght
n_channels = 20              # Int         | Number of cooling channels
theta_chan = 2*pi/n_channels # rad         | Angle a channel takes up
Tco = 293.15                 # K           | Coolant inlet temperature

# Viscosity of ethanol
# Source: https://en.wikipedia.org/wiki/Temperature_dependence_of_viscosity
visc_ethanol = lambda T: (0.00201 * exp(1614/T + 0.00618*T - 1.132E-5*T**2))/1000   # Pa-s
# Isobaric heat capacity of liquid ethanol
# Source (linear fit of straight portion):
Cp_ethanol = lambda T: (0.0109*T-0.6907)*1000   # J/kg-K
# Density of ethanol
# Source: https://www.engineeringtoolbox.com/ethanol-ethyl-alcohol-density-specific-weight-temperature-pressure-d_2028.html
rho_table_T = np.array([250, 280, 310,340,370,400,423.9,430,460,490,501.4,513.9,])  # K
rho_table_rho = np.array([825.2,800.5,775.3,747.8,716.7,680.6,647.2,637.7,584,507.1,461.3,276]) # kg/m^3
rho_ethanol = lambda T: interp(T,rho_table_T, rho_table_rho)
hg = None

def friction_factor(D_h, Re_b, Re_w, T_wc, T_co) -> float:
    if Re_b < 3000:
        return 64/Re_b
    epsiolon = microm2m(39.8) #source: Rz of 316L https://xometry.pro/en-eu/articles/3d-printing-surface-roughness/
    colebrook = lambda fr: 1/np.sqrt(fr)+2.0*np.log10(
        epsiolon/(3.72*D_h)+2.51/(Re_b*np.sqrt(fr))
    )
    fdr = fsolve(colebrook,np.array([.04]))[0]
    return fdr*pow(T_wc/T_co,-0.6-5.6*pow(Re_w,-0.38))

# How much do fins chnge heat transfer compared to a flat plate?
def fin_factor(x, h_c) -> float:
    R_wc = R_engine(x)+shell_thick_i
    R_wc2 = R_wc + channel_h
    s = (R_wc*theta_chan - channel_w)/2
    L2 = (R_wc2*theta_chan - channel_w)/2

    m1 = sqrt(2 * h_c / (k_316 * channel_w))
    m2 = sqrt(2 * h_c / (k_316 * shell_thick_o*2))
    c_1 = (
                  exp(2 * channel_h * m1) * (b * m2 * exp(2 * L2 * m2) - b * m2 + m1 * channel_w* exp(2 * L2 * m2) + m1 * channel_w)
          ) / (
                  -b * m2 * exp(2 * channel_h * m1) - b * m2 * exp(2 * L2 * m2) + b * m2 * exp(2 * channel_h * m1 + 2 * L2 * m2)
                  + b * m2 + m1 * channel_w* exp(2 * channel_h * m1) + m1 * channel_w* exp(2 * L2 * m2)
                  + m1 * channel_w* exp(2 * channel_h * m1 + 2 * L2 * m2) + m1 * channel_w
          )
    c_2 = (
              (-b * m2 * exp(2 * L2 * m2) + b * m2 + m1 * channel_w* exp(2 * L2 * m2) + m1 * channel_w)
          ) / (
                  -b * m2 * exp(2 * channel_h * m1) - b * m2 * exp(2 * L2 * m2) + b * m2 * exp(2 * channel_h * m1 + 2 * L2 * m2)
                  + b * m2 + m1 * channel_w* exp(2 * channel_h * m1) + m1 * channel_w* exp(2 * L2 * m2)
                  + m1 * channel_w* exp(2 * channel_h * m1 + 2 * L2 * m2) + m1 *channel_w
          )
    eta_f = 2 * s / (2 * s +channel_w) + k_316 * channel_w* (c_1 * m1 - c_2 * m2) / (h_c * (2 * s +channel_w))
    return eta_f

def heat_trans_eqs(variables, x, Tco, i):
    q, Twg, Twc = variables
    Taw = interp(x, x_list,Taw_list)
    A = interp(x, x_list, A_list)
    M = interp(x, x_list, M_list)

    # COOLANT CONVECTION COEFFECIENT
    #-------------------------------------------------------------------------------------------------------------------
    # Find the hydraulic diameter of the cooling channel
    area_channel = theta_chan/2*((R_engine(x)+shell_thick_i+channel_h)**2 - (R_engine(x)+shell_thick_i)**2) -channel_h* channel_w    # m^2
    perim_channel = 2*channel_h + theta_chan*(R_engine(x)+shell_thick_i+channel_h) + theta_chan*(R_engine(x)+shell_thick_i) - 2* channel_w # m
    d_hyd = 4 * area_channel/perim_channel              # m     Hydraulic diameter

    # Determine physical properties of coolant
    visc_co = visc_ethanol(Tco)                         # Pa-s   Coolant bulk Viscosity
    Cp_co = Cp_ethanol(Tco)                             # J/kg-K Coolant Isobaric heat capacity
    Pr_co = Cp_co * visc_co / k_co                      # N/A    Coolant Prandtl's number
    V_co = 2* m_dot_fuel/n_channels/rho_ethanol(Tco)/area_channel    # m/s    Coolant channel velocity
    Re_co = rho_ethanol(Tco)*V_co*d_hyd/visc_co                   # N/A    Coolant Reynold's number
    Re_w = rho_ethanol(Twc) * V_co * d_hyd / visc_ethanol(Twc)

    if Re_co < 3000:
        hc = k_co/d_hyd*4.36
    elif Re_co < 10000:
        f = friction_factor(d_hyd,Re_co,Re_w,Twc,Tco)
        hc = k_co/d_hyd*((f/8)*(Re_co-1000)*Pr_co)/(1+12.7*(f/8)**.5*(Pr_co**(2/3)-1))
    # Sieder–Tate equation to find coolant heat transfer
    else:
        hc = k_co/d_hyd * 0.027 * Re_co ** 0.8 * Pr_co ** 0.4 * (visc_co/visc_ethanol(Twc)) ** 0.14  # W/m^2-k

    # GAS SIDE CONVECTION COEFFECIENT
    # -------------------------------------------------------------------------------------------------------------------
    sigma = 1/((1/2*Twg/Tc_ns*(1+(gamma-1)/2*M**2)+1/2)**0.68*
               (1+(gamma-1)/2*M**2)**0.12
               )
    hg = h2SI((0.026/(m2in(D_t)**0.2)*
          (visc2Imp(visc_g)**0.2*Cp2Imp(Cp_g)/(Pr_g**0.6))*
               (kg2lb(m_dot)/m22in2(A_t))**0.8*
          (D_t/(0.941*R_t))**0.1    # Use average throat curviture
          )*(A_t/A)**0.9 * sigma)
    hg_list[i] = hg
    # Heat transfer equations
    eq1 = q - hg*(Taw - Twg)
    eq2 = q - k_316*(Twg-Twc)/shell_thick_i
    eq3 = q - fin_factor(x,hc)*hc*(Twc - Tco)
    return [eq1, eq2, eq3]

# Initial conditions
initial_guess = np.array([3.5e6,1100,700])
dx = abs(x_list[2]-x_list[1])

# Array initializtion
dA_list = np.zeros(x_list.shape)
q_list = np.zeros(x_list.shape)
Twg_list = np.zeros(x_list.shape)
Twc_list = np.zeros(x_list.shape)
Tco_list = np.zeros(x_list.shape)
Dco_hyd_list = np.zeros(x_list.shape)
Vco_list = np.zeros(x_list.shape)
deltaP_list = np.zeros(x_list.shape)
stress_list = np.zeros(x_list.shape)
hg_list = np.zeros(x_list.shape)

#solution = fsolve(heat_trans_eqs,initial_guess, args=(x_list[0],Tco,))
#print(f"q: {solution[0]}\tTwg: {solution[1]}\tTwc: {solution[2]}")

# Loop through each point starting from the nozzle and workin back to injector
for i, x in enumerate(x_list):
    # Calculate change in area
    epsilon = 0.000001
    dRdx = (R_engine(x + epsilon) - R_engine(x)) / epsilon
    L = dx * sqrt(1 + dRdx ** 2)
    dA = 2 * pi * (R_engine(x) + shell_thick_i) * L
    dA_list[i] = dA

    # Solve heat transfer steady state, store results
    solution = fsolve(heat_trans_eqs,initial_guess, args=(x,Tco,i,))
    Tco_list[i] = Tco
    q_list[i] = solution[0]
    Twg_list[i] = solution[1]
    Twc_list[i] = solution[2]

    # Update values for next loop
    initial_guess = np.array([q_list[i], Twg_list[i], Twc_list[i]])
    Tco += q_list[i] * dA / (m_dot_fuel * Cp_ethanol(Tco))

    # PRESSURE DROPS
    #------------------------------------
    # Find the hydraulic diameter of the cooling channel
    area_channel = theta_chan / 2 * (
                (R_engine(x) + shell_thick_i + channel_h) ** 2 - (R_engine(x) + shell_thick_i) ** 2) - channel_h * channel_w  # m^2
    perim_channel = 2 * channel_h + theta_chan * (R_engine(x) + shell_thick_i + channel_h) + theta_chan * (
                R_engine(x) + shell_thick_i) - 2 * channel_w  # m
    d_hyd = 4 * area_channel / perim_channel  # m     Hydraulic diameter

    # Determine physical properties of coolant
    visc_co = visc_ethanol(Tco)  # Pa-s   Coolant bulk Viscosity
    Cp_co = Cp_ethanol(Tco)  # J/kg-K Coolant Isobaric heat capacity
    Pr_co = Cp_co * visc_co / k_co  # N/A    Coolant Prandtl's number
    V_co = 2*m_dot_fuel / n_channels / rho_ethanol(Tco)/ area_channel  # m/s    Coolant channel velocity #rho_ethanol(Tco)
    Re_co = rho_ethanol(Tco) * V_co * d_hyd / visc_co  # N/A    Coolant Reynold's number
    Re_w = rho_ethanol(Twc_list[i]) * V_co * d_hyd / visc_ethanol(Twc_list[i])

    Dco_hyd_list[i] = d_hyd
    deltaP_list[i] = friction_factor(d_hyd,Re_co,Re_w,Twc_list[i],Tco) * (L/d_hyd) * rho_ethanol(Tco)/2 * V_co**2
    stress_list[i] = (pc_ns*1.1-p_list[i])*R_engine(x)/shell_thick_i + E_316*cte_316*q_list[i]*shell_thick_i/(2*(1-poisson_316)*k_316)
    Vco_list[i] = V_co
    print(hg_list[i])
    #print(Twg_list[i]-Twc_list[i])
    #print(Pa2MPa(stress_list[i]))
    #print

#print(np.sum(deltaP_list))
print(Pa2psia(np.sum(deltaP_list)))
#print(np.sum(dA_list))
T_boil = 200+273
T_i = 293
T_avg = (T_boil+T_i)/2
H_vap = 900000
eta = 0.5

indices = (x_list>=x_3)
A_cool = np.sum(dA_list[indices])
m_film = A_cool/eta/(Cp_ethanol(T_avg)*(T_boil-T_i)/(7000*(Tc_ns-T_avg)) + H_vap/(7000*(Tc_ns-T_boil)))
print(f"Film cooling mass flow: {m_film}")
print(f"Fuel mass flow {m_dot_fuel}")

plt.plot(x_list, hg_list/1000)
plt.ylabel("Stress (MPa)")
plt.show()

mass_frac = lambda m_film: (Tc_ns-(800+273))/(Tc_ns-450) - np.exp(-1400/(m_film/.012*7500*.25))
print("Coolant flow rate:")
print(fsolve(mass_frac,0.07)[0])
# Create a figure and axis
fig, ax1 = plt.subplots()

# Plot q_list vs x_list on the left axis
ax1.plot(x_list[:-1], q_list[:-1]/2943623, color='blue', label='Heat Flux')
ax1.set_xlabel('Nozzle Position (m)')
ax1.set_ylabel('Heat Flux (W/m^2)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

# Create a second y-axis for Twc_list
ax2 = ax1.twinx()
ax2.plot(x_list[:-1], Twc_list[:-1], color='red', label='Coolant Side Wall Temp')
ax2.set_ylabel('Coolant Side Wall Temp (K)', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.set_ylim(0, max(Taw_list)*1.05)

# Create a third y-axis for Twg_list
ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 50))  # Offset the third axis to the right
ax3.plot(x_list[:-1], Twg_list[:-1], color='green', label='Gas Side Wall Temp')
ax3.set_ylabel('Gas Side Wall Temp (K)', color='green')
ax3.tick_params(axis='y', labelcolor='green')
ax3.set_ylim(0, max(Taw_list)*1.05)

# Create a fourth y-axis for Tco_list
ax4 = ax1.twinx()
ax4.spines['right'].set_position(('outward', 100))  # Further offset the fourth axis
ax4.plot(x_list[:-1], Tco_list[:-1], color='orange', label='Coolant Temp')
ax4.set_ylabel('Coolant Temp (K)', color='orange')
ax4.tick_params(axis='y', labelcolor='orange')

# Create a fourth y-axis for Tco_list
ax5 = ax1.twinx()
ax5.spines['right'].set_position(('outward', 150))  # Further offset the fourth axis
ax5.plot(x_list[:-1], Taw_list[:-1], color='purple', label='Adiabatic Wall')
ax5.set_ylabel('Adiabatic Temp (K)', color='purple')
ax5.tick_params(axis='y', labelcolor='purple')
ax5.set_ylim(0, max(Taw_list)*1.05)

# Show the plot
fig.tight_layout()  # To ensure everything fits
plt.title("Without Fin Factor")
plt.show()