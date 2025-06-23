from math import floor
from math import exp, sqrt, log, tanh
from math import pi

visc_v = 1.0405e-4         # Pa-s        | Gas viscosity
visc_ethanol = lambda T: (0.00201 * exp(1614/T + 0.00618*T - 1.132E-5*T**2))/1000   # Pa-s
visc_l = visc_ethanol(50+273)

Gamma_cr = 1.01e5 * visc_v**2/visc_l
D_c = .0452
Circ = D_c * pi
e_t = .18

m_cr = Gamma_cr*Circ

m_dot = 0.12607 # kg/s
m_film = m_dot * .03
pho_l = 789
pho_g = 2.402 # kg/m^3
lambduh = 22500/0.04607
T_v = 473
T_i = 100+273
T_g = 3224
dT = T_g-T_v
Cpl = 4000  # J/kg-K
Cpg = 2000  # J/kg-K
Pr_g = 0.52  
lambda_star = lambduh + Cpl*(T_v-T_i)

K_t = 1 + 4*e_t

def Xe(x):
    m = 1.2
    return 3.53*D_c*(1+(x/(3.53*D_c))**-m)**(-1/m)

Gamma = m_film/Circ

U_g = m_dot/(pho_g*D_c**2*pi/4)
G_ch = pho_g*U_g
Rex = pho_g*U_g*.005/visc_v
Cf = .0592 * Rex**-.2

t_guess = .000045
U_liq = Gamma/(pho_l*t_guess)
U_l = 2*U_liq
print(U_l)
G_mean = G_ch*(T_g/((T_v+T_g)/2))*(U_g-U_l)/U_g

T_w = 1/2 *Cf*G_mean*(U_g - U_l)
t = sqrt(2*visc_l*Gamma/(pho_l*T_w))
print(f"Film thickness: (mm) {t*1000}")

Re_c = pho_l*U_l*D_c/visc_l
print(U_l)
Re_cfilm = 250*log(Re_c)-1265
E_m = 1 - Re_cfilm/Re_c
a = 2.31e-4 *Re_c**-.35
d_pho = pho_l-pho_g
sigma = 0.022
St0 = 1/2*Cf*Pr_g**.8

We = pho_g*U_g**2*D_c/(sigma*(d_pho/pho_g)**.25)

Lc = Gamma*lambda_star/(G_mean*Cpg*dT*St0)
E = E_m * tanh(a*We**1.25)
print(E)
print(a*We**1.25)
