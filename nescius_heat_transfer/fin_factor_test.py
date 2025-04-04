from math import sqrt, exp

import matplotlib.pyplot as plt
import numpy as np
from mpmath.rational import mpq_1

s = 1.5
w = 1
L1 = 3
L2 = s*1.1
b = w*2
h = 5
k = 15
T_wg = 700
T_co = 400
theta = T_wg-T_co

m1 = sqrt(2*h/(k*w))
m2 = sqrt(2*h/(k*b))



c_1 =(
    theta*exp(2 * L1 * m1) * (b * m2 * exp(2*L2 * m2) - b * m2 + m1 * w * exp(2 * L2 * m2) + m1 * w)
) / (
    -b * m2 * exp(2 * L1 * m1) - b * m2 * exp(2*L2 * m2) + b * m2 * exp(2 * L1 * m1 + 2*L2 * m2)
    + b * m2 + m1 * w * exp(2 * L1 * m1) + m1 * w * exp(2 * L2 * m2)
    + m1 * w * exp(2 * L1 * m1 + 2 * L2 * m2) + m1 * w
)
c_2 = (
    theta*(-b * m2 * exp(2*L2 * m2) + b * m2 + m1 * w * exp(2 * L2 * m2) + m1 * w)
) / (
    -b * m2 * exp(2 * L1 * m1) - b * m2 * exp(2*L2 * m2) + b * m2 * exp(2 * L1 * m1 + 2*L2 * m2)
    + b * m2 + m1 * w * exp(2 * L1 * m1) + m1 * w * exp(2 * L2 * m2)
    + m1 * w * exp(2 * L1 * m1 + 2 * L2 * m2) + m1 * w
)

eta_f = 2*s/(2*s+w) + k*w*(c_1*m1-c_2*m2)/(h*(2*s+w)*theta)
print(eta_f)

c_4 = (c_1*exp(-m1*L1)+c_2*exp(m1*L1))/(1+exp(2*m2*L2))
c_3 = c_4*exp(2*m2*L2)

res1 = theta-c_1-c_2
res2 = w*m1*(-c_1*exp(-m1*L1)+c_2*exp(m1*L1))-m2*b*(-c_4*exp(2*m2*L2)+c_4)
res3 = (c_1*exp(-m1*L1)+c_2*exp(m1*L1))/(1+exp(2*m2*L2))-c_4

x1 = np.linspace(0,L1,10000)
x2 = np.linspace(L1,L1+L2,10000)
T1 = c_1 * np.exp(-m1*x1) + c_2 * np.exp(m1*x1)+T_co
T2 = c_3 * np.exp(-m2*(x2-L1)) + c_4 * np.exp(m2*(x2-L1))+T_co

q1 =-k*(T1[-1]-T1[-2])/(x1[-1]-x1[-2]) * w
q2 =-k*(T2[1]-T2[0])/(x2[1]-x2[0]) * b
q2_end =-k*(T2[-1]-T2[-2])/(x2[-1]-x2[-2]) * b
print(q1-q2)
print(q2_end)


plt.plot(np.concatenate((x1,x2)),np.concatenate((T1,T2)))
plt.xlabel("Position along fin (mm)")
plt.ylabel("Temperature (K)")
plt.show()



