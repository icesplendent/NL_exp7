from scipy.optimize import fsolve
import numpy as np
import math
import matplotlib.pyplot as plt

A_1 = 5.35583
A_2 = 0.100473
A_3 = 0.20692
A_4 = 100.0
A_5 = 11.34927
A_6 = 1.5334e-2
B_1 = 4.629e-7
B_2 = 3.862e-8
B_3 = -0.89e-8
B_4 = 2.657e-5
T_0 = 24.5
lam_p = 1.064
lam_s = 1.428067
lam_i = 4.173585
I_p = 10**12
epsilon = 8.854 * 10**(-12)
c = 3 * 10**8
h = 6.626 * 10**(-34)
L = 2.5 * 10**(-2)
deff_test = 14 * 10**(-12)

def F(T):
    return (T - T_0) * (T + 570.82)

def n_e(lam, T):
    return np.sqrt(A_1 + B_1 * F(T) + 
                   (A_2 + B_2 * F(T)) / (lam ** 2 - (A_3 + B_3 * F(T)) ** 2 ) +
                   (A_4 + B_4 * F(T)) / (lam ** 2 - A_5 ** 2) - (A_6 * lam ** 2))

def omega(lam):
    return 2 * np.pi * c / (lam * 10 **(-6))  # unit of wavelength in um

def Gamma(deff):
    return np.sqrt(2 * omega(lam_s) * omega(lam_i) * deff**2 * I_p / 
                   (n_e(lam_s, 51) * n_e(lam_i, 51) * n_e(lam_p, 51) * epsilon *  c**3)
                   )

def vacuum_energy(lam):
    return h * c / (lam * 10 **(-6))

# generated signal energy after crystal length
def signal_energy():
    return vacuum_energy(lam_s) * 0.25 * np.exp(2 * Gamma(deff_test) * L)


# print(n_e(lam_p, 51))
# print(n_e(lam_s, 51))
# print(n_e(lam_i, 51))
print(vacuum_energy(lam_s))
print(Gamma(deff_test))
print(signal_energy())

print(signal_energy() * lam_s / lam_p)