# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 19:30:18 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
import BEM_Functions as BEM

V0 = 10  # Wind Speed (m/s)
R = 20.94
r = 0.3103  # Radial Position (m)
c = 1.29  # Chord Length (m)
B = 3  # Number of Blades
omega = 2.83  # Angular Velocity (rad/s)
theta = 4.85  # Pitch Angle (degrees)
ac = 1/3

xi = (omega * r) / V0  # Local Velocity Ratio
s = (c * B) / (2 * np.pi * r)  # Solidity

aa = 0  # Initial Induction Factor
ar = 0  # Initial Radial Induction Factor

# Initialise the lists for the variables
phi_list = []
alpha_list = []
Cl_list = []
Cd_list = []
Cn_list = []
Cr_list = []
F_list = []
aa_list = []
ar_list = []

iterations = 20
# Iterate to a constant value (1D Actuator theory also invalid for a > 0.5)
for i in range(iterations):
    phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle (radians)

    alpha = phi * (180 / np.pi) - theta  # Angle of Attack (degrees)

    fcl, fcd = BEM.cld_func('Aerofoil-data\\CLD.csv')  # Import Cl, Cd functions

    Cl = fcl(alpha)  # Lift Coefficient
    Cd = fcd(alpha)  # Drag Coefficient

    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Force Coefficient
    Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangential Force Coefficient

    # Prandtl Loss Factor
    F = (2 / np.pi) * np.arccos(np.exp(- (B * (1 - (r / R))) /
                                       (2 * (r / R) * np.sin(phi))))

    K = (4 * F * (np.sin(phi) ** 2)) / (s * Cn)  # Useful Coefficient

    phi_list.append(phi)
    alpha_list.append(alpha)
    Cl_list.append(Cl)
    Cd_list.append(Cd)
    Cn_list.append(Cn)
    Cr_list.append(Cr)
    F_list.append(F)
    aa_list.append(aa)
    ar_list.append(ar)

    # Calc New Induction Factor Using Calculation Given in Slides
    if K > (ac ** -1) - 1:
        aa = 1 / (K + 1)  # Calc New Induction Factor
    # From Lecture Calculation
    if K <= (ac ** -1) - 1:
        aa = 1 - ((K * (1 - (2 * ac))) / 2) * \
            (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

    # # Calc New Induction Factor Using Matlab Lecture Code
    # aa = 1 / (K + 1)
    # # From Converting Matlab Code
    # if aa > ac:
    #     aa = 1 - np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2))
    #     aa = 1 + ((K * (1 - (2 * ac))) / 2) * aa

    # # Calc New Induction Factor Using Converted Matlab Code
    # aa = 1 / (K + 1)
    # # From Lecture Calculation
    # if aa > ac:
    #     aa = 1 - ((K * (1 - (2 * ac))) / 2) * (np.sqrt(1 + (4 / K) * (
    #         ((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

    ar = 1 / ((4 * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

angle = np.arctan((ar * xi) / aa)  # Angle of Induction (if != phi,
                                   # justifies using BEM)

Ct = 4 * aa * (1 - aa)  # Thrust Coefficient

print(f"phi = {round(phi, 2)}, Angle of Induction Flow = {round(angle, 2)}")

if angle != phi:
    print("Angle of Induction != phi, therefore using BEM is justified")

print(f"The empirical thrust coefficient is {round(Ct, 2)}")
print(f"The tip loss factor (Prandtl loss factor) is {round(F, 2)}")

plt.plot(range(iterations), aa_list)
plt.show()
plt.plot(range(iterations), ar_list)
