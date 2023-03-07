# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 19:06:46 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
import BEM_Functions as BEM

V0 = 10  # Wind Speed (m/s)
B = 3  # Number of Blades
omega = 2.83  # Angular Velocity (rad/s)
theta = 4.85  # Pitch Angle (degrees)

r = 10.5  # Radial Position (m) (Lecture uses 10.5 and 40.5)
c = 1.29  # Chord Length (m) (Lecture uses 1.29 and 3.256)

xi = (omega * r) / V0  # Local Velocity Ratio
s = (c * B) / (2 * np.pi * r)  # Solidity

aa = 0  # Initial Induction Factor
ar = 0  # Initial Radial Induction Factor

aa_list = []
ar_list = []

iterations = 10
# Iterate to a constant value (1D Actuator theory also invalid for a > 0.5)
for i in range(iterations):
    phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle (radians)

    alpha = phi * (180 / np.pi) - theta  # Angle of Attack (degrees)

    fcl, fcd = BEM.cld_func('Aerofoil-data\\CLD.csv')  # Import Cl, Cd functions

    Cl = fcl(alpha)  # Lift Coefficient
    Cd = fcd(alpha)  # Drag Coefficient

    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Force Coefficient
    Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangential Force Coefficient

    # Final Induction Factors
    aa = 1 / (((4 * np.sin(phi) * np.sin(phi)) / (s * Cn)) + 1)
    ar = 1 / (((4 * np.sin(phi) * np.cos(phi)) / (s * Cr)) - 1)

    aa_list.append(aa)
    ar_list.append(ar)

angle = np.arctan((ar * xi) / aa)  # Angle of Induction (if != phi,
                                   # justifies using BEM)

print(f"phi = {round(phi, 2)}, Angle of Induction Flow = {round(angle, 2)}")

if angle != phi:
    print("Angle of Induction != phi, therefore using BEM is justified")

plt.plot(range(iterations), aa_list)
plt.show()
plt.plot(range(iterations), ar_list)
