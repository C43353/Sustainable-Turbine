# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 15:44:30 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# Constants (currently)
R = 50  # Radius (m)
segments = np.linspace(0.6, R-0.1, 17)  # Radial Position of Segments (m)
c = 3.256  # Aerofoil Chord Length (m) (Depends on radial position)
B = 3  # Number of Blades
omega = 2.83  # Angular Veolcity (rad/s)
theta = 4.188  # Pitch Angle (degree)
V0 = 10  # Wind Speed (m/s)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

aa = 0.0  # Induction Factor
ar = 0.0  # Angular Induction Factor

# Open CSV containing aerofoil CLD profile
file = open('Aerofoil-data\\CLD.csv')
CLD = file.read()
file.close()

# Split CSV into lines and convert to numeric
lines = CLD.split('\n')
lines.pop(0)
info = []
for line in lines:
    if len(line) != 0:
        info.append([float(i) for i in line.split(',')])

# Convert aerofoil info into np.array to allow easier calcuations
info = np.array(info)

# Create a function to interpolate values for Cl and Cd
fcl = interpolate.interp1d(info[:, 0], info[:, 1])
fcd = interpolate.interp1d(info[:, 0], info[:, 2])

# final = []

phi_list = []
alpha_list = []
Cl_list = []
Cd_list = []
Cn_list = []
Cr_list = []
F_list = []
aa_list = []
ar_list = []

for r in segments:
    xi = (omega * r) / V0  # Local Velocity Ratio
    s = (c * B) / (2 * np.pi * r)  # Solidity

    for i in range(100):

        phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle
        # phi_list.append(phi)

        alpha = phi * (180 / np.pi) - theta  # Angle of Attack
        # alpha_list.append(alpha)

        Cl = float(fcl(alpha))  # Lift Coefficient
        # Cl_list.append(Cl)
        Cd = float(fcd(alpha))  # Drag Coefficient
        # Cd_list.append(Cd)

        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Coefficient
        # Cn_list.append(Cn)
        Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangent Coefficient
        # Cr_list.append(Cr)

        F = (2 / np.pi) * np.arccos(np.exp(- (B * (1 - (r / R))) / (2 * (r / R) * np.sin(phi) * r)))  # Prandtl Loss Factor
        # F_list.append(F)

        K = (4 * F * (np.sin(phi) ** 2)) / (s * Cn)  # Useful Coefficient

        # aa_list.append(aa)
        # ar_list.append(ar)

        # if aa < ac:
        aa = 1 / (K + 1)  # Calc New Induction Factor
        # From Lecture Calculation
        if aa > ac:
            aa = 1 - ((K * (1 - (2 * ac))) / 2) * (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

        # From Converting Matlab Code
        # if aa > ac:
            # aa = 1 + np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2))
            # aa = 1 - (K * ((1 - (2 * ac)) / (2)))

        ar = 1 / ((4 * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    # An array of the iterations of the code (unnecessary?)
    # outputs = np.array([phi_list, alpha_list, Cl_list, Cd_list, Cn_list, Cr_list, F_list, aa_list, ar_list])

    # Final Output Values
    phi_list.append(phi)
    alpha_list.append(alpha)
    Cl_list.append(Cl)
    Cd_list.append(Cd)
    Cn_list.append(Cn)
    Cr_list.append(Cr)
    F_list.append(F)
    aa_list.append(aa)
    ar_list.append(ar)

    # final.append([phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar])

plt.figure(1, figsize=(6, 6))
plt.plot(segments, aa_list, marker='o')
plt.plot(segments, ar_list, marker='o')
plt.title("Induction Factors")
plt.xlabel("Radial Position (m)")
plt.legend(("a", "a'"))
plt.show()

plt.figure(1, figsize=(6, 6))
plt.plot(segments, F_list, marker='o')
plt.title("Prandtl Loss Factor")
plt.xlabel("Radial Position (m)")
plt.show()

plt.figure(1, figsize=(6, 6))
plt.plot(segments, alpha_list, marker='o')
plt.plot(segments, phi_list, marker='o')
plt.title("Angles")
plt.xlabel("Radial Position (m)")
plt.legend((r"$\alpha$", r"$\phi$"))
plt.show()
