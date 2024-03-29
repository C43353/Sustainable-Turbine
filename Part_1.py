# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:33:53 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions import cld_func

P = 8E6  # Desired Power Output (MW)
V0 = 10  # Nominal Wind Speed (m/s)
tsr = 7  # Nominal Tip Speed Ratio
B = 3  # Number of Turbine Blades
r = 9.46875  # Radial Position (m)
c = 1.597  # Chord Length (m)
theta = 16.3  # Pitch Angle (degrees)

rho = 1.225  # Air Density (kg/m^3)
ac = 1/3  # Critical Induction Factor


"""1D Initial Analysis to Find Minimum Radius"""
# Calculate the minimum raidus to obtain power output
R = np.sqrt((P * 27) / ((8 * np.pi * rho * (V0 ** 3))))
print(f"The Minimum Radius Required to Obtain {P*10**-6} MW Output is "
      f"{round(R, 2)} m (for a Wind Speed of {V0} m/s)")

omega = (tsr * V0) / R  # Angular Velocity (rad/s)
print(f"And if the Nominal Tip Speed Ratio is {tsr} Then the Angular "
      f"Velocity is {round(omega, 2)} rad/s")

R = 85
omega = 2.83

xi = (tsr * r) / R  # Local Velocity Ratio
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
K_list = []

# Import Cl, Cd functions
fcl, fcd = cld_func("CLD.csv")

# Iterate to a constant value (1D Actuator theory also invalid for a > 0.5)
for i in range(20):
    # Relative Wind Angle (radians)
    phi = np.arctan((1 - aa) / ((1 + ar) * xi))

    # Angle of Attack (degrees)
    alpha = phi * (180 / np.pi) - theta

    # if alpha > 19.25 or alpha < -15.5:
    #     alpha = 5

    Cl = float(fcl(alpha))  # Lift Coefficient
    Cd = float(fcd(alpha))  # Drag Coefficient

    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Force Coefficient
    Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangential Force Coefficient

    # Prandtl Loss Factor (not sure if theres a difference anymore (lecture
    # notes has an extra r that causes issues compared to example code))
    F = (2 / np.pi) * np.arccos(np.exp(- (B * (1 - (r / R))) /
                                       (2 * (r / R) * np.sin(phi))))
    # F = (2 / np.pi) * np.arccos(np.exp(- (((B / 2) * (R - r)) / r)
    #                                    / np.sin(phi)))

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
    K_list.append(K)

    # Calc New Induction Factor Using Calculation Given in Slides
    aa = 1 / (K + 1)  # Calc New Induction Factor
    # From Lecture Calculation
    if aa > ac:
        aa = 1 - ((K * (1 - (2 * ac))) / 2) * \
            (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

    ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

angle = np.arctan((ar * xi) / aa)  # Angle of Induction (if != phi,
                                   # justifies using BEM)

Ct = 4 * aa * (1 - aa)  # Thrust Coefficient

print(f"phi = {round(phi, 2)}, Angle of Induction Flow = {round(angle, 2)}")

if angle != phi:
    print("Angle of Induction != phi, Therefore Using BEM is Justified")

print(f"The Empirical Thrust Coefficient is {round(Ct, 2)}")
print(f"The Tip Loss Factor (Prandtl Loss Factor) is {round(F, 2)}")

plt.plot(aa_list)
plt.show()
plt.plot(ar_list)
