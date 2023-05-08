# -*- coding: utf-8 -*-
"""
Created on Mon May  8 21:30:37 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import zipfile
from scipy import interpolate

R = 85
r = 20
V0 = 10
profile = 20
B = 3  # Number of Blades

# omega = 2.83  # Angular Veolcity (rad/s) (Constant for varying wind speeds)
tsr = 7  # Tip Speed Ratio (Used to define the angular velocity)
omega = (tsr * 10) / R  # Angular Velocity (dependent on tip speed ratio)
rpm = (omega * 60) / (2 * np.pi)

method = 3

""" Rank Airfoils, Find Optimum Angle of Attack """
# Create dictionaries to store data
data = {}
maxcld1 = {}

# Define the zipfile that stores all the cld data CSVs
zf = zipfile.ZipFile("Aerofoil-data.zip")

for number in range(51):
    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"
    name2 = "Profile-" + filenumb + "-CLD.csv"
    df2 = pd.read_csv(zf.open(name2), header=None)
    df2 = df2[pd.to_numeric(df2[0], errors="coerce").notnull()]
    df2 = {0: pd.to_numeric(df2[0]),
           1: pd.to_numeric(df2[1]),
           2: pd.to_numeric(df2[2])}
    data[number] = pd.DataFrame(df2)

    data[number][3] = data[number][1] / data[number][2]

    maxcld1[number] = data[number].loc[data[number][3].idxmax()]

maxcld = {key: val for key, val in sorted(maxcld1.items(),
                                          key=lambda ele: ele[1][3])}

alpha = list(maxcld1.items())[profile][1][0]

fcl = interpolate.interp1d(data[profile][0], data[profile][1])
fcd = interpolate.interp1d(data[profile][0], data[profile][2])


ac = 1/3  # Critical Induction Factor

aa = 0.0  # Induction Factor
ar = 0.0  # Angular Induction Factor

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
c_list = []
theta_list = []


xi = (omega * r) / V0  # Local Velocity Ratio
tsr = (omega * R) / V0  # Tip Speed Ratio

for i in range(20):
    phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

    theta = phi * (180 / np.pi) - alpha  # Twist Angle

    Cl = float(fcl(alpha))  # Lift Coefficient
    Cd = float(fcd(alpha))  # Drag Coefficient

    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Coefficient
    Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangent Coefficient

    if method == 1:
        # https://www.ehow.co.uk/how_7697179_calculate-along-wind-turbine-blade.html
        c = (5.6 * R**2) / (B * Cl * r * (tsr**2))

    elif method == 2:
        # https://www.mdpi.com/1996-1073/13/9/2320
        c = (8 * np.pi * r * np.sin(phi)) / (3 * B * Cl * tsr)

    elif method == 3:
        # https://ieeexplore.ieee.org/abstract/document/7884538
        c = ((8 * np.pi * r) / (B * Cl)) * (1 - np.cos(phi))

    else:
        print("Must select a suitable method of calculating chord length"
              "1, 2 or 3, default=1")

    s = (c * B) / (2 * np.pi * r)  # Solidity

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
    K_list.append(K)
    c_list.append(c)
    theta_list.append(theta)

    # Calc New Induction Factor Using Calculation Given in Slides
    aa = 1 / (K + 1)
    if aa > ac:
        aa = 1 - ((K * (1 - (2 * ac))) / 2) * \
            (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

    ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

plt.plot(aa_list)
plt.show()
plt.plot(ar_list)
