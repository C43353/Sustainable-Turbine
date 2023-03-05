# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:41:32 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate

number = 11  # Set the file number containing aerofoil data
filenumb = f"{number:02d}"  # Force the number to be 0x for 0-9

# Constants (currently)
R = 80  # Radius (m)
r = 40.5  # Radial Position (m) (will have to iterate over)
c = 3.256  # Aerofoil Chord Length (m) (Depends on radial position)
B = 3  # Number of Blades
omega = 2.83  # Angular Veolcity (rad/s)
theta = 4.188  # Pitch Angle (degree)
V0 = 10  # Wind Speed (m/s)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

xi = (omega * r) / V0  # Local Velocity Ratio
s = (c * B) / (2 * np.pi * r)  # Solidity

aa = 0.0  # Induction Factor
ar = 0.0  # Angular Induction Factor

# Open CSV containing aerofoil CLD profile
file = open('C:\\Users\\C43353\\OneDrive - University of Southampton\\Year 3\\'
            'Technology Fundamentals for Sustainable Energy\\Group Coursework'
            '\\Code\\Aerofoil-data\\Profile-' + filenumb + '-CLD.csv')
data = file.read()
file.close()

# Split CSV into lines and convert to numeric
lines = data.split('\n')
info = []
for line in lines[0:200]:
    info.append([float(i) for i in line.split(',')])

# Convert aerofoil info into np.array to allow easier calcuations
info = np.array(info)

# Create a function to interpolate values for Cl and Cd
fcl = interpolate.interp1d(info[:, 0], info[:, 1])
fcd = interpolate.interp1d(info[:, 0], info[:, 2])

# Create plots to check function fits data
plt.figure(1, figsize=(18, 6))
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.plot(info[:, 0], info[:, 1])
plt.title("CSV Cl")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcl(np.linspace(-180, 180, 1000)))
plt.title("Func Cl")
plt.show()

plt.figure(1, figsize=(18, 6))
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.plot(np.linspace(-180, 180, 1000), fcd(np.linspace(-180, 180, 1000)))
plt.title("CSV Cd")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcd(np.linspace(-180, 180, 1000)))
plt.title("Func Cd")
plt.show()

for i in range(10):

    phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

    alpha = phi * (180 / np.pi) - theta  # Angle of Attack

    Cl = fcl(alpha)
    Cd = fcd(alpha)

    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
    Cr = Cl * np.sin(phi) - Cd * np.cos(phi)

    F = (2 / np.pi) * np.arccos(np.exp(- (B * (1 - (r / R))) / (2 * (r / R) * np.sin(phi) * r)))

    K = (4 * F * (np.sin(phi) ** 2)) / (s * Cn)

    print(aa, ar, phi, alpha, Cl, Cd, Cn, Cr)

    aa = 1 / (K + 1)
    # From Lecture Calculation
    if aa > ac:
        print("over")
        aa = 1 - ((K * (1 - (2 * ac))) / 2) * (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

    # From Converting Matlab Code
    # if aa > ac:
    #     aa = 1 + np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2))
    #     aa = 1 - (K * ((1 - (2 * ac)) / (2)))

    ar = 1 / ((4 * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)
