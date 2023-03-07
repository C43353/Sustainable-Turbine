# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:41:32 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from BEM_Functions import nodal

number = 11  # Set the file number containing aerofoil data
filenumb = f"{number:02d}"  # Force the number to be 0x for 0-9

# Constants (currently)
R = 80  # Radius (m)
r = 0.5  # Radial Position (m) (will have to iterate over)
c = 3.256  # Aerofoil Chord Length (m) (Depends on radial position)
B = 3  # Number of Blades
omega = 2.83  # Angular Veolcity (rad/s)
theta = 4.188  # Pitch Angle (degree)
V0 = 10  # Wind Speed (m/s)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

# Open CSV containing aerofoil CLD profile
file = open('Aerofoil-data\\Profile-' + filenumb + '-CLD.csv')
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
plt.plot(info[:, 0], info[:, 2])
plt.title("CSV Cd")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcd(np.linspace(-180, 180, 1000)))
plt.title("Func Cd")
plt.show()

phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr = nodal(R, r, V0, c, theta,
                                                      omega, B, fcl, fcd)
