# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:28:31 2023

@author: C43353
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

"""
NACA airfoil only goes from -15.5-19.25 degrees
"""
"""Aerofoil-data\\NACA 63-415 AIRFOIL Aerodynamic Data.csv"""
"""Aerofoil-data\\CLD.csv"""

# Open CSV containing aerofoil CLD profile
file = open("Aerofoil-data\\CLD.csv")
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

fcl1 = interpolate.CubicSpline(info[:, 0], info[:, 1])
fcd1 = interpolate.CubicSpline(info[:, 0], info[:, 2])

fcl2 = interpolate.interp1d(info[:, 0], info[:, 1])
fcd2 = interpolate.interp1d(info[:, 0], info[:, 2])

# Create plots to check function fits data
plt.figure(1, figsize=(18, 6))
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.plot(info[:, 0], info[:, 1])
plt.title("CSV Cl")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcl1(np.linspace(-180, 180, 1000)))
plt.plot(np.linspace(-180, 180, 1000), fcl2(np.linspace(-180, 180, 1000)))
#plt.plot(np.linspace(-180, 180, 1000), np.interp(np.linspace(-180, 180, 1000), xp=info[:, 0], fp=info[:, 1]))
plt.title("Func Cl")
plt.show()

plt.figure(1, figsize=(18, 6))
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.plot(info[:, 0], info[:, 2])
plt.title("CSV Cd")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcd1(np.linspace(-180, 180, 1000)))
plt.plot(np.linspace(-180, 180, 1000), fcd2(np.linspace(-180, 180, 1000)))
#plt.plot(np.linspace(-180, 180, 1000), np.interp(np.linspace(-180, 180, 1000), xp=info[:, 0], fp=info[:, 2]))
plt.title("Func Cd")
plt.show()
