# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:18:33 2023

@author: C43353
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import zipfile
import pandas as pd

"""
NACA airfoil only goes from -15.5-19.25 degrees
"""
"""NACA 63-415 AIRFOIL Aerodynamic Data.csv"""
"""CLD.csv"""
"""Profile-00-CLD.csv"""

# Open zip file
zf = zipfile.ZipFile("Aerofoil-data.zip")

df = pd.read_csv(zf.open("Profile-30-CLD.csv"), header=None)
df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
df = {0: pd.to_numeric(df[0]),
      1: pd.to_numeric(df[1]),
      2: pd.to_numeric(df[2])}
df = pd.DataFrame(df)

fcl1 = interpolate.CubicSpline(df[0], df[1])
fcd1 = interpolate.CubicSpline(df[0], df[2])

fcl2 = interpolate.interp1d(df[0], df[1])
fcd2 = interpolate.interp1d(df[0], df[2])

# Create plots to check function fits data
plt.figure(1, figsize=(18, 6))
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.plot(df[0], df[1])
plt.title("CSV Cl")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcl1(np.linspace(-180, 180, 1000)))
plt.plot(np.linspace(-180, 180, 1000), fcl2(np.linspace(-180, 180, 1000)))
#plt.plot(np.linspace(-180, 180, 1000), np.interp(np.linspace(-180, 180, 1000), xp=df[0], fp=df[1]))
plt.title("Func Cl")
plt.show()

plt.figure(1, figsize=(18, 6))
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.plot(df[0], df[2])
plt.title("CSV Cd")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcd1(np.linspace(-180, 180, 1000)))
plt.plot(np.linspace(-180, 180, 1000), fcd2(np.linspace(-180, 180, 1000)))
#plt.plot(np.linspace(-180, 180, 1000), np.interp(np.linspace(-180, 180, 1000), xp=df[0], fp=df[2]))
plt.title("Func Cd")
plt.show()
