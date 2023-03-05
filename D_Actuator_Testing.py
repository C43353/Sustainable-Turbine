# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 11:59:02 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt

P = 0.5E6  # Power (W)
rho = 1.225  # Air Density (kg/m^3)
v0 = np.linspace(5, 20, 1000)  # Wind Speed (m/s)
a = np.linspace(0.1, 0.6, 1000)  # Induction Factor

# Ct = 4 * a * (1 - a)  # Thrust Coefficient
# Cp = 4 * a * ((1 - a) ** 2)  # Power Coefficient

# A = P / (Cp * rho * (v0 ** 3))  # Swept Area of Turbine

# r = np.sqrt(A / np.pi)  # Radius of Turbine


def radius(P, rho, v0, a):

    Cp = 4 * a * ((1 - a) ** 2)  # Power Coefficient

    # Combine variables into a matrix to allow contour plotting
    [X, Y] = np.meshgrid(v0, Cp)

    r = np.sqrt((P / (Y * rho * (X ** 3))) / np.pi)  # Radius of Turbine

    return r, Cp, v0


r, Cp, v0 = radius(P, rho, v0, a)

plt.contourf(v0, a, r, 1000)

plt.colorbar(label="Turbine Radius / m")

c = plt.contour(v0, a, r, [10], colors="k")
plt.clabel(c, fmt="Radius = 10 m", inline=1, manual=[(16, 0.3)])

plt.xlabel("Wind Speed / ms$^-$$^1$")
plt.xlim(min(v0), max(v0))

plt.ylabel("Induction Coefficient")
plt.ylim(min(a), max(a))

# plt.savefig("test.png", dpi=1000)
