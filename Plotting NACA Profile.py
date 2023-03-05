# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 11:57:00 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt

file = open('Aerofoil-data\\NACA 63-415 AIRFOIL Geometry.csv')
data = file.read()
file.close()

# Split CSV into lines and convert to numeric
lines = data.split('\n')
lines.pop(0)
info = []
for line in lines:
    if len(line) != 0 and len(line) != 1:
        info.append([float(i) for i in line.split(',')])

info = np.array(info)

x = info[:, 0]
y = info[:, 1]

plt.plot(x, y)
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 0.1)
plt.axhline(y=0, color="black", linestyle="--")
