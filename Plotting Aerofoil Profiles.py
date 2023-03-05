# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 23:09:09 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt

number = 0  # Set the file number containing aerofoil data
filenumb = f"{number:02d}"  # Force the number to be 0x for 0-9

file = open('Aerofoil-data\\Profile-' + filenumb + '-Geom.csv')
data = file.read()
file.close()

# Split CSV into lines and convert to numeric
lines = data.split('\n')
lines = lines[:-1]
info = []
for line in lines:
    if len(line) != 0:
        info.append([float(i) for i in line.split(',')])

info = np.array(info)

x = info[:, 0]
y = info[:, 1]

plt.plot(x, y)
plt.xlim(-0.1, 1.1)
plt.ylim(-0.6, 0.6)
plt.axhline(y=0, color="black", linestyle="--")
