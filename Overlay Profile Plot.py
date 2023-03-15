# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 17:09:52 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt

plt.figure(1, figsize=(6, 6))

for i in [0, 0, 5, 46, 19, 27, 10, 24, 28, 22, 20, 39, 43, 37, 45, 44, 35]:
    number = i  # Set the file number containing aerofoil data
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
    plt.title("Overlay Plot of Profiles")
    plt.axhline(y=0, color="black", linestyle="--")
