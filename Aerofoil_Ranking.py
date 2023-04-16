# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 18:57:35 2023

@author: C43353
"""

import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Create a dictionaries to store data
data = {}
maxcl1 = {}
maxcl = {}
maxcd = {}
maxcld = {}

zf = zipfile.ZipFile("Aerofoil-data.zip")

"""Creating dictionaries of the Aerofoil CLDs"""
for i in range(51):
    # Set the file number containing aerofoil data
    number = i

    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"
    name2 = "Profile-" + filenumb + "-CLD.csv"
    df2 = pd.read_csv(zf.open(name2), header=None)
    df2 = df2[pd.to_numeric(df2[0], errors="coerce").notnull()]
    df2 = {0: pd.to_numeric(df2[0]),
           1: pd.to_numeric(df2[1]),
           2: pd.to_numeric(df2[2])}
    df2 = pd.DataFrame(df2)

    data[number] = df2
    data[number][3] = data[number][1] / data[number][2]
    # print(filenumb)
    # print("max Cl", round(max(df2[1]), 2))
    # print("max Cd", round(max(df2[2]), 2))

    # maxcl1[number] = max(df2[1])
    maxcl[number] = df2.loc[df2[1].idxmax()]

    # maxcd[number] = max(df2[2])
    maxcd[number] = df2.loc[df2[2].idxmax()]

    # Calc lift to drag ratio
    maxcld[number] = df2.loc[df2[3].idxmax()]


# Orders dictionary in number size order
# maxcl2 = {key: val for key, val in sorted(maxcl1.items(),
#                                           key=lambda ele: ele[1])}

maxcl3 = {key: val for key, val in sorted(maxcl.items(),
                                          key=lambda ele: ele[1][1])}

maxcd = {key: val for key, val in sorted(maxcd.items(),
                                         key=lambda ele: ele[1][2])}

maxcld = {key: val for key, val in sorted(maxcld.items(),
                                          key=lambda ele: ele[1][3])}

profilescld = list(maxcld.items())[2::3]
profilesnumb = [x[0] for x in profilescld]

aoa = [x[1][0] for x in profilescld]

# Plot all CLDs overlayed
for i in range(len(data)):
    plt.plot(data[i][0], data[i][1])
plt.show()
for i in range(len(data)):
    plt.plot(data[i][0], data[i][2])
plt.show()

# Plot CLDs for selected profiles
for i in profilesnumb:
    plt.plot(data[i][0], data[i][1])
plt.show()
for i in profilesnumb:
    plt.plot(data[i][0], data[i][2])
plt.show()
