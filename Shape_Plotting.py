# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:44:19 2023

@author: C43353
"""

import matplotlib.pyplot as plt
import os
import zipfile
import pandas as pd
import numpy as np

# Change default dpi for plots (used for png (svg is vector based))
# plt.rcParams['figure.dpi'] = 600

# Change default saved figure format to svg
# (smaller file size than high resolution png but better quality)
plt.rcParams['savefig.format'] = "svg"

"""Check ""Shape Profiles"" file exists, if not create it"""
path = "Shape Profiles"
isExist = os.path.isdir(path)
if not isExist:
    os.makedirs(path)
    print("Path Created")

zf = zipfile.ZipFile("Aerofoil-data.zip")

"""Plotting the 00-50 profiles"""
for i in range(51):
    # Set the file number containing aerofoil data
    number = i

    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"
    name = "Profile-" + filenumb + "-Geom.csv"

    # Store CSV data as pandas dataframe
    df = pd.read_csv(zf.open(name), header=None)
    df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
    df = {0: pd.to_numeric(df[0]), 1: pd.to_numeric(df[1])}
    df = pd.DataFrame(df)

    # Plot the blade profiles on separate figures
    plt.figure(1, figsize=(12, 12))
    plt.plot(df[0], df[1])
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.6, 0.6)
    plt.title(f"{number:02d}", fontsize=30)
    plt.axhline(y=0, color="black", linestyle="--")
    plt.savefig(path + f"/{number:02d}")
    plt.show()


"""Plotting the NACA Airfoil profile"""
# Store CSV data as pandas dataframe
name = "NACA 63-415 AIRFOIL Geometry.csv"
df = pd.read_csv(zf.open(name), header=None)
df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
df = {0: pd.to_numeric(df[0]), 1: pd.to_numeric(df[1])}
df = pd.DataFrame(df)

# Plot the blade profiles on separate figures
plt.figure(1, figsize=(12, 12))
plt.plot(df[0], df[1])
plt.ylim(-0.6, 0.6)
plt.title("NACA 63-415 AIRFOIL", fontsize=30)
plt.axhline(y=0, color="black", linestyle="--")
plt.savefig(path + "/NACA 63-415 AIRFOIL")
plt.show()


"""Plotting the profiles overlayed"""
# Initialise the figure for the overlay plot
plt.figure(1, figsize=(12, 12))
colour = iter(plt.cm.rainbow(np.linspace(0, 1, 17)))
for i in [42, 46, 4, 24, 12, 36, 26, 48, 13, 34, 33, 3, 11, 45, 15, 23, 38]:
    # Set the file number containing aerofoil data
    number = i

    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"
    name = "Profile-" + filenumb + "-Geom.csv"

    # Store CSV data as pandas dataframe
    df = pd.read_csv(zf.open(name), header=None)
    df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
    df = {0: pd.to_numeric(df[0]), 1: pd.to_numeric(df[1])}
    df = pd.DataFrame(df)

    c = next(colour)
    # Plot the blade profiles overlayed on eachother
    plt.plot(df[0], df[1], c=c)
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.6, 0.6)
    plt.title("Overlay Plot of Profiles", fontsize=30)
    plt.axhline(y=0, color="black", linestyle="--")

plt.savefig(path + "/Overlay")
plt.show()
