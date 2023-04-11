# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 14:24:18 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import zipfile
import pandas as pd
from scipy import interpolate

"""Check ""Compare CLD Profiles"" file exists, if not create it"""
path = "Compare CLD Profiles"
isExist = os.path.isdir(path)
if not isExist:
    os.makedirs(path)
    print("Path Created")

"""Check ""CLD Profiles"" file exists, if not create it"""
path2 = "CLD Profiles"
isExist = os.path.isdir(path2)
if not isExist:
    os.makedirs(path2)
    print("Path2 Created")

zf = zipfile.ZipFile("Aerofoil-data.zip")

"""Plotting the 00-50 CLD profiles"""
for i in range(51):
    # Set the file number containing aerofoil data
    number = i

    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"
    name = "Profile-" + filenumb + "-CLD.csv"

    df = pd.read_csv(zf.open(name), header=None)
    df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
    df = {0: pd.to_numeric(df[0]),
          1: pd.to_numeric(df[1]),
          2: pd.to_numeric(df[2])}
    df = pd.DataFrame(df)

    # Create a function to interpolate values for Cl and Cd
    fcl = interpolate.interp1d(df[0], df[1])
    fcd = interpolate.interp1d(df[0], df[2])

    # Plot the comparison CLD profiles on separate figures (interpolation)
    plt.figure(1, figsize=(24, 12))
    plt.suptitle(filenumb, fontsize=40)
    plt.subplots_adjust(wspace=0.25)
    plt.subplot(121)
    plt.scatter(df[0], df[1])
    plt.scatter(df[0], df[2])
    plt.title("No Interpolation", fontsize=30)
    plt.xlim(-180, 180)
    plt.ylim(-1.5, 2.5)
    plt.legend(labels=["Lift", "Drag"], fontsize=20)
    plt.axhline(y=0, color="black", linestyle="--")
    plt.subplot(122)
    plt.plot(np.linspace(-180, 180, 1000), fcl(np.linspace(-180, 180, 1000)))
    plt.plot(np.linspace(-180, 180, 1000), fcd(np.linspace(-180, 180, 1000)))
    plt.title("Interpolated", fontsize=30)
    plt.xlim(-180, 180)
    plt.ylim(-1.5, 2.5)
    plt.legend(labels=["Lift", "Drag"], fontsize=20)
    plt.axhline(y=0, color="black", linestyle="--")
    plt.savefig(path + f"/{number:02d}")
    plt.show()

    # Plot the CLD profiles on separate figures (csv points)
    plt.figure(1, figsize=(12, 12))
    plt.scatter(df[0], df[1])
    plt.scatter(df[0], df[2])
    plt.title(filenumb, fontsize=30)
    plt.xlim(-180, 180)
    plt.ylim(-1.5, 2.5)
    plt.legend(labels=["Lift", "Drag"], fontsize=20)
    plt.axhline(y=0, color="black", linestyle="--")
    plt.savefig(path2 + f"/{number:02d}")
    plt.show()


"""Plotting the NACA Airfoil CLD profile"""
# Store CSV data as pandas dataframe
name = "NACA 63-415 AIRFOIL Aerodynamic Data.csv"
df = pd.read_csv(zf.open(name), header=None)
df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
df = {0: pd.to_numeric(df[0]),
      1: pd.to_numeric(df[1]),
      2: pd.to_numeric(df[2])}
df = pd.DataFrame(df)

# Create a function to interpolate values for Cl and Cd
fcl = interpolate.interp1d(df[0], df[1])
fcd = interpolate.interp1d(df[0], df[2])

# Plot the comparison CLD profiles on separate figures (interpolation)
plt.figure(1, figsize=(24, 12))
plt.suptitle("NACA 63-415 AIRFOIL", fontsize=40)
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.scatter(df[0], df[1])
plt.scatter(df[0], df[2])
plt.title("No Interpolation", fontsize=30)
plt.xlim(-15.5, 19.25)
plt.ylim(-1.5, 2.5)
plt.legend(labels=["Lift", "Drag"], fontsize=20)
plt.axhline(y=0, color="black", linestyle="--")
plt.subplot(122)
plt.plot(np.linspace(-15.5, 19.25, 1000), fcl(np.linspace(-15.5, 19.25, 1000)))
plt.plot(np.linspace(-15.5, 19.25, 1000), fcd(np.linspace(-15.5, 19.25, 1000)))
plt.title("Interpolated", fontsize=30)
plt.xlim(-15.5, 19.25)
plt.ylim(-1.5, 2.5)
plt.legend(labels=["Lift", "Drag"], fontsize=20)
plt.axhline(y=0, color="black", linestyle="--")
plt.savefig(path + "/NACA 63-415 AIRFOIL")
plt.show()

# Plot the CLD profiles on separate figures (csv points)
plt.figure(1, figsize=(12, 12))
plt.scatter(df[0], df[1])
plt.scatter(df[0], df[2])
plt.title("NACA 63-415 AIRFOIL", fontsize=30)
plt.xlim(-15.5, 19.25)
plt.ylim(-1.5, 2.5)
plt.legend(labels=["Lift", "Drag"], fontsize=20)
plt.axhline(y=0, color="black", linestyle="--")
plt.savefig(path2 + "/NACA 63-415 AIRFOIL")
plt.show()


"""Plotting the Lecture CLD Profile"""
name = "CLD.csv"

df = pd.read_csv(zf.open(name), header=None)
df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
df = {0: pd.to_numeric(df[0]),
      1: pd.to_numeric(df[1]),
      2: pd.to_numeric(df[2])}
df = pd.DataFrame(df)

# Create a function to interpolate values for Cl and Cd
fcl = interpolate.interp1d(df[0], df[1])
fcd = interpolate.interp1d(df[0], df[2])

# Plot the comparison CLD profiles on separate figures (interpolation)
plt.figure(1, figsize=(24, 12))
plt.suptitle("CLD", fontsize=40)
plt.subplots_adjust(wspace=0.25)
plt.subplot(121)
plt.scatter(df[0], df[1])
plt.scatter(df[0], df[2])
plt.title("No Interpolation", fontsize=30)
plt.xlim(-180, 180)
plt.ylim(-1.5, 2.5)
plt.legend(labels=["Lift", "Drag"], fontsize=20)
plt.axhline(y=0, color="black", linestyle="--")
plt.subplot(122)
plt.plot(np.linspace(-180, 180, 1000), fcl(np.linspace(-180, 180, 1000)))
plt.plot(np.linspace(-180, 180, 1000), fcd(np.linspace(-180, 180, 1000)))
plt.title("Interpolated", fontsize=30)
plt.xlim(-180, 180)
plt.ylim(-1.5, 2.5)
plt.legend(labels=["Lift", "Drag"], fontsize=20)
plt.axhline(y=0, color="black", linestyle="--")
plt.savefig(path + "/CLD")
plt.show()

# Plot the CLD profiles on separate figures (csv points)
plt.figure(1, figsize=(12, 12))
plt.scatter(df[0], df[1])
plt.scatter(df[0], df[2])
plt.title("CLD", fontsize=30)
plt.xlim(-180, 180)
plt.ylim(-1.5, 2.5)
plt.legend(labels=["Lift", "Drag"], fontsize=20)
plt.axhline(y=0, color="black", linestyle="--")
plt.savefig(path2 + "/CLD")
plt.show()

# """Plotting the profiles overlayed"""
# # Initialise the figure for the overlay plot
# plt.figure(1, figsize=(6, 6))

# for i in [0, 0, 5, 46, 19, 27, 10, 24, 28, 22, 20, 39, 43, 37, 45, 44, 35]:
#     # Set the file number containing aerofoil data
#     number = i

#     # Force the number to be 0x for 0-9
#     filenumb = f"{number:02d}"
#     name = "Profile-" + filenumb + "-Geom.csv"

#     # Store CSV data as pandas dataframe
#     df = pd.read_csv(zf.open(name), header=None)
#     df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
#     df = {0: pd.to_numeric(df[0]), 1: pd.to_numeric(df[1])}
#     df = pd.DataFrame(df)

#     # Plot the blade profiles overlayed on eachother
#     plt.plot(df[0], df[1])
#     plt.xlim(-0.1, 1.1)
#     plt.ylim(-0.6, 0.6)
#     plt.title("Overlay Plot of Profiles")
#     plt.axhline(y=0, color="black", linestyle="--")

# plt.savefig(path + "/Overlay")
# plt.show()
