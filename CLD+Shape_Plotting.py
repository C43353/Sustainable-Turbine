# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 16:24:21 2023

@author: C43353
"""

import matplotlib.pyplot as plt
import os
import zipfile
import pandas as pd

"""Check ""Compare CLD Profiles"" file exists, if not create it"""
path = "CLD + Shape Profiles"
isExist = os.path.isdir(path)
if not isExist:
    os.makedirs(path)
    print("Path Created")

zf = zipfile.ZipFile("Aerofoil-data.zip")

"""Plotting the 00-50 CLD profiles"""
for i in range(51):
    # Set the file number containing aerofoil data
    number = i

    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"

    name1 = "Profile-" + filenumb + "-Geom.csv"
    name2 = "Profile-" + filenumb + "-CLD.csv"

    # Store aerofoil shape CSV data as pandas dataframe
    df1 = pd.read_csv(zf.open(name1), header=None)
    df1 = df1[pd.to_numeric(df1[0], errors="coerce").notnull()]
    df1 = {0: pd.to_numeric(df1[0]), 1: pd.to_numeric(df1[1])}
    df1 = pd.DataFrame(df1)

    # Store the aerofoil CLD CSV data as pandas dataframe
    df2 = pd.read_csv(zf.open(name2), header=None)
    df2 = df2[pd.to_numeric(df2[0], errors="coerce").notnull()]
    df2 = {0: pd.to_numeric(df2[0]),
           1: pd.to_numeric(df2[1]),
           2: pd.to_numeric(df2[2])}
    df2 = pd.DataFrame(df2)

    # Create the subplots
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(24, 12))
    plt.suptitle(filenumb, fontsize=40)
    # Plot the aerofoil shape profile
    ax1.set_box_aspect(1)
    ax1.scatter(df1[0], df1[1])
    ax1.set_title("Shape", fontsize=30)
    ax1.set_xlim(-0.1, 1.1)
    ax1.set_ylim(-0.6, 0.6)
    ax1.axhline(y=0, color="black", linestyle="--")
    # Plot the aerofoil CLD profile
    ax2.set_box_aspect(1)
    ax2.scatter(df2[0], df2[1])
    ax2.scatter(df2[0], df2[2])
    ax2.set_title("CLD", fontsize=30)
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(-1.5, 2.5)
    ax2.legend(labels=["Lift", "Drag"], fontsize=20)
    ax2.axhline(y=0, color="black", linestyle="--")
    # Save the figure with the named after the aerofoil number
    plt.savefig(path + f"/{number:02d}")
    plt.show()
