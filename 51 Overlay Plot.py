# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 16:30:58 2023

@author: C43353
"""

import matplotlib.pyplot as plt
import pandas as pd
import zipfile
import numpy as np

zf = zipfile.ZipFile("Aerofoil-data.zip")

"""Plotting the profiles overlayed"""
# Initialise the figure for the overlay plot
plt.figure(1, figsize=(12, 12))

colour = iter(plt.cm.rainbow(np.linspace(0, 1, 51)))
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

    c = next(colour)
    # Plot the blade profiles overlayed on eachother
    plt.plot(df[0], df[1], c=c)
    plt.xlim(-0.1, 1.1)
    plt.ylim(-0.6, 0.6)
    plt.title("Overlay Plot of Profiles", fontsize=30)
    plt.axhline(y=0, color="black", linestyle="--")
    plt.legend(range(51))
