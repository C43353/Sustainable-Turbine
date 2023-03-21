# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:01:57 2023

@author: C43353
"""

import numpy as np
from scipy import interpolate
from pathlib import Path
import zipfile
import pandas as pd


def cld_func(filename):
    # Open zip file
    zf = zipfile.ZipFile("Aerofoil-data.zip")

    df = pd.read_csv(zf.open(filename), header=None)
    df = df[pd.to_numeric(df[0], errors="coerce").notnull()]
    df = {0: pd.to_numeric(df[0]),
          1: pd.to_numeric(df[1]),
          2: pd.to_numeric(df[2])}
    df = pd.DataFrame(df)

    # Create a function to interpolate values for Cl and Cd
    fcl = interpolate.interp1d(df[0], df[1])
    fcd = interpolate.interp1d(df[0], df[2])

    return [fcl, fcd]


def nodal(R, r, V0, c, theta, omega, B, fcl, fcd):
    rho = 1.225  # Air Density (kg/m^3)
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    s = (c * B) / (2 * np.pi * r)  # Solidity
    # tsr = (omega * R) / V0  # Tip Speed Ratio
    xi = (omega * r) / V0  # Local Velocity Ratio

    for i in range(100):
        phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

        alpha = phi * (180 / np.pi) - theta  # Angle of Attack

        Cl = float(fcl(alpha))  # Lift Coefficient
        Cd = float(fcd(alpha))  # Drag Coefficient

        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Coefficient
        Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangent Coefficient

        # Prandtl Loss Factor
        F = (2 / np.pi) * np.arccos(np.exp(- (B * (1 - (r / R))) /
                                           (2 * (r / R) * np.sin(phi))))

        K = (4 * F * (np.sin(phi) ** 2)) / (s * Cn)  # Useful Coefficient

        # Calc New Induction Factor Using Calculation Given in Slides
        if K > (ac ** -1) - 1:
            aa = 1 / (K + 1)
        if K <= (ac ** -1) - 1:
            aa = 1 - ((K * (1 - (2 * ac))) / 2) * \
                (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

        ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    # Relative Wind Speed (Both equations near identical output)
    Vrel = ((1 - aa) / (np.sin(phi))) * V0
    # Vrel = ((1 + ar) / (np.cos(phi))) * omega * r

    fn = (1/2) * Cn * rho * (Vrel ** 2) * c  # Normal Force at Node
    fr = (1/2) * Cr * rho * (Vrel ** 2) * c  # Rotational Force at Node

    Ct = 4 * aa * (1 - aa)
    Cp = 4 * aa * ((1 - aa) ** 2)

    return [phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, Vrel, Ct, Cp]


def forces(segments, fn_list, fr_list):
    T = []
    tau = []
    # Calculate Normal Force and Torque on Each Segment
    for i in range(len(segments)-1):
        j = i + 1
        T.append((1/2) * (fn_list[j] + fn_list[i]) * (
            segments[j] - segments[i]))

        tau.append((1/6) * (((fr_list[j] + fr_list[i]) * ((
            segments[j] ** 2) - (segments[i] ** 2))) + ((fr_list[j] * (
                segments[j] ** 2)) - (fr_list[i] * (segments[i] ** 2))) - ((
                    fr_list[j] - fr_list[i]) * (segments[j] * segments[i]))))

    return [T, tau]
