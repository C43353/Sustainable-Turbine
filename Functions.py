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


def nodal(R, r, V0, c, theta, omega, B, rho, fcl, fcd):
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
        aa = 1 / (K + 1)
        if aa > ac:
            aa = 1 - ((K * (1 - (2 * ac))) / 2) * \
                (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

        ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    # Relative Wind Speed (Both equations near identical output)
    Vrel = ((1 - aa) / (np.sin(phi))) * V0
    # Vrel = ((1 + ar) / (np.cos(phi))) * omega * r

    fn = (1/2) * Cn * rho * (Vrel ** 2) * c  # Normal Force at Node
    fr = (1/2) * Cr * rho * (Vrel ** 2) * c  # Rotational Force at Node

    Ct = 4 * aa * (1 - aa) * F
    if aa > ac:
        Ct = 4 * ((ac**2) + ((1 - (2 * ac)) * aa)) * F  # Linearisation
        # Ct = 4 * aa * (1 - (((5 - (3 * aa)) * aa)/4)) * F  # Glauert's
    # Ct = 4 * aa * (1 - aa)
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


def nodal_twist(R, r, V0, alpha, omega, B, rho, fcl, fcd):
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    xi = (omega * r) / V0  # Local Velocity Ratio
    tsr = (omega * R) / V0  # Tip Speed Ratio

    for i in range(20):
        phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

        theta = phi * (180 / np.pi) - alpha  # Twist Angle

        Cl = float(fcl(alpha))  # Lift Coefficient
        Cd = float(fcd(alpha))  # Drag Coefficient

        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Coefficient
        Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangent Coefficient

        # c = (5.6 * R**2) / (B * Cl * r * (tsr**2))  # https://www.ehow.co.uk/how_7697179_calculate-along-wind-turbine-blade.html
        c = ((8 * np.pi * r) / (B * Cl)) * (1 - np.cos(phi))  # https://ieeexplore.ieee.org/abstract/document/7884538
        # c = (8 * np.pi * r * np.sin(phi)) / (3 * B * Cl * tsr)  # https://www.mdpi.com/1996-1073/13/9/2320
        s = (c * B) / (2 * np.pi * r)  # Solidity

        # Prandtl Loss Factor
        F = (2 / np.pi) * np.arccos(np.exp(- (B * (1 - (r / R))) /
                                           (2 * (r / R) * np.sin(phi))))

        K = (4 * F * (np.sin(phi) ** 2)) / (s * Cn)  # Useful Coefficient

        # Calc New Induction Factor Using Calculation Given in Slides
        aa = 1 / (K + 1)
        if aa > ac:
            aa = 1 - ((K * (1 - (2 * ac))) / 2) * \
                (np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

        ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    return theta, c
