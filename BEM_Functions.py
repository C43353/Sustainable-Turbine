# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 10:32:22 2023

@author: C43353
"""

import numpy as np


def nodal(R, r, V0, c, theta, omega, B, fcl, fcd):
    rho = 1.225  # Air Density (kg/m^3)
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    xi = (omega * r) / V0  # Local Velocity Ratio
    s = (c * B) / (2 * np.pi * r)  # Solidity

    for i in range(100):
        phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

        alpha = phi * (180 / np.pi) - theta  # Angle of Attack

        Cl = float(fcl(alpha))  # Lift Coefficient
        Cd = float(fcd(alpha))  # Drag Coefficient

        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Coefficient
        Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangent Coefficient

        # Prandtl Loss Factor
        F = (2 / np.pi) * np.arccos(np.exp(- (
            B * (1 - (r / R))) / (2 * (r / R) * np.sin(phi) * r)))

        K = (4 * F * (np.sin(phi) ** 2)) / (s * Cn)  # Useful Coefficient

        aa = 1 / (K + 1)  # Calc New Induction Factor
        # From Lecture Calculation
        if aa > ac:
            aa = 1 - ((K * (1 - (2 * ac))) / 2) * (
                np.sqrt(1 + (4 / K) * (((1 - ac) / (
                    1 - (2 * ac))) ** 2)) - 1)

        # From Converting Matlab Code
        # if aa > ac:
            # aa = 1 + np.sqrt(
            # 1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2))
            # aa = 1 - (K * ((1 - (2 * ac)) / (2)))

        ar = 1 / ((4 * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    # Relative Wind Speed (Both equations near identical output)
    Vrel = ((1-aa)/(np.sin(phi))) * V0
    # Vrel = ((1 + ar) / (np.cos(phi))) * omega * r

    fn = (1 / 2) * Cn * rho * (Vrel ** 2) * c  # Normal Force at Node
    fr = (1/2) * Cr * rho * (Vrel ** 2) * c  # Rotational Force at Node

    return [phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr]
