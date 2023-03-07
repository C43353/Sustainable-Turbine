# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 10:32:22 2023

@author: C43353
"""

import numpy as np
from scipy import interpolate

"""
nodal - function that outputs the data for given radius
nodal_plot - function that gives the progression over the 100 iterations
"""


def cld_func(filename):
    # Open CSV containing aerofoil CLD profile
    file = open(filename)
    CLD = file.read()
    file.close()

    # Split CSV into lines and convert to numeric
    lines = CLD.split('\n')
    lines.pop(0)
    info = []
    for line in lines:
        if len(line) != 0:
            info.append([float(i) for i in line.split(',')])

    # Convert aerofoil info into np.array to allow easier calcuations
    info = np.array(info)

    # Create a function to interpolate values for Cl and Cd
    fcl = interpolate.interp1d(info[:, 0], info[:, 1])
    fcd = interpolate.interp1d(info[:, 0], info[:, 2])

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

        # # Calc New Induction Factor Using Matlab Lecture Code
        # aa = 1 / (K + 1)
        # if aa > ac:
        #     aa = 1 - np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2))
        #     aa = 1 + ((K * (1 - (2 * ac))) / 2) * aa

        # # Calc New Induction Factor Using Converted Matlab Code
        # aa = 1 / (K + 1)
        # if aa > ac:
        #     aa = 1 - ((K * (1 - (2 * ac))) / 2) * (np.sqrt(1 + (4 / K) * (
        #         ((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

        ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    # Relative Wind Speed (Both equations near identical output)
    Vrel = ((1-aa)/(np.sin(phi))) * V0
    # Vrel = ((1 + ar) / (np.cos(phi))) * omega * r

    fn = (1 / 2) * Cn * rho * (Vrel ** 2) * c  # Normal Force at Node
    fr = (1/2) * Cr * rho * (Vrel ** 2) * c  # Rotational Force at Node

    Ct = 4 * aa * (1 - aa)
    Cp = 4 * aa * ((1 - aa) ** 2)

    return [phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, Vrel, Ct, Cp]


def nodal_c(R, r, V0, theta, omega, B, fcl, fcd):
    rho = 1.225  # Air Density (kg/m^3)
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    tsr = (omega * R) / V0  # Tip Speed Ratio
    xi = (omega * r) / V0  # Local Velocity Ratio

    for i in range(100):
        phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

        alpha = phi * (180 / np.pi) - theta  # Angle of Attack

        Cl = float(fcl(alpha))  # Lift Coefficient
        Cd = float(fcd(alpha))  # Drag Coefficient

        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  # Normal Coefficient
        Cr = Cl * np.sin(phi) - Cd * np.cos(phi)  # Tangent Coefficient

        c = (5.6 * (R ** 2)) / (B * Cl * r * (tsr ** 2))  # Chord Length
        s = (c * B) / (2 * np.pi * r)  # Solidity

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

        # # Calc New Induction Factor Using Matlab Lecture Code
        # aa = 1 / (K + 1)
        # if aa > ac:
        #     aa = 1 - np.sqrt(1 + (4 / K) * (((1 - ac) / (1 - (2 * ac))) ** 2))
        #     aa = 1 + ((K * (1 - (2 * ac))) / 2) * aa

        # # Calc New Induction Factor Using Converted Matlab Code
        # aa = 1 / (K + 1)
        # if aa > ac:
        #     aa = 1 - ((K * (1 - (2 * ac))) / 2) * (np.sqrt(1 + (4 / K) * (
        #         ((1 - ac) / (1 - (2 * ac))) ** 2)) - 1)

        ar = 1 / ((4 * F * np.sin(phi) * np.cos(phi)) / (s * Cr) - 1)

    # Relative Wind Speed (Both equations near identical output)
    Vrel = ((1-aa)/(np.sin(phi))) * V0
    # Vrel = ((1 + ar) / (np.cos(phi))) * omega * r

    fn = (1 / 2) * Cn * rho * (Vrel ** 2) * c  # Normal Force at Node
    fr = (1/2) * Cr * rho * (Vrel ** 2) * c  # Rotational Force at Node

    Ct = 4 * aa * (1 - aa)
    Cp = 4 * aa * ((1 - aa) ** 2)

    return [phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, Vrel, Ct, Cp, c]


def nodal_plot(R, r, V0, c, theta, omega, B, fcl, fcd):
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    xi = (omega * r) / V0  # Local Velocity Ratio
    s = (c * B) / (2 * np.pi * r)  # Solidity

    phi_list = []
    alpha_list = []
    Cl_list = []
    Cd_list = []
    Cn_list = []
    Cr_list = []
    F_list = []
    aa_list = []
    ar_list = []

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

        phi_list.append(phi)
        alpha_list.append(alpha)
        Cl_list.append(Cl)
        Cd_list.append(Cd)
        Cn_list.append(Cn)
        Cr_list.append(Cr)
        F_list.append(F)
        aa_list.append(aa)
        ar_list.append(ar)

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

    return [phi_list, alpha_list, Cl_list, Cd_list, Cn_list, Cr_list,
            F_list, aa_list, ar_list]


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
