# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:01:57 2023

@author: C43353
"""

import numpy as np
from scipy import interpolate
import zipfile
import pandas as pd


def cld_func(filename):
    """
    This function takes one of the number files file names and extracts the
    lift and drag coefficient data from them to create an interpolation
    function.

    Parameters
    ----------
    filename : string
        A string describing a csv file found in Aerofoil-data.zip such as:
            "CLD.csv" or "Profile-00-CLD.csv"

    Returns
    -------
    fcl : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil lift characteristics.
    fcd : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil drag characteristics.
    As the functions are interpolation functions they can only be used for
    inputs between maximum and minimum values in the csv.

    """
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

    return fcl, fcd


def nodal(R, r, V0, c, theta, omega, B, rho, fcl, fcd):
    """
    This function performs BEM theory calculations on an individual node of
    the turbine as described in Dr Yifeng Yangs lecture slides.

    Parameters
    ----------
    R : float
        The outer radius of the wind turbine being analysed in metres.
    r : float
        The nodal position in the wind turbine being analysed in metres.
        Must be 0 < r < R, r cannot be equal to 0 or R as this causes
        calculation errors.
    V0 : float
        The wind speed for the turbine that is being analysed in metres per
        second.
    c : float
        The chord length of the nodal aerofoil in metres.
    theta : float
        The twist angle of the nodal aerofoil in degrees.
    omega : float
        The angular velocity in radians per second.
    B : int
        The number of blades in the turbine.
    rho : float
        The air density for the turbine location in kilograms per metre cubed.
    fcl : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil lift characteristics.
    fcd : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil drag characteristics.

    Returns
    -------
    phi : float
        The final relative wind angle calculated for the aerofoil at the node
        in radians.
    alpha : float
        The final angle of attack for the aerofoil at the node in degrees.
    Cl : float
        The final lift coefficient for the aerofoil at the node.
    Cd : float
        The final drag coefficient for the aerofoil at the node.
    Cn : float
        The final normal force coefficient for the aerofoil at the node.
    Cr : float
        The final tangential force coefficient for the aerofoil at the
        node.
    F : float
        The final Prandtl loss factor for the node.
    aa : float
        The final induction factor for the node.
    ar : float
        The final angular induction factor for the node.
    fn : float
        The normal force at the node in Newtons.
    fr : float
        The rotational force at the node in Newtons.
    Vrel : float
        The relative wind speed at the node in metres per second.
    Ct : float
        Glauert's coeffiient.
    Cp : float
        An initial calculation of the power coefficient nodally (not final)

    """
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    s = (c * B) / (2 * np.pi * r)  # Solidity

    xi = (omega * r) / V0  # Local Velocity Ratio
    # tsr = (omega * R) / V0  # Tip Speed Ratio

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

    return phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, Vrel, Ct, Cp


def nodal_twist(R, r, V0, c, alpha, omega, B, rho, fcl, fcd):
    """
    This function performs BEM calculations to calculate an optimum twist
    angle for the nodes aerofoil.
    The main BEM method is performed as descibed in Dr Yifeng Yangs lecture
    slides.

    Parameters
    ----------
    R : float
        The outer radius of the wind turbine being optimised in metres.
    r : float
        he nodal position in the wind turbine being analysed in metres.
        Must be 0 < r < R, r cannot be equal to 0 or R as this causes
        calculation errors.
    V0 : float
        The nominal wind speed for the turbine that is being optimised in
        metres per second.
    c : float
        The chord length of the nodal aerofoil in metres.
    alpha : float
        The optimum/desired angle of attack for the aerofoil at nominal wind
        speed in degrees.
    omega : float
        The angular velocity in radians per second.
    B : int
        The number of blades in the turbine.
    rho : float
        The air density for the turbine location in kilograms per metre cubed.
    fcl : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil lift characteristics.
    fcd : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil drag characteristics.

    Returns
    -------
    theta : float
        The twist angle required to obtain the desired angle of attack for the
        aerofoil at the node in degrees.

    """
    ac = 1/3  # Critical Induction Factor

    aa = 0.0  # Induction Factor
    ar = 0.0  # Angular Induction Factor

    s = (c * B) / (2 * np.pi * r)  # Solidity

    xi = (omega * r) / V0  # Local Velocity Ratio
    # tsr = (omega * R) / V0  # Tip Speed Ratio

    for i in range(20):
        phi = np.arctan((1 - aa) / ((1 + ar) * xi))  # Relative Wind Angle

        theta = phi * (180 / np.pi) - alpha  # Twist Angle

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

    return theta


def nodal_chord(R, r, V0, alpha, omega, B, rho, fcl, fcd, method=3):
    """
    This function performs BEM calculations to calculate an optimum twist
    angle and chord length for the nodes aerofoil.
    The main BEM method is performed as descibed in Dr Yifeng Yangs lecture
    slides.

    Parameters
    ----------
    R : float
        The outer radius of the wind turbine being optimised in metres.
    r : float
        he nodal position in the wind turbine being analysed in metres.
        Must be 0 < r < R, r cannot be equal to 0 or R as this causes
        calculation errors.
    V0 : float
        The nominal wind speed for the turbine that is being optimised in
        metres per second.
    alpha : float
        The optimum/desired angle of attack for the aerofoil at nominal wind
        speed in degrees.
    omega : float
        The angular velocity in radians per second.
    B : int
        The number of blades in the turbine.
    rho : float
        The air density for the turbine location in kilograms per metre cubed.
    fcl : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil lift characteristics.
    fcd : scipy.interpolate.interpolate.interp1d
        The interpolation funcion describing the aerofoil drag characteristics.
    method : int, optional
        The method used to calculate the chord length, there are three methods.
        1 - Uses the method described here:
            https://www.ehow.co.uk/how_7697179_calculate-along-wind-turbine-blade.html
        2 - Uses the method described here:
            https://www.mdpi.com/1996-1073/13/9/2320
        3 - Uses the method described here:
            https://ieeexplore.ieee.org/abstract/document/7884538
        The default is 3.

    Returns
    -------
    theta : float
        The twist angle required to obtain the desired angle of attack for the
        aerofoil at the node in degrees.
    c : float
        The optimum chord length found for the aerofoli at the node in metres.

    """
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

        if method == 1:
            # https://www.ehow.co.uk/how_7697179_calculate-along-wind-turbine-blade.html
            c = (5.6 * R**2) / (B * Cl * r * (tsr**2))

        elif method == 2:
            # https://www.mdpi.com/1996-1073/13/9/2320
            c = (8 * np.pi * r * np.sin(phi)) / (3 * B * Cl * tsr)

        elif method == 3:
            # https://ieeexplore.ieee.org/abstract/document/7884538
            c = ((8 * np.pi * r) / (B * Cl)) * (1 - np.cos(phi))

        else:
            return("Must select a suitable method of calculating chord length"
                   "1, 2 or 3, default=1")

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


def forces(segments, fn_list, fr_list):
    """
    This function integrates the forces between each node to output lists of
    the segmental normal forces and torque using the method described in
    Dr Yifeng Yangs lecture slides.
    All input lists should be the same length.

    Parameters
    ----------
    segments : list
        This is a numeric list of the nodal positions for the turbine blade.
    fn_list : list
        This is a numeric list of the nodal normal forces for the turbine
        blade.
    fr_list : list
        This is a numeri clist of the nodal rotational forces for the turbine
        blade.

    Returns
    -------
    T : list
        A list of the normal segmental forces between each node.
        List is one element shorter than input lists.
    tau : list
        A list of the segmental torque between each node.
        List is one element shorter than input lists.

    """
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

    return T, tau
