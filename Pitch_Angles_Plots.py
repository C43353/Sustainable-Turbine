# -*- coding: utf-8 -*-
"""
Created on Fri May 19 22:22:29 2023

@author: C43353
"""

import zipfile
import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from Functions import nodal, forces, nodal_chord
import os


""" Inputs """
R = 85  # Radius (m)

V0 = 10  # m/s Nominal Wind Speed

B = 3  # Number of Blades

tsr = 7  # Nominal Tip Speed Ratio (Used to define the angular velocity)


""" Constants """
rho = 1.225  # Air Density (kg/m^3)
N = 17  # Number of Elements

omega = (tsr * 10) / R  # Angular Velocity (dependent on tip speed ratio)
rpm = (omega * 60) / (2 * np.pi)


""" Sweep Parameters """
# Wind Speeds (m/s)
speeds = np.linspace(5, 20, 31)

# Radial Position of Nodes (m)
segments = np.linspace(4.5, R-0.5, N)

# Global Pitch Angle
thetaps = np.linspace(22, 0, 23)


""" Rank Airfoils, Find Optimum Angle of Attack """
# Create dictionaries to store data
data = {}
maxcld1 = {}

# Define the zipfile that stores all the cld data CSVs
zf = zipfile.ZipFile("Aerofoil-data.zip")

for number in range(51):
    # Force the number to be 0x for 0-9
    filenumb = f"{number:02d}"
    name2 = "Profile-" + filenumb + "-CLD.csv"
    df2 = pd.read_csv(zf.open(name2), header=None)
    df2 = df2[pd.to_numeric(df2[0], errors="coerce").notnull()]
    df2 = {0: pd.to_numeric(df2[0]),
           1: pd.to_numeric(df2[1]),
           2: pd.to_numeric(df2[2])}
    data[number] = pd.DataFrame(df2)

    data[number][3] = data[number][1] / data[number][2]

    maxcld1[number] = data[number].loc[data[number][3].idxmax()]

maxcld = {key: val for key, val in sorted(maxcld1.items(),
                                          key=lambda ele: ele[1][3])}

# Select every 3rd profile starting at element 3 (profile 42)
profilescld = list(maxcld.items())[2::3]

# Create a list of the profiles for the turbine
profiles = [x[0] for x in profilescld]

aoa = [list(maxcld1.items())[x][1][0] for x in profiles]


"""
Calculate Twist Angles and Chord Length for Blade
Using Optimum Angle of Attack
"""
thetas = []
chords = []
for m, r in enumerate(segments):
    alpha = aoa[m]  # Twist Angle from list

    profile = profiles[m]  # Profile Cross section from list
    # Create a function to interpolate to find Cl and Cd
    fcl = interpolate.interp1d(data[profile][0], data[profile][1])
    fcd = interpolate.interp1d(data[profile][0], data[profile][2])

    theta, chord = nodal_chord(R, r, V0, alpha, omega, B, rho, fcl, fcd, 3)

    thetas.append(theta)
    chords.append(chord)

# Mirror the chord lengths about the third in the list for realistic sizes
original_chords = list(chords)
chords[0] = 0.9 * chords[2]
chords[1] = (chords[2] + chords[0]) / 2
thetas[0] = 20


""" Perform Calculations Over Varying Global Pitch Angles """
# Initialise the lists for pitch angle output lists
phi_final = []
alpha_final = []
Cl_final = []
Cd_final = []
Cn_final = []
Cr_final = []
F_final = []
aa_final = []
ar_final = []
fn_final = []
fr_final = []
Vrel_final = []
Ct_final = []
Cpinit_final = []
T_final = []
tau_final = []
P_final = []
Cp_final = []

# Calculations for all global pitch angles (thetaps)
for i, thetap in enumerate(thetaps):

    # Initialise the lists for the wind speed output lists
    phi_out = []
    alpha_out = []
    Cl_out = []
    Cd_out = []
    Cn_out = []
    Cr_out = []
    F_out = []
    aa_out = []
    ar_out = []
    fn_out = []
    fr_out = []
    Vrel_out = []
    Ct_out = []
    Cpinit_out = []
    T_out = []
    tau_out = []
    P_out = []
    Cp_out = []
    tsr_out = []

    # Add fixed pitch angle to varying pitch angle
    pitch = thetas + thetap

    """ Perform Calculations Over Varying Wind Speeds """
    # Perform calculations over wind speeds (speeds)
    for n, V0 in enumerate(speeds):

        # Initialise the lists for the radial outputs
        phi_list = []
        alpha_list = []
        Cl_list = []
        Cd_list = []
        Cn_list = []
        Cr_list = []
        F_list = []
        aa_list = []
        ar_list = []
        fn_list = []
        fr_list = []
        Vrel_list = []
        Ct_list = []
        Cpinit_list = []
        T = []
        tau = []

        """ Perform Calculations Over Radial Position on Blade"""
        # Perform the calculations over the radial positions (segments)
        for m, r in enumerate(segments):
            profile = profiles[m]  # Profile Cross section from list

            # Create a function to interpolate to find Cl and Cd
            fcl = interpolate.interp1d(data[profile][0], data[profile][1])
            fcd = interpolate.interp1d(data[profile][0], data[profile][2])

            # Use nodal function to calculate outputs for radial position
            phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, Vrel, Ct, Cpinit = \
                nodal(R, r, V0, chords[m], pitch[m], omega, B, rho, fcl, fcd)

            # Append the outputs for radial position to lists
            phi_list.append(phi)
            alpha_list.append(alpha)
            Cl_list.append(Cl)
            Cd_list.append(Cd)
            Cn_list.append(Cn)
            Cr_list.append(Cr)
            F_list.append(F)
            aa_list.append(aa)
            ar_list.append(ar)
            fn_list.append(fn)
            fr_list.append(fr)
            Vrel_list.append(Vrel)
            Ct_list.append(Ct)
            Cpinit_list.append(Cpinit)

        """ Calculate Turbine Outputs """
        # Use forces function to calculate the normal force and torque on blade
        T, tau = forces(segments, fn_list, fr_list)

        # Calculate the power output of the turbine
        P = B * omega * sum(tau)

        # Calculate the power coefficient of the turbine
        Cp = P / ((1/2) * np.pi * rho * (R ** 2) * (V0 ** 3))

        # Append the outputs for wind speed to lists
        phi_out.append(phi_list)
        alpha_out.append(alpha_list)
        Cl_out.append(Cl_list)
        Cd_out.append(Cd_list)
        Cn_out.append(Cn_list)
        Cr_out.append(Cr_list)
        F_out.append(F_list)
        aa_out.append(aa_list)
        ar_out.append(ar_list)
        fn_out.append(fn_list)
        fr_out.append(fr_list)
        Vrel_out.append(Vrel_list)
        Ct_out.append(Ct_list)
        Cpinit_out.append(Cpinit_list)
        T_out.append(T)
        tau_out.append(tau)
        P_out.append(P)
        Cp_out.append(Cp)
        tsr_out.append((omega * R) / V0)  # Tip Speed Ratio)

    phi_final.append(phi_out)
    alpha_final.append(alpha_out)
    Cl_final.append(Cl_out)
    Cd_final.append(Cd_out)
    Cn_final.append(Cn_out)
    Cr_final.append(Cr_out)
    F_final.append(F_out)
    aa_final.append(aa_out)
    ar_final.append(ar_out)
    fn_final.append(fn_out)
    fr_final.append(fr_out)
    Vrel_final.append(Vrel_out)
    Ct_final.append(Ct_out)
    Cpinit_final.append(Cpinit_out)
    T_final.append(T_out)
    tau_final.append(tau_out)
    P_final.append(P_out)
    Cp_final.append(Cp_out)


""" PLotting """
# Change default saved figure format to svg
# (smaller file size than high resolution png but better quality)
plt.rcParams['savefig.format'] = "svg"


"""Check file path exists, if not create it"""
path = os.path.join("Iterations",
                    os.path.basename(__file__).replace('.py', ''))
isExist = os.path.isdir(path)
if not isExist:
    os.makedirs(path)
    print("Path Created")


P_func = {"thetas": [],
          "P_func": [],
          "P_func_rev": [],
          "P0": [],
          "speed": []}

for i, t in enumerate(thetaps):
    P_func["thetas"].append(t)
    P_func["P_func"].append(interpolate.interp1d(speeds, P_final[i]))
    P_func["P_func_rev"].append(interpolate.interp1d(P_final[i], speeds))

for i in range(len(P_func["P_func_rev"])):
    P_func["speed"].append(P_func["P_func_rev"][i](8E6))

for i in range(len(P_func["thetas"])):
    P_func["P0"].append(P_func["P_func_rev"][i](8E6))


plt.figure(1, figsize=(6, 6))
plt.plot(P_func["P_func_rev"][22](np.linspace(P_final[22][0], 8E6, 12)),
         np.linspace(P_final[22][0], 8E6, 12)*1E-6,
         )  # marker="o")
plt.plot(P_func["speed"],
         np.linspace(8, 8, len(P_func["speed"])),
         c="tab:blue",
         )  # marker="o")
plt.title("Power Against Wind Speed with Pitch Control")
plt.xlabel(r"$V_0$, m/s")
plt.xlim(min(speeds), max(speeds))
plt.ylabel("P, MW")
plt.ylim(0, 10)
# plt.axhline(8, color="black", linestyle="--")
plt.savefig(os.path.join(path, "Pitch Control Power"))
plt.show()


plt.figure(1, figsize=(6, 6))
plt.plot(P_func["P_func_rev"][22](np.linspace(P_final[22][0], 8E6, 12)),
         np.linspace(0, 0, 12),
         )  # marker="o")
plt.plot(P_func["speed"],
         P_func["thetas"],
         c="tab:blue",
         )  # marker="o")
plt.title("Pitch Control Angle Against Wind Speed")
plt.xlabel(r"$V_0$, m/s")
plt.xlim(min(speeds), max(speeds))
plt.ylabel(u"Pitch Control Angle, \N{DEGREE SIGN}")
plt.savefig(os.path.join(path, "Pitch Control Angle"))
plt.show()
