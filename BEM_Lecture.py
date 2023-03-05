# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 18:55:18 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

"""
BEM Using CLD from lectures folder to match with lectures data
code has -
Constants:
    angular velocity
    pitch angle
    chord length
Varies:
    radial position
    wind speed

Kind of works, power output doesn't have a peak and Cp is about half what it
should be
"""

# Constants (currently, may change to dependent on r)
B = 3  # Number of Blades
R = 20.5  # Radius (m)
rho = 1.225  # Air Density (kg/m^3)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

# Variable Constants
segments = np.linspace(1, R-0.1, 17)  # Radial Position of Nodes (m)
                                        # (for some reason doesn't work below
                                        # 0.5m and tip causes divide by 0)
speeds = np.linspace(5, 20, 20)  # Wind Speed (m/s)
c = 1.29  # Aerofoil Chord Length (m) (Depends on radial position)
theta = 4.85  # Pitch Angle (degree) (Depends on radial position)
omega = 2.83  # Angular Veolcity (rad/s) (may need to vary with wind speed?)


# Open CSV containing aerofoil CLD profile
file = open('Aerofoil-data\\CLD.csv')
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

data = []

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
T_out = []
tau_out = []
P_out = []
Cp_out = []

for n, V0 in enumerate(speeds):

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
    T = []
    tau = []

    for r in segments:
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

        # Final Output Values
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

    # Calculate Normal Foce and Torque on Each Segment
    for i in range(len(segments)-1):
        j = i + 1
        T.append((1/2) * (fn_list[j] + fn_list[i]) * (
            segments[j] - segments[i]))

        tau.append((1/6) * (((fr_list[j] + fr_list[i]) * ((
            segments[j] ** 2) - (segments[i] ** 2))) + ((fr_list[j] * (
                segments[j] ** 2)) - (fr_list[i] * (segments[i] ** 2))) - ((
                    fr_list[j] - fr_list[i]) * (segments[j] * segments[i]))))

    P = B * omega * sum(tau)
    Cp = P / ((np.pi / 2) * rho * (R ** 2) * (V0 ** 3))

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
    T_out.append(T)
    tau_out.append(tau)
    P_out.append(P)
    Cp_out.append(Cp)

plt.figure(1, figsize=(6, 6))
plt.plot(speeds, np.array(P_out) * 1E-3, marker='o')
plt.title("Power")
plt.xlabel("Wind Speed (m/s)")
plt.ylabel("Power Output (kW)")
plt.show()

plt.figure(1, figsize=(6, 6))
plt.plot(speeds, Cp_out, marker='o')
plt.title("Power Coefficient")
plt.xlabel("Wind Speed (m/s)")
plt.ylabel("Power Coefficient")
plt.ylim(0, 0.5)
plt.show()

plt.figure(1, figsize=(6, 6))
plt.plot(((omega*R)/np.array(speeds)), (np.array(Cp_out)*(27/16)), marker='o')
plt.title("Normalised Coefficients")
plt.xlabel(r"$\lambda$ = $\Omega$R/V$_0$")
plt.ylabel("C$_p$ $\\times$ 27/16")
plt.ylim(0, 1)
plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments, aa_list, marker='o')
# plt.plot(segments, ar_list, marker='o')
# plt.title("Induction Factors")
# plt.xlabel("Radial Position (m)")
# plt.legend(("a", "a'"))
# plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments, F_list, marker='o')
# plt.title("Prandtl Loss Factor")
# plt.xlabel("Radial Position (m)")
# plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments, alpha_list, marker='o')
# plt.plot(segments, phi_list, marker='o')
# plt.title("Angles")
# plt.xlabel("Radial Position (m)")
# plt.legend((r"$\alpha$", r"$\phi$"))
# plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments, np.array(fn_list)/1000, marker='o')
# plt.plot(segments, np.array(fr_list)/1000, marker='o')
# plt.title("Nodal Force")
# plt.xlabel("Radial Position (m)")
# plt.ylabel("f$_{N, i}$, f$_{R, i}$ (kN/m)")
# plt.legend(("F$_{N, i}$", "F$_{R, i}$"))
# plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments[1:], np.array(T)/1000, marker='o')
# plt.plot(segments[1:], np.array(tau)/1000, marker='o')
# plt.title("Segmental Force")
# plt.xlabel("Radial Position (m)")
# plt.ylabel(r"$\tau_i$ (kNm), T$_i$ (kN)")
# plt.legend(("T$_i$", r"$\tau_i$"))
# plt.show()
