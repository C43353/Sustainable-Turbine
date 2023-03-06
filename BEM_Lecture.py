# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 10:49:42 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt

import BEM_Functions as BEM

# imports specific functions from file
# from BEM_Functions import nodal
# from BEM_Functions import cld_func
# from BEM_Functions import forces

"""
BEM Using CLD from lectures folder to match with lectures data

Code has -
Constants:
    angular velocity (with wind)
    pitch angle (with radial position)
    chord length (with radial position)
Variables:
    radial position
    wind speed

Notes -
Kind of works, power output doesn't have a peak and Cp is about half what it
should be according to lecture notes graphs (may be due to having a constant
                                             chord length?)

Appears to work for 20 m radius but not for much bigger without decreasing the
outer radial node

For any given outer radius the minimum radius has to be over 0.5m

Not sure if final segmental force and torque plots are correct (have just
removed the lowest radial node to allow plotting)
"""

# Constants (currently, may change to dependent on r)
B = 3  # Number of Blades
R = 20.5  # Radius (m)
rho = 1.225  # Air Density (kg/m^3)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

# Variable Constants
segments = np.linspace(1, R-0.1, 20)  # Radial Position of Nodes (m)
                                        # (for some reason doesn't work below
                                        # 0.5m and tip causes divide by 0)
speeds = np.linspace(5, 20, 20)  # Wind Speed (m/s)
c = 1.29  # Aerofoil Chord Length (m) (Depends on radial position)
theta = 4.85  # Pitch Angle (degree) (Depends on radial position)
omega = 2.83  # Angular Veolcity (rad/s) (may need to vary with wind speed?)

fcl, fcd = BEM.cld_func('Aerofoil-data\\CLD.csv')

# Initialise the lists for the variable lists
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

    # Initialise the lists for the variables
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

    # Perform the calculations over the radial positions
    phi_list, alpha_list, Cl_list, Cd_list, Cn_list, Cr_list, F_list, \
        aa_list, ar_list, fn_list, fr_list = \
        zip(*[BEM.nodal(R, r, V0, c, theta, omega, B,
                        fcl, fcd) for r in segments])

    T, tau = BEM.forces(segments, fn_list, fr_list)

    P = B * omega * sum(tau)
    Cp = P / ((1/2) * np.pi * rho * (R ** 2) * (V0 ** 3))

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

# for n in range(len(speeds)):
#     plt.figure(1, figsize=(6, 6))
#     plt.plot(segments, aa_out[n], marker='o')
#     plt.title("Induction Factor")
#     plt.xlabel("Radial Position (m)")
# plt.legend((speeds))
# plt.show()

# for n in range(len(speeds)):
#     plt.figure(1, figsize=(6, 6))
#     plt.plot(segments, ar_out[n], marker='o')
#     plt.title("Radial Induction Factor")
#     plt.xlabel("Radial Position (m)")
# plt.legend((speeds))
# plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, aa_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Induction Factor")
plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, ar_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Radial Induction Factor")
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

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, F_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Prandtl Loss Factor")
plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments, alpha_list, marker='o')
# plt.plot(segments, phi_list, marker='o')
# plt.title("Angles")
# plt.xlabel("Radial Position (m)")
# plt.legend((r"$\alpha$", r"$\phi$"))
# plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, alpha_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label=r"Angle of Attack ($\alpha$)")
plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, phi_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label=r"Pitch Angle ($\phi$)")
plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments, np.array(fn_out[7])/1000, marker='o')
# plt.plot(segments, np.array(fr_out[7])/1000, marker='o')
# plt.title("Nodal Force")
# plt.xlabel("Radial Position (m)")
# plt.ylabel("f$_{N, i}$, f$_{R, i}$ (kN/m)")
# plt.legend(("f$_{N, i}$", "f$_{R, i}$"))
# plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, np.array(fn_out)/1000, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Normal Nodal Force (kN/m)")
plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, np.array(fr_out)/1000, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Rotational Nodal Force (kN/m)")
plt.show()

# plt.figure(1, figsize=(6, 6))
# plt.plot(segments[1:], np.array(T_out[7])/1000, marker='o')
# plt.plot(segments[1:], np.array(tau_out[7])/1000, marker='o')
# plt.title("Segmental Force")
# plt.xlabel("Radial Position (m)")
# plt.ylabel(r"$\tau_i$ (kNm), T$_i$ (kN)")
# plt.legend(("T$_i$", r"$\tau_i$"))
# plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments[-(len(segments)-1):], speeds, np.array(T_out)/1000,
             50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Normal Segment Force (kN)")
plt.show()

plt.figure(1, figsize=(12, 6))
plt.contourf(segments[-(len(segments)-1):], speeds, np.array(tau_out)/1000,
             50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Segmental Torque (kNm)")
plt.show()

# Attempting to do a 3D plot like in lecture notes
# ax = plt.axes(projection="3d")
# ax.invert_xaxis()
# ax.invert_yaxis()
# ax.contour3D(segments, speeds, aa_out)
