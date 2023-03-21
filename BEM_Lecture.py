# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:40:33 2023

@author: C43353
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions import cld_func, nodal, forces


"""
BEM Using CLD from lectures folder to match with lectures data

Code has -
Constants:
    angular velocity (stay constant to keep power const)

Variables:
    (Using information provided from lectures)
    radial position
    chord length
    pitch angle
    wind speed

Notes -
Not sure if final segmental force and torque plots are correct (have just
removed the lowest radial node to allow plotting)
"""

# Variable Constants
R = 20.5  # Radius (m)
omega = 2.83  # Angular Veolcity (rad/s) (Constant for varying wind speeds)
B = 3  # Number of Blades

# Constant Constants
rho = 1.225  # Air Density (kg/m^3)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

# Variables
# Wind Speeds (m/s)
speeds = np.linspace(5, 20, 20)

# Radial Position of Nodes (m)
segments = [4.5, 5.5, 6.5, 7.5, 8.5,
            9.5, 10.5, 11.5, 12.5, 13.5,
            14.5, 15.5, 16.5, 17.5, 18.5,
            19.5, 20.3]

# Chord Length for Radial Position (m)
chords = [1.63, 1.597, 1.54, 1.481, 1.42,
          1.356, 1.294, 1.229, 1.163, 1.095,
          1.026, 0.955, 0.881, 0.806, 0.705,
          0.545, 0.265]

# Pitch Angle for Radial Position (degree)
thetas = [20, 16.3, 13, 10.05, 7.45,
          5.85, 4.85, 4, 3.15, 2.6,
          2.02, 1.36, 0.77, 0.33, 0.14,
          0.05, 0.02]

fcl, fcd = cld_func("CLD.csv")


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

# Perform calculations over wind speeds
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

    # Perform the calculations over the radial positions
    for m, r in enumerate(segments):
        c = chords[m]  # Chord Length from list
        theta = thetas[m]  # Twist Angle from list

        # Use nodal function to calculate outputs for radial position
        phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, Vrel, Ct, Cpinit = \
            nodal(R, r, V0, c, theta, omega, B, fcl, fcd)

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


# Plot the power output against wind speed
plt.figure(1, figsize=(6, 6))
plt.plot(speeds, np.array(P_out) * 1E-3, marker='o')
plt.title("Power")
plt.xlabel("Wind Speed (m/s)")
plt.xlim(min(speeds), max(speeds))
plt.ylabel("Power Output (kW)")
plt.ylim(0, 600)
plt.axhline(450, color="black", linestyle="--")
plt.show()

# Plot the power coefficient against wind speed
plt.figure(1, figsize=(6, 6))
plt.plot(speeds, Cp_out, marker='o')
plt.title("Power Coefficient")
plt.xlabel("Wind Speed (m/s)")
plt.xlim(min(speeds), max(speeds))
plt.ylabel("Power Coefficient")
plt.ylim(0, 0.5)
plt.show()

# Plot the normalised power coefficient against tip speed ratio
plt.figure(1, figsize=(6, 6))
plt.plot(((omega*R)/np.array(speeds)), (np.array(Cp_out)*(27/16)), marker='o')
plt.title("Normalised")
plt.xlabel(r"$\lambda$ = $\Omega$R/V$_0$ (Tip Speed Ratio)")
plt.xlim(2, 12)
plt.ylabel("C$_p$ $\\times$ 27/16 (Normalised Power Coefficient)")
plt.ylim(0, 1)
plt.show()

# Plot Prandtl Loss Factor on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, F_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Prandtl Loss Factor")
plt.show()

# Plot Angle of Attack on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, alpha_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label=r"Angle of Attack ($\alpha$) / deg")
plt.show()

# Plot Relative Wind Angle on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, phi_out, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label=r"Relative Wind Angle ($\phi$) / rad")
plt.show()

# Plot Normal Nodal Force on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, np.array(fn_out)/1000, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Normal Nodal Force (kN/m)")
plt.show()

# Plot Rotational Nodal Force on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments, speeds, np.array(fr_out)/1000, 50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Rotational Nodal Force (kN/m)")
plt.show()

# Plot Normal Segment Force on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments[-(len(segments)-1):], speeds, np.array(T_out)/1000,
             50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Normal Segment Force (kN)")
plt.show()

# Plot Segmental Torque on contour
plt.figure(1, figsize=(12, 6))
plt.contourf(segments[-(len(segments)-1):], speeds, np.array(tau_out)/1000,
             50, cmap="gist_earth_r")
plt.xlabel("Radial Position / m")
plt.ylabel("Wind Speed / ms$^-$$^1$")
plt.colorbar(label="Segmental Torque (kNm)")
plt.show()
