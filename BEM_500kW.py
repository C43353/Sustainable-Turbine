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

# Variables
# Wind Speeds (m/s)
speeds = np.linspace(5, 20, 31)

# Radial Position of Nodes (m)
segments = np.array([4.5, 5.5, 6.5, 7.5, 8.5,
                     9.5, 10.5, 11.5, 12.5, 13.5,
                     14.5, 15.5, 16.5, 17.5, 18.5,
                     19.5, 20.3])

# Chord Length for Radial Position (m)
chords = np.array([1.63, 1.597, 1.54, 1.481, 1.42,
                   1.356, 1.294, 1.229, 1.163, 1.095,
                   1.026, 0.955, 0.881, 0.806, 0.705,
                   0.545, 0.265])

# Global Pitch Angle
thetaps = np.array([20, 16, 12, 8, 5, 0])

fcl, fcd = cld_func("CLD.csv")


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

# Perform calculations over varying global pitch angles
for n, thetap in enumerate(thetaps):

    # Pitch Angle for Radial Position (degree)
    thetas = np.array([20, 16.3, 13, 10.05, 7.45,
                       5.85, 4.85, 4, 3.15, 2.6,
                       2.02, 1.36, 0.77, 0.33, 0.14,
                       0.05, 0.02])

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
    thetas = thetas + thetap

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
                nodal(R, r, V0, c, theta, omega, B, rho, fcl, fcd)

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


""" Plots To Compare to Lectures """

# Plot the nodal forces against radial position (V0 = 9.5 m/s)
plt.figure(1, figsize=(6, 6))
plt.plot(segments, np.array(fn_out[9])/1000, marker='o')
plt.plot(segments, np.array(fr_out[9])/1000, marker='o')
plt.title("Nodal Force (V0 = 9.5 m/s)")
plt.xlabel(r"$r_i$, m")
plt.xlim(4.5, 20.5)
plt.ylabel(r"$f_{N,i}$, $f_{R,i}$, kN/m")
plt.ylim(0, 1.8)
plt.legend(labels=[r"$f_{N,i}$", r"$f_{R,i}$"])
plt.show()

# Plot the normal force and torque against segmental speed ratio (V0 = 9.5 m/s)
plt.figure(1, figsize=(6, 6))
plt.plot((omega * np.array(segments[-(len(segments)-1):])) / speeds[9],
         np.array(T_out[9])/1000, marker='o')
plt.plot((omega * np.array(segments[-(len(segments)-1):])) / speeds[9],
         np.array(tau_out[9])/1000, marker='o')
plt.title("Segmental Forces (V0 = 9.5 m/s)")
plt.xlabel(r"$\xi$$_i$ = $\Omega$$r_i$/$V_0$")
plt.xlim(1, 6.05)
plt.ylabel(r"$\tau$$_i$, kNm; $T_i$, kN")
plt.ylim(0, 3.5)
plt.legend(labels=["$T_i$", r"$\tau$$_i$"])
plt.show()

# Plot Induction Factor in 3D
X, Y = np.meshgrid(segments, speeds)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(Y, X, aa_out, 100)
ax.set_xlabel('Speeds (m/s)')
ax.set_ylabel('Segments (m)')
ax.set_zlabel('Induction Factor')
plt.show()

# Plot Angular Induction Factor in 3D
X, Y = np.meshgrid(segments, speeds)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(Y, X, ar_out, 100)
ax.set_xlabel('Speeds (m/s)')
ax.set_ylabel('Segments (m)')
ax.set_zlabel('Radial Induction Factor')
ax.invert_xaxis()
ylim = ax.get_ylim()
ax.set_yticks(ax.get_yticks())
ax.set_ylim(ylim[::-1])
plt.show()

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

# Plot the power output against wind speed for all global pitch angles
plt.figure(1, figsize=(6, 6))
for i in range(len(P_final)):
    plt.plot(speeds, np.array(list(reversed(P_final)))[i]/1000)
plt.title("Power Against Wind Speed")
plt.xlabel(r"$V_0$, m/s")
plt.xlim(5, 20)
plt.ylabel("P, kW")
plt.ylim(0, 1000)
plt.axhline(450, color="black", linestyle="--")
plt.legend(labels=[r"$\theta$$_p$ = 0",
                   r"$\theta$$_p$ = 5",
                   r"$\theta$$_p$ = 8",
                   r"$\theta$$_p$ = 12",
                   r"$\theta$$_p$ = 16",
                   r"$\theta$$_p$ = 20"])
plt.show()

# Plot the power coefficient against wind speed for all global pitch angles
plt.figure(1, figsize=(6, 6))
for i in range(len(P_final)):
    plt.plot(speeds, np.array(list(reversed(Cp_final)))[i])
plt.title("Power Coefficient Against Wind Speed")
plt.xlabel(r"$V_0$, m/s")
plt.xlim(5, 20)
plt.ylabel("Cp")
plt.ylim(0, 0.5)
plt.legend(labels=[r"$\theta$$_p$ = 0",
                   r"$\theta$$_p$ = 5",
                   r"$\theta$$_p$ = 8",
                   r"$\theta$$_p$ = 12",
                   r"$\theta$$_p$ = 16",
                   r"$\theta$$_p$ = 20"])
plt.show()

# Plot the normal force against power output for all global pitch angles
plt.figure(1, figsize=(6, 6))
for i in range(len(P_final)):
    plt.plot(np.array(list(reversed(P_final)))[i]/1000,
             np.sum(list(reversed(T_final))[i], 1)/1000)
plt.title("Normal Force Against Power Output")
plt.xlabel("P, kW")
plt.xlim(0, 900)
plt.ylabel("T, kN")
# plt.ylim(0, 8)
plt.axvline(450, color="black", linestyle="--")
plt.legend(labels=[r"$\theta$$_p$ = 0",
                   r"$\theta$$_p$ = 5",
                   r"$\theta$$_p$ = 8",
                   r"$\theta$$_p$ = 12",
                   r"$\theta$$_p$ = 16",
                   r"$\theta$$_p$ = 20"])
plt.show()

# """ Non Lecture Plots """

# # Plot Prandtl Loss Factor on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments, speeds, F_out, 50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label="Prandtl Loss Factor")
# plt.show()

# # Plot Angle of Attack on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments, speeds, alpha_out, 50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label=r"Angle of Attack ($\alpha$) / deg")
# plt.show()

# # Plot Relative Wind Angle on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments, speeds, phi_out, 50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label=r"Relative Wind Angle ($\phi$) / rad")
# plt.show()

# # Plot Normal Nodal Force on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments, speeds, np.array(fn_out)/1000,
#              50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label="Normal Nodal Force (kN/m)")
# plt.show()

# # Plot Rotational Nodal Force on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments, speeds, np.array(fr_out)/1000,
#              50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label="Rotational Nodal Force (kN/m)")
# plt.show()

# # Plot Normal Segment Force on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments[-(len(segments)-1):], speeds, np.array(T_out)/1000,
#              50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label="Normal Segment Force (kN)")
# plt.show()

# # Plot Segmental Torque on contour
# plt.figure(1, figsize=(12, 6))
# plt.contourf(segments[-(len(segments)-1):], speeds, np.array(tau_out)/1000,
#              50, cmap="gist_earth_r")
# plt.xlabel("Radial Position / m")
# plt.ylabel("Wind Speed / ms$^-$$^1$")
# plt.colorbar(label="Segmental Torque (kNm)")
# plt.show()
