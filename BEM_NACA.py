# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 15:00:11 2023

@author: C43353
"""

import numpy as np
from scipy import interpolate

from BEM_Functions import nodal

# Constants (currently)
R = 20  # Radius (m)
segments = np.linspace(9, 19.5, 17)  # Radial Position of Segments (m)
c = 3.256  # Aerofoil Chord Length (m) (Depends on radial position)
B = 3  # Number of Blades
omega = 2.83  # Angular Veolcity (rad/s)
theta = 4.188  # Pitch Angle (degree)
V0 = 10  # Wind Speed (m/s)
ac = 1/3  # Critical Induction Factor (Just use 1/3 as stated in lecture)

# Open CSV containing aerofoil CLD profile
file = open('Aerofoil-data\\NACA 63-415 AIRFOIL Aerodynamic Data.csv')
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

# final = []

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

for r in segments:
    phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar, fn, fr, _, _, _ = nodal(R, r, V0, c,
                                                                   theta,
                                                                   omega, B,
                                                                   fcl, fcd)

    # An array of the iterations of the code (unnecessary?)
    # outputs = np.array([phi_list, alpha_list, Cl_list, Cd_list, Cn_list,
    # Cr_list, F_list, aa_list, ar_list])

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

    # final.append([phi, alpha, Cl, Cd, Cn, Cr, F, aa, ar])
