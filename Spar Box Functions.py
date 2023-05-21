# -*- coding: utf-8 -*-
"""
Created on Tue May  9 00:08:25 2023

@author: Marcus Lawrence
"""

import matplotlib.pyplot as plt
import numpy as np

#Geometry
t = 0.08
b = 0.6
h = 0.4
A = 2*(t*b) + 2*(t*h)
R = 85

'''Function to Calculate Bending Moments'''
def moments(fr_list, fn_list, segments):
    M_fn = []
    M_fr = []
    rw = []
    for i in range(len(segments)-1):
        j = i + 1
        rw.append( (segments[j] + segments[i])/2 )
    #M_fn
    mr_fn = 0.0 #Reaction Bending Moment
    m_tally_fn = 0.0 #Sum of Forces up to point r_i
    mr_fn = sum(np.multiply(fn_list[0:16],rw)) #Total Bending Moment
    #M_fr
    mr_fr = 0.0 #Reaction Bending Moment
    m_tally_fr = 0.0 #Sum of Forces up to point r_i
    mr_fr = sum(np.multiply(fr_list[0:16],rw)) #Total Bending Moment

    for i in range(len(segments)-1):
        j = i + 1
        ri = (segments[j] + segments[i])/2
        
        m_tally_fn += fn_list[i]*ri
        M_fn.append(mr_fn - m_tally_fn) # Reaction Moment - Moments of each Individual Segment

        m_tally_fr += fr_list[i]*ri
        M_fr.append(mr_fr - m_tally_fr)        
    return M_fn, M_fr

'''Function to Calculate Moment due to gravitationaly loading '''
def M_edge(rho,A,R,segments):
    g = 9.81
    M_edge = []
    rw = []
    for i in range(len(segments)-1):
        j = i + 1
        rw.append( (segments[j] + segments[i])/2 )
    for i in range(len(segments)-1):
        j = i + 1
        M_edge.append(0.5 * rho * A * g * (R-rw[i])**2)
    return M_edge

'''Function to Calculate Spar Cap Thickness'''
def thicknessSpar(F,R,b,d,segments,defl,E):
    t_s = F * R ** 3 / (16 * b * d ** 2 * defl * E)
    return t_s

'''Function to Calculate Web Thickness'''
def thicknessWeb(F,R,b,d,segments,defl,E):
    t_s = F * R ** 3 / (16 * b ** 2 * d * defl * E)
    return t_s

'''
W - load [N]
L - Length of blade [m]
d - depth of Spar box [m]
b - breath of Spar box [m]
E - Young's Modulus [N/m^2]
defl - Maximum allowed deflection [m]

This function returns the required thickness for the given parameters
'''
def thickness(W,L,b,d,E,defl):
    return ( (W * L ** 3) / (16 * b * d ** 2 * E * defl))

'''
W - load [N]
L - Length of blade [m]
d - depth of Spar box [m]
b - breath of Spar box [m]
of - fracture stress of material [N/m^2]

This function returns the required thickness for the given parameters
'''
def thickness_stress(W, L, of, b, d):
    return ((W * L) / (4 * of * b * d))

'''
t - thickness [m]
b - breath [m]
d - depth [m]
This function returns the 2nd moment of area for a given set of parameters
'''
def Inertia(t,b,d):
    return (2*t*b*d**2)

m_fn = np.multiply(fn_list,segments)
m_fr = np.multiply(fr_list,segments)

plt.figure(1, figsize=(8, 8))
plt.axes
plt.plot((np.array(segments[0:16])), np.array(M_fn[0:16])/1000, marker='o', color = 'blue', label = 'Beam Bending Moment')
plt.plot((np.array(segments[0:16])), np.array(m_fn[0:16])/1000, marker='o', color = 'orangered', label = 'Elemental Bending Moments')
plt.legend()
plt.xlim(4.5,80)
plt.ylim(0)
plt.title('Plot showing Flapwise Bending Moment')
plt.xlabel('Radial position r [m]')
plt.ylabel("Moment [kNm]")
plt.show()

plt.figure(1, figsize=(8, 8))
plt.axes
plt.plot((np.array(segments[0:16])), np.array(M_fr[0:16])/1000, marker='o', color = 'blue', label = 'Beam Bending Moment')
plt.plot((np.array(segments[0:16])), np.array(m_fr[0:16])/1000, marker='o', color = 'orangered', label = 'Elemental Bending Moments')
plt.legend()
plt.xlim(4.5,80)
plt.ylim(0)
plt.title('Plot showing Edgewise Bending Moment vs. Radial Position')
plt.xlabel('Radial position r [m]')
plt.ylabel("Moment [kNm]")
plt.show()


'''
Geometry of SparBox
'''
tb = 0.08 #Thickness of aerofoil
b = 0.6 
h = 0.4
g = 9.81
A = 2*(tb*b) + 2*(tb*h) #Area
M_ed = M_edge(1800, A, R, segments) #Moment due to gravity 

ts = 0.07
tw = 0.11
d2 = np.array([4.99,4.155,2.56,1.98,1.32,1,0.725,0.605,0.525,0.505,0.465,0.49,0.485,0.43,0.405,0.31])
b2 = np.array([3.49,3.3,3.11,2.91,2.72,2.53,2.34,2.14,1.95,1.76,1.56,1.37,1.18,0.99,0.79,0.6])

b3 = 2
d3 = 1.2

'''
Stresses
'''
s_fn = np.divide(M_fn,(ts*b2*d2))
s_fr = np.divide(M_fr,(tw*b2*d2))
s_fn2 = np.divide(M_fn,(ts*b3*d3))
s_fr2 = np.divide(M_fr,(tw*b3*d3))


'''Plots for stresses'''
plt.figure(1, figsize=(8, 8))
plt.axes
plt.plot((np.array(segments[0:16])), s_fn/1000000, marker='o', color = 'fuchsia')
#plt.plot((np.array(segments[0:16])), np.array(M_fr[0:16])/(0.07*b2*d2)*10e-6, marker='o', color = 'fuchsia')
plt.legend()
plt.xlim(4.5,80)
plt.ylim(0)
plt.title('Plot showing Flapwise Stress Distribution')
plt.xlabel('Radial position [m]')
plt.ylabel("Stress [MPa]")
plt.show()

plt.figure(1, figsize=(8, 8))
plt.axes
plt.plot((np.array(segments[0:16])), s_fr/1000000, marker='o', color = 'deepskyblue')
#plt.plot((np.array(segments[0:16])), (np.array(M_fr[0:16])/(0.11*b2*d2))*10e-6, marker='o', color = 'deepskyblue')
plt.legend()
plt.xlim(4.5,80)
plt.ylim(0)
plt.title('Plot showing Edgewise Stress Distribution')
plt.xlabel('Radial position r [m]')
plt.ylabel("Stress [MPa]")
plt.show()

#d = np.divide(np.array([3.3,2.8,2.5,2,1.7,1.5,1.4,1.4,1.4,1.2,0.9,0.7,0.6,0.6,0.6,0.6,0.6]),2)
defl = 3.5
d = 0.65 #d half the height of the spar box (distance to nuetral axis)
b = 2 #b
E = 40*1e9

'''Finding Thicknesses'''
tspar = thicknessSpar((sum(fn_list))*2,R,b,d,segments,defl,E)
tweb = thicknessWeb((sum(fr_list)+(1500*A*g*R))*2, R, 1, 1.2, segments, defl, E)

print(tspar)
print(tweb)
