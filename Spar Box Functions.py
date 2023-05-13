# -*- coding: utf-8 -*-
"""
Created on Tue May  9 00:08:25 2023

@author: Marcus Lawrence
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May  9 00:08:25 2023

@author: Marcus Lawrence
"""

#place function in functions file when working
def moment(segments,R,r, fr_list, fn_list):
    M_fr = []
    M_fn = []
    for i in range(len(segments)-1):
        # j = i + 1
        # ri = (segments[i] + segments[j])/2
        M_fr.append( ((fr_list[i] * (R-segments[i]) *(R-segments[i])/2))) 
        M_fn.append( ((fn_list[i] * (R-segments[i]) *(R-segments[i])/2))) 
    return [M_fr, M_fn]
'''
rho - density
A - Area of load bearing part of blade
R - Blade Length
segments - From aerofoil geometries
'''
def M_edge(rho,A,R,segments):
    g = 9.81
    M_edge = []
    for i in range(len(segments)-1):
        M_edge.append(0.5 * rho * A * g * (R-segments[i])**2)
    return M_edge

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
