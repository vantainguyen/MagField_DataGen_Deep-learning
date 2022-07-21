# -*- coding: utf-8 -*-
"""
Created on Tue May  3 15:16:36 2022

@author: s4612643
"""
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 11:19:43 2021

@author: s4612643
"""
import random
import time
import pickle
from numpy import pi
import MagField_3
import matplotlib.pyplot as plt
import numpy as np

mu0 = 4*pi*10**-7

random.seed(123)

para = []
    
output_field_axial = []
output_field_azimuthal = []
output_field_radial = []

start_time = time.process_time()

# Geometrical parameters

phi = pi/6
h = 6
z = 12
#beta = pi/6
J_Y = 1
J_Z = 0
r = 50
ra = []
for beta in np.linspace(0, 2*pi, 100):
   
    
    para.append([phi, h, r, beta, z, J_Y, J_Z])
    #ra.append(r)
    # Computation of axial component
    H_diam_axial = MagField_3.B_cone_diam_axial(J_Y, h/1000, phi, r/1000, z/1000, beta)/mu0
    #H_axial_axial = MagField_3.B_cone_axial_axial(J_Z, h/1000, phi, r/1000, z/1000, beta)/mu0
    H_axial = H_diam_axial #+ H_axial_axial
    
    # Computation of azimuthal component
    H_diam_azimuthal = MagField_3.B_cone_diam_azimuthal(J_Y, h/1000, phi, r/1000, z/1000, beta)/mu0
    #H_axial_azimuthal = MagField_3.B_cone_axial_azimuthal(J_Z, h/1000, phi, r/1000, z/1000, beta)/mu0
    H_azimuthal = H_diam_azimuthal #+ H_axial_azimuthal
    
    # Computation of radial component
    H_diam_radial = MagField_3.B_cone_diam_radial(J_Y, h/1000, phi, r/1000, z/1000, beta)/mu0
    #H_axial_radial = MagField_3.B_cone_axial_radial(J_Z, h/1000, phi, r/1000, z/1000, beta)/mu0
    H_radial = H_diam_radial #+ H_axial_radial
    
    output_field_axial.append(H_axial)
    output_field_azimuthal.append(H_azimuthal)
    output_field_radial.append(H_radial)
    #output_field.append([H_axial, H_azimuthal, H_radial])
    #, H_axial/mu0, H_radial/mu0 
    
used_time = time.process_time() - start_time

with open('para', 'wb') as filename:
    pickle.dump(para, filename)
    
with open('output_field_cone_axial', 'wb') as filename:
    pickle.dump(output_field_axial, filename)

with open('output_field_cone_azimuthal', 'wb') as filename:
    pickle.dump(output_field_azimuthal, filename)
    
with open('output_field_cone_radial', 'wb') as filename:
    pickle.dump(output_field_radial, filename)
    
#plt.plot(ra, output_field_axial)
#plt.show()
#plt.plot(ra, output_field_azimuthal)
#plt.plot(ra, output_field_radial)