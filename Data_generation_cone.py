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

mu0 = 4*pi*10**-7

random.seed(123)

magnetization = [0, 1]
para = []
output_field = []
    
output_field_axial = []
output_field_azimuthal = []
output_field_radial = []

alpha_low = 5*pi/180.0
alpha_up = 80*pi/180.0

start_time = time.process_time()

for i in range(int(2e6)):
    #R = random.uniform(3, 50)
    phi = random.uniform(alpha_low, alpha_up)
    h = random.uniform(3, 50)
    
    r = random.uniform(0, 150)
    z = random.uniform(-75, 75)
    beta = random.uniform(0, 2*pi)
    
    random.shuffle(magnetization)
    J_Y, J_Z = magnetization
    
    para.append([phi, h, r, beta, z, J_Y, J_Z])
    
    # Computation of axial component
    H_diam_axial = MagField_3.B_cone_diam_axial(J_Y, phi, h/1000, r/1000, z/1000, beta)/mu0
    H_axial_axial = MagField_3.B_cone_axial_axial(J_Z, phi, h/1000, r/1000, z/1000, beta)/mu0
    H_axial = H_diam_axial + H_axial_axial
    
    # Computation of azimuthal component
    H_diam_azimuthal = MagField_3.B_cone_diam_azimuthal(J_Y, phi, h/1000, r/1000, z/1000, beta)/mu0
    H_axial_azimuthal = MagField_3.B_cone_axial_azimuthal(J_Z, phi, h/1000, r/1000, z/1000, beta)/mu0
    H_azimuthal = H_diam_azimuthal + H_axial_azimuthal
    
    # Computation of radial component
    H_diam_radial = MagField_3.B_cone_diam_radial(J_Y, phi, h/1000, r/1000, z/1000, beta)/mu0
    H_axial_radial = MagField_3.B_cone_axial_radial(J_Z, phi, h/1000, r/1000, z/1000, beta)/mu0
    H_radial = H_diam_radial + H_axial_radial
    
    output_field_axial.append([H_axial])
    output_field_azimuthal.append([H_azimuthal])
    output_field_radial.append([H_radial])
    output_field.append([H_axial, H_azimuthal, H_radial])
    #, H_axial/mu0, H_radial/mu0 
    
used_time = time.process_time() - start_time

with open('n_para_cone_2m_dataset', 'wb') as filename:
    pickle.dump(para, filename)
    
with open('n_H_cone_2m_dataset', 'wb') as filename:
    pickle.dump(output_field, filename)