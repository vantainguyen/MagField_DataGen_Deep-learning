# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 15:10:12 2021

@author: s4612643
"""
import random
import time
from numpy import pi
import pickle
import MagField_2


mu0 = 4*pi*10**-7

random.seed(123)
#para = [[R, h, r, \phi, z, J_X, J_Y, J_Z]] 
R_l = []
h_l = []
r_l = []
phi_l = []
z_l = []
J_X_l = []
J_Y_l = []
J_Z_l = []
magnetization = [0, 0, 1]
para = []
output_field = []
    
output_field_axial = []
output_field_azimuthal = []
output_field_radial = []

H_axial_axial_l = []
H_diam_axial_l = []

start_time = time.time()
for i in range(int(200000)):
    #R = random.uniform(3, 50)
    a = random.uniform(55, 100)
    b = random.uniform(55, 100)
    #R_l.append(R)
    h = random.uniform(55, 100)
    h_l.append(h)
    r = random.uniform(155, 200)
    r_l.append(r)
    phi = random.uniform(0, 2*pi)
    phi_l.append(phi)
    z = random.uniform(80, 150)
    z_l.append(z)
    random.shuffle(magnetization)
    J_X, J_Y, J_Z = magnetization
    J_X_l.append(J_X)
    J_Y_l.append(J_Y)
    J_Z_l.append(J_Z)
    para.append([a, b, h, r, phi, z, J_X, J_Y, J_Z])
    # Computation of axial component
    H_diam_axial = MagField_2.B_ellip_cylinder_diam_axial(J_X, J_Y, a/1000, b/1000, h/1000, r/1000, z/1000, phi)/mu0
    H_diam_axial_l.append(H_diam_axial)
    H_axial_axial = MagField_2.B_ellip_cylinder_axial_axial(J_Z, a/1000, b/1000, h/1000, r/1000, z/1000, phi)/mu0
    H_axial_axial_l.append(H_axial_axial)
    H_axial = H_diam_axial + H_axial_axial
    # Computation of azimuthal component
    H_diam_azimuthal = MagField_2.B_ellip_cylinder_diam_azimuthal(J_X, J_Y, a/1000, b/1000, h/1000, r/1000, z/1000, phi)/mu0
    H_axial_azimuthal = MagField_2.B_ellip_cylinder_axial_azimuthal(J_Z, a/1000, b/1000, h/1000, r/1000, z/1000, phi)/mu0
    H_azimuthal = H_diam_azimuthal + H_axial_azimuthal
    # Computation of radial component
    H_diam_radial = MagField_2.B_ellip_cylinder_diam_radial(J_X, J_Y, a/1000, b/1000, h/1000, r/1000, z/1000, phi)/mu0
    H_axial_radial = MagField_2.B_ellip_cylinder_axial_radial(J_Z, a/1000, b/1000, h/1000, r/1000, z/1000, phi)/mu0
    H_radial = H_diam_radial + H_axial_radial
    
    output_field_axial.append([H_axial])
    output_field_azimuthal.append([H_azimuthal])
    output_field_radial.append([H_radial])
    output_field.append([H_axial, H_azimuthal, H_radial])
    #, H_axial/mu0, H_radial/mu0 
    
with open('para_200k_outtest', 'wb') as filename:
    pickle.dump(para, filename)
    
with open('output_field_200k_outtest', 'wb') as filename:
    pickle.dump(output_field, filename)
        
used_time = time.time() - start_time