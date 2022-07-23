#!/usr/bin/env python
# coding: utf-8

# ## Computation of the magnetic field of a paraboloid with axial magnetization

# In[22]:


# Revised # Import necessary library
from numpy import sqrt, sin, cos, pi, log, tan, pi
from scipy.integrate import quad, dblquad, tplquad

mu0 = 4*pi*10**-7

"""
Parameter definition

The magnetization is axially oriented 

The function to compute the magnetic field is denoted as: 
B_paraboloid_axial_axial (azimuthal, radial) where the meaning of the abbreviations as B - magnetic flux density; 
paraboloid - paraboloic geometry; axial - axial orientation of the magnetization; 
axial (azimuthal, radial) - the axial (azimuthal, radial) components of the magnetic flux density

Note that B = mu0*H

J (T or 10**4 Gauss): the remanence in Tesla of the paraboloid
H (m): the height of the paraboloid
a: the slope coefficient of the paraboloid
(r, z, delta) in (m, m, rad): the coordinates of the computed point in the Cylindrical coordinate system
"""

"""
Computation of the axial component

""" 

# Auxilliary function

def main_func_Ha_1(H, a, r, z, delta, r1, theta):
    

    z1 = a*(r1*1000)**2/1000
        
    main_11 = -(z - z1)*r1
    main_21 = (r1**2 + r**2 - 2*r1*r*cos(delta - theta) + (z - z1)**2)**1.5
    
    return main_11/main_21

def main_func_Ha_2(H, a, r, z, delta, r2, theta):
    
    main_12 = (z - H)*r2
    main_22 = (r2**2 + r**2 - 2*r2*r*cos(delta - theta) + (z - H)**2)**1.5
    
    return main_12/main_22


def B_paraboloid_axial_axial(J, H, a, r, z, delta):
    
    try:
        r_b = sqrt(H*1000/a)/1000
    except ZeroDivisionError:
        raise ZeroDivisionError('The slope coefficient is equal to zero')
        
    integral_B_1 = dblquad(lambda r1, theta: main_func_Ha_1(H, a, r, z, delta, r1, theta), 0, 2*pi, lambda theta: 0,
                        lambda theta: r_b)[0]

    integral_B_2 = dblquad(lambda r2, theta: main_func_Ha_2(H, a, r, z, delta, r2, theta), 0, 2*pi, lambda theta: 0,
                        lambda theta: r_b)[0]
          
    return J*(integral_B_1 + integral_B_2)/4/pi

"""
Computation of the azimuthal component

"""
# Auxilliary function

def main_func_Haz_1(H, a, r, z, delta, r1, theta):
    

    z1 = a*(r1*1000)**2/1000
    
    main_11 = -sin(delta - theta)*r1**2
    main_21 = (r1**2 + r**2 - 2*r1*r*cos(delta - theta) + (z - z1)**2)**1.5
    
    return main_11/main_21

def main_func_Haz_2(H, a, r, z, delta, r2, theta):
    
    main_12 = sin(delta - theta)*r2**2
    main_22 = (r2**2 + r**2 - 2*r2*r*cos(delta - theta) + (z - H)**2)**1.5
    
    return main_12/main_22


def B_paraboloid_axial_azimuthal(J, H, a, r, z, delta):
    
    try:
        r_b = sqrt(H*1000/a)/1000
    except ZeroDivisionError:
        raise ZeroDivisionError('The slope coefficient is equal to zero')
        
    integral_B_1 = dblquad(lambda r1, theta: main_func_Haz_1(H, a, r, z, delta, r1, theta), 0, 2*pi, lambda theta: 0,
                        lambda theta: r_b)[0]

    integral_B_2 = dblquad(lambda r2, theta: main_func_Haz_2(H, a, r, z, delta, r2, theta), 0, 2*pi, lambda theta: 0,
                        lambda theta: r_b)[0]
          
    return J*(integral_B_1 + integral_B_2)/4/pi


"""
Computation of the radial component

"""
# Auxilliary function

def main_func_Hra_1(H, a, r, z, delta, r1, theta):
    

    z1 = a*(r1*1000)**2/1000
    main_11 = -(r - r1*cos(delta - theta))*r1
    main_21 = (r1**2 + r**2 - 2*r1*r*cos(delta - theta) + (z - z1)**2)**1.5
    
    return main_11/main_21

def main_func_Hra_2(H, a, r, z, delta, r2, theta):
    
    main_12 = (r - r2*cos(delta - theta))*r2
    main_22 = (r2**2 + r**2 - 2*r2*r*cos(delta - theta) + (z - H)**2)**1.5
    
    return main_12/main_22


def B_paraboloid_axial_radial(J, H, a, r, z, delta):
    
    try:
        r_b = sqrt(H*1000/a)/1000
    except ZeroDivisionError:
        raise ZeroDivisionError('The slope coefficient is equal to zero')
        
    integral_B_1 = dblquad(lambda r1, theta: main_func_Hra_1(H, a, r, z, delta, r1, theta), 0, 2*pi, lambda theta: 0,
                        lambda theta: r_b)[0]

    integral_B_2 = dblquad(lambda r2, theta: main_func_Hra_2(H, a, r, z, delta, r2, theta), 0, 2*pi, lambda theta: 0,
                        lambda theta: r_b)[0]
          
    return J*(integral_B_1 + integral_B_2)/4/pi


# In[299]:






# In[300]:






# In[ ]:





# In[23]:






# In[27]:


#Plotting with z = -


# In[17]:




# In[330]:






# In[271]:

#fig.show()



# In[209]:

# In[335]:


# In[328]:



# In[46]:


        
        

            


# In[347]:



# In[346]:




# In[67]:



# In[478]:




# In[462]:



# In[294]:



# In[316]:




# In[80]:



# In[93]:




# In[68]:



# ## Computation of the magnetic field of a cone with axial magnetization

# In[38]:


# Import necessary library
from numpy import sqrt, sin, cos, pi,log,tan
from scipy.integrate import quad, dblquad, tplquad
# B_axial = B_KUZ + B_HCZ
mu0 = 4*pi*10**-7
# Computation of H_KUZ

# Auxilliary function 1
"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, “Magnetic field 
distribution of a conical permanent magnet with an application in magnetic resonance imaging,” 
J. Magn. Magn. Mater. 498, 166136 (2019).)

The function to compute the magnetic field is denoted as: 
B_cone_axial_axial where the meaning of the abbreviations as B - magnetic flux density; 
cone - conical geometry; axial - axial orientation of the magnetization; 
axial - the axial component of the magnetic flux density

J: the remanence in Tesla of the lower cone
J_1: the remanance in Tesla of the upper cone
r_b: the base radius of the lower cone
L: the height of the lower cone
2*alpha: the outer apex angle of the lower annular cone 
2*beta: the inner apex angle of the lower annular cone
2*phi: the apex angle of the upper cone
K: the height of the upper cone
(r, z, delta): the coordinates of the computed point in the Cylindrical coordinate system
pxi: the distance along Z axis between the two apexes 
(auxiliary parameter khi denotes the slant height of the cone)
"""
# Computation of the axial component
def main_func_KUZ_axial_axial(theta, r_b, L, r, z, delta):
    pxi = 2*r*cos(delta-theta)
    khi = r**2 + (z-L)**2
    main_1 = 2*(pxi*r_b - 2*khi)/(4*khi - pxi**2)/sqrt(r_b*(r_b - pxi)+khi)
    main_2 = 4*sqrt(khi)/(4*khi - pxi**2)
    return (main_1+main_2)*(z-L)

# B_KUZ

def B_KUZ_axial_axial(J,r_b, L, r, z, delta):
    return J*quad(lambda theta: main_func_KUZ_axial_axial(theta, r_b, L, r, z, delta),-pi,pi)[0]/4/pi

# Auxilliary function 2

def main_func_HCZ_axial_axial(theta, alpha, r_b, L, r, z, delta):
    var_1 = sqrt(L**2+r_b**2)
    var_2 = 2*r*sin(alpha)*cos(delta - theta)+2*z*cos(alpha)
    var_3 = r**2 + z**2
    main_1 = (2*(z*var_2+2*var_3*cos(alpha)-var_2**2*cos(alpha))*var_1-
              4*var_3*z + 2*var_2*var_3*cos(alpha))/(4*var_3-var_2**2)/sqrt(var_1*(var_1-var_2)+var_3)
    main_2 = (-4*var_3*z+2*var_2*var_3*cos(alpha))/(4*var_3-var_2**2)/sqrt(var_3)
    main_3 = cos(alpha)*log((2*(sqrt(var_1*(var_1-var_2)+var_3)+var_1)-var_2)/(2*sqrt(var_3)-var_2))
    return main_1 - main_2 - main_3

# B_HCZ

def B_HCZ_axial_axial(J, alpha, r_b, L, r, z, delta):
    return -J*sin(alpha)**2*quad(lambda theta: main_func_HCZ_axial_axial(theta, alpha, r_b, L, r, z, delta),
                                 -pi,pi)[0]/4/pi

# B_cone_axial
# phi is the half convex angle of the cone
# delta is the coordinate angle of the computed point
def B_cone_axial_axial(J, alpha, L, r, z, delta):
    r_b = L*tan(alpha)
    return B_KUZ_axial_axial(J,r_b, L, r, z, delta) + B_HCZ_axial_axial(J, alpha, r_b, L, r, z, delta)

# Computation of azimuthal component

def B_cone_axial_azimuthal(J, alpha, L, r, z, delta):
    return 0.

# Computation of radial component

def main_func_KUZ_axial_radial(theta, r_b, L, r, z, delta):
    v = cos(delta - theta)
    psi = 2*r*cos(delta - theta)
    chi = r**2 + (z - L)**2
    
    fradi_1 = (2*(psi*r + 2*chi*v - psi**2*v)*r_b - 4*chi*r + 2*psi*chi*v)/((4*chi - psi**2)*sqrt(
    r_b*(r_b - psi) + chi)) - v*log(2*(sqrt(r_b*(r_b - psi) + chi) + r_b) - psi)
    
    fradi_2 = -(-4*chi*r + 2*psi*chi*v)/(4*chi - psi**2)/sqrt(chi) + v*log(2*sqrt(chi) - psi)
    
    return fradi_1 + fradi_2

def B_KUZ_axial_radial(J, r_b, L, r, z, delta):
    return J*quad(lambda theta: main_func_KUZ_axial_radial(theta, r_b, L, r, z, delta),-pi,pi)[0]/4/pi
    
def main_func_HCZ_axial_radial(theta, alpha, r_b, L, r, z, delta):
    
    ksi = sqrt(L**2 + r_b**2)
    ups = 2*r*sin(alpha)*cos(delta - theta) + 2*z*cos(alpha)
    sigma = r**2 + z**2
    tau = sin(alpha)*cos(delta - theta)
    
    fradi_1 = (2*(2*sigma*tau - ups**2*tau + ups*r)*ksi + 2*ups*sigma*tau - 4*sigma*r)/((4*sigma - ups**2)*
                                                                                       sqrt(ksi*(ksi - ups) + sigma))
    fradi_2 = -tau*log((2*(sqrt(ksi*(ksi - ups) + sigma) + ksi) - ups)/(2*sqrt(sigma) - ups)) - (2*ups*sigma*tau -
                                                                                                4*sigma*r)/((4*sigma - ups**2)*sqrt(sigma))
    
    return fradi_1 + fradi_2

def B_HCZ_axial_radial(J, alpha, r_b, L, r, z, delta):
    return -J*sin(alpha)**2*quad(lambda theta: main_func_HCZ_axial_radial(theta, alpha, r_b, L, r, z, delta), 0, 2*pi)[0]/4/pi


def B_cone_axial_radial(J, alpha, L, r, z, delta):
    r_b = L*tan(alpha)
    
    return B_KUZ_axial_radial(J, r_b, L, r, z, delta) + B_HCZ_axial_radial(J, alpha, r_b, L, r, z, delta)


# In[39]:




# In[33]:




# ## Computation of magnetic field of a cone with diametrical magnetization

# In[44]:


# Import necessary library
from numpy import sqrt, sin, cos, pi, log, tan
from scipy.integrate import quad, dblquad, tplquad

mu0 = 4*pi*10**-7

"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, “Analytical computation of the magnetic field
of a conical permanent magnet with arbitrarily
uniform magnetization,” AIP Advances 10, 045208 (2020).)

The function to compute the magnetic field is denoted as: 
B_cone_diam_axial where the meaning of the abbreviations as B - magnetic flux density; 
cone - conical geometry; diam - diametrical orientation of the magnetization; 
axial - the axial component of the magnetic flux density

J: the remanence in Tesla of the cone
h: the height of the cone
2*phi: the apex angle of the cone 
(r, z, beta): the coordinates of the computed point in the Cylindrical coordinate system
"""

# The axial component of the magnetic field

def main_func_diam_axial(h, R, phi, r, z, beta, alpha):
    
    # declaring the auxilliary parameters
    
    khi = sqrt(h**2 + R**2)
    psi = 2*r*sin(phi)*cos(beta - alpha) + 2*z*cos(phi)
    theta = r**2 + z**2
    
    # denoting some components in the main functions
    # first component
    main_1_num = 2*(z*psi + 2*theta*cos(phi) - psi**2*cos(phi))*khi - 4*theta*z + 2*psi*theta*cos(phi)
    main_1_deno = (4*theta - psi**2)*sqrt(khi*(khi - psi) + theta)
    main_1 = main_1_num/main_1_deno
    
    # second component
    main_2_num = -4*theta*z + 2*psi*theta*cos(phi)
    main_2_deno = (4*theta - psi**2)*sqrt(theta)
    main_2 = main_2_num/main_2_deno
    
    # third component
    main_3_num = 2*(sqrt(khi*(khi - psi) + theta) + khi) - psi
    main_3_deno = 2*sqrt(theta) - psi
    main_3 = cos(phi)*log(main_3_num/main_3_deno)
    
    return (main_1 - main_2 - main_3)*cos(alpha)

def B_cone_diam_axial(J, h, phi, r, z, beta):
    R = h*tan(phi)
    integral_B = quad(lambda alpha: main_func_diam_axial(h, R, phi, r, z, beta, alpha), 0, 2*pi)[0]
    return J*sin(2*phi)/8/pi*integral_B

def main_func_diam_azimuthal(h, R, phi, r, z, beta, alpha):
    
    # declaring the auxilliary parameters
    
    ksi = sqrt(h**2 + R**2)
    psi = 2*r*sin(phi)*cos(beta - alpha) + 2*z*cos(phi)
    theta = r**2 + z**2
    
    # first component
    main_1_num = 2*(sqrt(ksi*(ksi - psi) + theta) + ksi) - psi
    main_1_deno = 2*sqrt(theta) - psi
    
    main_1 = log(main_1_num/main_1_deno)
    
    # second component
    main_2_num = 2*(2*theta - psi**2)*ksi + 2*psi*theta
    main_2_deno = (4*theta - psi**2)*sqrt(ksi*(ksi - psi) + theta)
    
    main_2 = main_2_num/main_2_deno
    
    # third component
    main_3_num = 2*psi*theta
    main_3_deno = (4*theta - psi**2)*sqrt(theta)
    
    main_3 = main_3_num/main_3_deno
    
    # main in parentheses
    
    main_4 = main_1 - main_2 + main_3
    
    return main_4*sin(beta - alpha)*cos(alpha)


def B_cone_diam_azimuthal(J, h, phi, r, z, beta):
    R = h*tan(phi)
    integral_B = quad(lambda alpha: main_func_diam_azimuthal(h, R, phi, r, z, beta, alpha), 0, 2*pi)[0]
    return J*sin(phi)*sin(2*phi)/8/pi*integral_B


def main_func_diam_radial(h, R, phi, r, z, beta, alpha):
    
    
    # declaring the auxilliary parameters
    
    ksi = sqrt(h**2 + R**2)
    psi = 2*r*sin(phi)*cos(beta - alpha) + 2*z*cos(phi)
    theta = r**2 + z**2
    lamb = sin(phi)*cos(beta - alpha)
    
    # first component
    
    main_1_num = 2*(2*theta*lamb - psi**2*lamb + psi*r)*ksi + 2*psi*theta*lamb - 4*theta*r
    main_1_deno = (4*theta -psi**2)*sqrt(ksi*(ksi - psi) + theta)
    
    main_1 = main_1_num/main_1_deno
    
    # second component
    
    main_2_num = 2*(sqrt(ksi*(ksi - psi) + theta) + ksi) - psi
    main_2_deno = 2*sqrt(theta) - psi
    
    main_2 = lamb*log(main_2_num/main_2_deno)
    
    # third component
    
    main_3_num = 2*psi*theta*lamb - 4*theta*r
    main_3_deno = (4*theta - psi**2)*sqrt(theta)
    
    main_3 = main_3_num/main_3_deno
    
    # main in parenthesese
    
    main_4 = main_1 - main_2 - main_3
    
    return main_4*cos(alpha)


def B_cone_diam_radial(J, h, phi, r, z, beta):
    R = h*tan(phi)
    integral_B = quad(lambda alpha: main_func_diam_radial(h, R, phi, r, z, beta, alpha), 0, 2*pi)[0]
    return J*sin(2*phi)/8/pi*integral_B
    


# In[45]:



# ## Computation of the magnetic field from cylinder with axial magnetization

# In[3]:


# Importing relevant libraries

import matplotlib as plt
import math as m
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, tplquad

"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, T - F. Lu, W. Robertson, P. Grimshaw, “Magnetic Field Distribution of an Elliptical Permanent Magnet,” 
Progress In Electromagnetics Research C, Vol. 97, 69–82, 2019.)

The function to compute the magnetic field is denoted as: 
B_cylinder_axial_axial where the meaning of the abbreviations as B - magnetic flux density; 
cylinder - cylindrical geometry; axial - axial orientation of the magnetization; 
axial - the axial component of the magnetic flux density

J: the remanence in Tesla of the cylinder
h: the height of the cylinder
R: the radius of the cylinder
(r, z, alpha): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability
K = 1/4/pi/mu0 # auxiliary constant

# Delta funcion
def delta(theta,alpha,r):
    return 2*r*cos(alpha-theta)

# Pxi plus
def pxi_plus(r,z,h):
    return r**2 + (z-h)**2

# Pxi minus
def pxi_minus(r,z,h):
    return r**2 + z**2

# Mag function
def mag(J):
    return J/(4*pi)

# Auxiliary radius
def r0(theta,R):
    a = R
    b = R
     
    return a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)

# Computation of the magnetic field from the upper (plus) surface (dB_plus)
def dB_plus(theta,r, alpha, z, R, h):
    db_func = (2*(delta(theta,alpha,r)*r0(theta,R) - 
                  2*pxi_plus(r,z,h))/(4*pxi_plus(r,z,h) - 
                                      delta(theta,alpha,r)**2)/sqrt((r0(theta,R)*
                                                                     (r0(theta,R)-delta(theta,alpha,r))+
                                                                     pxi_plus(r,z,h))) + 
               4*sqrt(pxi_plus(r,z,h))/\
               (4*pxi_plus(r,z,h)-delta(theta,alpha,r)**2))*(z-h)
    return db_func
    
# Computation of the magnetic field from the lower (minus) surface (dB_minus)
def dB_minus(theta,r, alpha, z, R, h):
    db_func = -(2*(delta(theta,alpha,r)*r0(theta,R) - 
                  2*pxi_minus(r,z,h))/(4*pxi_minus(r,z,h) - 
                                      delta(theta,alpha,r)**2)/sqrt((r0(theta,R)*
                                                                     (r0(theta,R)-delta(theta,alpha,r))+
                                                                     pxi_minus(r,z,h))) + 
                4*sqrt(pxi_minus(r,z,h))/\
                (4*pxi_minus(r,z,h)-delta(theta,alpha,r)**2))*(z)
    return db_func

# Computation of the total magnetic field from both surfaces (dB)
def dB(theta,r, alpha, z, R, h):
    return dB_plus(theta,r, alpha, z, R, h) + dB_minus(theta,r, alpha, z, R, h)

def B_cylinder_axial_axial(J, h, R, r, z, alpha):
    integral_B = quad(lambda theta: dB(theta,r, alpha, z, R, h), -pi, pi)[0]
    return J*integral_B/4/pi


# In[471]:





# ## Computation of magnetic field from elliptical cylinder with axial magnetization

# In[36]:


# Importing relevant libraries

import matplotlib as plt
from numpy import sqrt, sin, cos, pi, log
from scipy.integrate import quad, tplquad

"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, T - F. Lu, W. Robertson, P. Grimshaw, “Magnetic Field Distribution of an Elliptical Permanent Magnet,” 
Progress In Electromagnetics Research C, Vol. 97, 69–82, 2019.)

The function to compute the magnetic field is denoted as: 
B_ellip_cylinder_axial_axial (azimuthal, radial) where the meaning of the abbreviations as B - magnetic flux density; 
ellip_cylinder - elliptical cylinder; axial - axial orientation of the magnetization; 
axial (azimuthal, radial) - the axial (azimuthal, radial) component of the magnetic flux density

J: the remanence in Tesla of the elliptical cylinder
h: the height of the elliptical cylinder
R: the radius of the elliptical cylinder
(r, z, alpha): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability


""" Computation of the axial component B_ellip_cylinder_axial_axial"""

# B_ellip_cylinder_axial_axial = B_plus + B_minus

# Computation of B_plus
def main_plus_ellip_axial_axial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_plus = r**2 + (z - h)**2
    main_1_num = 2*(delta*r0 - 2*xi_plus)
    main_1_deno = (4*xi_plus - delta**2)*sqrt(r0*(r0 - delta) + xi_plus)
    main_1 = main_1_num/main_1_deno
    main_2_num = 4*sqrt(xi_plus)
    main_2_deno = 4*xi_plus - delta**2
    main_2 = main_2_num/main_2_deno
    return (main_1 + main_2)*(z - h)

def B_plus_ellip_axial_axial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_plus_ellip_axial_axial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return J*integral_B/4/pi 

# Computation of B_minus
def main_minus_ellip_axial_axial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_minus = r**2 + z**2
    main_1_num = 2*(delta*r0 - 2*xi_minus)
    main_1_deno = (4*xi_minus - delta**2)*sqrt(r0*(r0 - delta) + xi_minus)
    main_1 = main_1_num/main_1_deno
    main_2_num = 4*sqrt(xi_minus)
    main_2_deno = 4*xi_minus - delta**2
    main_2 = main_2_num/main_2_deno
    return (main_1 + main_2)*(z)

def B_minus_ellip_axial_axial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_minus_ellip_axial_axial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return -J*integral_B/4/pi


def B_ellip_cylinder_axial_axial(J, a, b, h, r, z, alpha):
    return B_plus_ellip_axial_axial(J, a, b, h, r, z, alpha) + B_minus_ellip_axial_axial(J, a, b, h, r, z, alpha)

""" Computation of the azimuthal component B_ellip_cylinder_axial_azimuthal"""

# B_ellip_cylinder_axial_azimuthal = B_azi_plus + B_azi_minus

def main_plus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_plus = r**2 + (z - h)**2
    main_1 = log(2*(sqrt(r0*(r0 - delta) + xi_plus) + r0) - delta)
    main_2_nomi = 2*(2*xi_plus - delta**2)*r0 + 2*delta*xi_plus
    main_2_deno = (4*xi_plus - delta**2)*sqrt(r0*(r0 - delta) + xi_plus)
    main_2 = main_2_nomi/main_2_deno
    main_3 = log(2*sqrt(xi_plus) - delta)
    main_4_nomi = 2*delta*xi_plus
    main_4_deno = (4*xi_plus - delta**2)*sqrt(xi_plus)
    main_4 = main_4_nomi/main_4_deno
    return (main_1 - main_2 - main_3 + main_4)*sin(alpha - theta)

def B_plus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_plus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return J*integral_B/4/pi 


def main_minus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_minus = r**2 + z**2
    main_1 = log(2*(sqrt(r0*(r0 - delta) + xi_minus) + r0) - delta)
    main_2_nomi = 2*(2*xi_minus - delta**2)*r0 + 2*delta*xi_minus
    main_2_deno = (4*xi_minus - delta**2)*sqrt(r0*(r0 - delta) + xi_minus)
    main_2 = main_2_nomi/main_2_deno
    main_3 = log(2*sqrt(xi_minus) - delta)
    main_4_nomi = 2*delta*xi_minus
    main_4_deno = (4*xi_minus - delta**2)*sqrt(xi_minus)
    main_4 = main_4_nomi/main_4_deno

    return (main_1 - main_2 - main_3 + main_4)*sin(alpha - theta)

def B_minus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_minus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return -J*integral_B/4/pi

def B_ellip_cylinder_axial_azimuthal(J, a, b, h, r, z, alpha):
    return B_plus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha) + B_minus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha)

""" Computation of the radial component B_ellip_cylinder_axial_radial"""

# B_ellip_cylinder_axial_radial = B_radi_plus + B_radi_minus

def main_plus_ellip_axial_radial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_plus = r**2 + (z - h)**2
    gamma = cos(alpha - theta)
    main_1_nomi = 2*(delta*r + 2*xi_plus*gamma - delta**2*gamma)*r0 - 4*xi_plus*r + 2*delta*xi_plus*gamma
    main_1_deno = (4*xi_plus - delta**2)*sqrt(r0*(r0 - delta) + xi_plus)
    main_1 = main_1_nomi/main_1_deno
    main_2 = gamma*log(2*(sqrt(r0*(r0 - delta) + xi_plus) + r0) - delta)
    main_3_nomi = -4*xi_plus*r + 2*delta*xi_plus*gamma
    main_3_deno = (4*xi_plus - delta**2)*sqrt(xi_plus)
    main_3 = main_3_nomi/main_3_deno
    main_4 = gamma*log(2*sqrt(xi_plus) - delta)
    return main_1 - main_2 - main_3 + main_4

def B_plus_ellip_axial_radial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_plus_ellip_axial_radial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return J*integral_B/4/pi

def main_minus_ellip_axial_radial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_minus = r**2 + z**2
    gamma = cos(alpha - theta)
    main_1_nomi = 2*(delta*r + 2*xi_minus*gamma - delta**2*gamma)*r0 - 4*xi_minus*r + 2*delta*xi_minus*gamma
    main_1_deno = (4*xi_minus - delta**2)*sqrt(r0*(r0 - delta) + xi_minus)
    main_1 = main_1_nomi/main_1_deno
    main_2 = gamma*log(2*(sqrt(r0*(r0 - delta) + xi_minus) + r0) - delta)
    main_3_nomi = -4*xi_minus*r + 2*delta*xi_minus*gamma
    main_3_deno = (4*xi_minus - delta**2)*sqrt(xi_minus)
    main_3 = main_3_nomi/main_3_deno
    main_4 = gamma*log(2*sqrt(xi_minus) - delta)
    return main_1 - main_2 - main_3 + main_4
    
def B_minus_ellip_axial_radial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_minus_ellip_axial_radial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return -J*integral_B/4/pi

def B_ellip_cylinder_axial_radial(J, a, b, h, r, z, alpha):
    return B_plus_ellip_axial_radial(J, a, b, h, r, z, alpha) + B_minus_ellip_axial_radial(J, a, b, h, r, z, alpha)


# In[ ]:





# In[ ]:





# In[6]:





# In[ ]:





# In[41]:


# code verification


# In[13]:


# magnetic field from elliptical cylinder with diametrical magnetization

# In[11]:


# Importing relevant libraries

import matplotlib as plt
import math as m
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, tplquad

"""
Parameter definition

The magnetization is diametrically oriented (V. T. Nguyen, T - F. Lu, "Modelling of magnetic field distributions of elliptical cylinder permanent
magnets with diametrical magnetization", Journal of magnetism and magnetic materials, vol. 491, 2019)

The function to compute the magnetic field is denoted as: 
B_ellip_cylinder_diam_axial (azimuthal, radial) where the meaning of the abbreviations as B - magnetic flux density; 
ellip_cylinder - elliptical cylinder; diam - diametrical orientation of the magnetization; 
axial (azimuthal, radial) - the axial (azimuthal, radial) component of the magnetic flux density

J: the remanence in Tesla of the cylinder
h: the height of the elliptical cylinder
R: the radius of the elliptical cylinder
(r, z, alpha): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability

"""Computing the axial component of the magnetic field B_ellip_cylinder_diam_axial"""
# B_ellip_cylinder_diam_axial = B_diam_plus + B_diam_minus

# Computation of B_plus
def main_ellip_diam_axial(J_X, J_Y, h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    di = sqrt(r0**2 + r**2 - 2*r0*r*cos(alpha - theta))
    main_1 = 1/sqrt((z - h)**2 + di**2)
    main_2 = 1/sqrt(z**2 + di**2)
    main_3_num = a*b**3*J_X*cos(theta) + b*a**3*J_Y*sin(theta)
    main_3_deno = (b**2*cos(theta)**2 + a**2*sin(theta)**2)**(3/2)
    main_3 = main_3_num/main_3_deno
    
    return (main_1 - main_2)*main_3

def B_ellip_cylinder_diam_axial(J_X, J_Y, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_ellip_diam_axial(J_X, J_Y, h, a, b, r, z, alpha, theta), 0, 2*pi)[0]
    return integral_B/4/pi 

"""Computing the azimuthal component of the magnetic field B_ellip_cylinder_diam_azimuthal"""

def main_ellip_diam_azimuthal(J_X, J_Y, h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    di = sqrt(r0**2 + r**2 - 2*r0*r*cos(alpha - theta))
    main_1 = (h - z)/(di**2*sqrt(di**2 + (h - z)**2))
    main_2 = z/(di**2*sqrt(di**2 + z**2))
    main_3 = r0*sin(alpha - theta)
    main_4_num = a*b**3*J_X*cos(theta) + b*a**3*J_Y*sin(theta)
    main_4_deno = (b**2*cos(theta)**2 + a**2*sin(theta)**2)**(3/2)
    main_4 = main_4_num/main_4_deno
    return (main_1 + main_2)*main_3*main_4

def B_ellip_cylinder_diam_azimuthal(J_X, J_Y, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_ellip_diam_azimuthal(J_X, J_Y, h, a, b, r, z, alpha, theta), 0, 2*pi)[0]
    return integral_B/4/pi


"""Computing the radial component of the magnetic field B_ellip_cylinder_diam_radial"""

def main_ellip_diam_radial(J_X, J_Y, h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    di = sqrt(r0**2 + r**2 - 2*r0*r*cos(alpha - theta))
    main_1 = (h - z)/(di**2*sqrt(di**2 + (h - z)**2))
    main_2 = z/(di**2*sqrt(di**2 + z**2))
    main_3 = (r - r0*cos(alpha - theta))
    main_4_num = a*b**3*J_X*cos(theta) + b*a**3*J_Y*sin(theta)
    main_4_deno = (b**2*cos(theta)**2 + a**2*sin(theta)**2)**(3/2)
    main_4 = main_4_num/main_4_deno
    return (main_1 + main_2)*main_3*main_4

def B_ellip_cylinder_diam_radial(J_X, J_Y, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_ellip_diam_radial(J_X, J_Y, h, a, b, r, z, alpha, theta), 0, 2*pi)[0]
    return integral_B/4/pi


# In[34]:



# In[35]:


# In[12]:


# #### verify for the cylinder

# In[13]:




# ## Computation of magnetic field from a sphere with axial (Z coordinate) magnetization

# In[3]:


# Importing relevant libraries

import matplotlib as plt
import math as m
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, dblquad

"""
Parameter definition

The magnetization is axially oriented 

The function to compute the magnetic field is denoted as: 
H_sphere_axial_axial (azimuthal, radial) where the meaning of the abbreviations as H - magnetic 
field strength; 
sphere - spherical geometry; axial - axial orientation of the magnetization; 
axial, azimuthal and radial - the axial, azimuthal and radial components of the magnetic 
field strength

J: the remanence in Tesla of the sphere
R: the radius of the sphere
(r, z, phi): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability


# Computation of H_sphere_axial_axial
def main_sphere_axial_axial(J, R, r, z, phi, zO, alpha):
    rO = sqrt(R**2 - zO**2)
    main_num = (z - zO)*zO
    main_deno = (rO**2 + r**2 - 2*rO*r*cos(phi - alpha) + (z - zO)**2)**(3/2)
    return main_num/main_deno

def H_sphere_axial_axial(J, R, r, z, phi):
    integral_H = dblquad(lambda alpha, zO: main_sphere_axial_axial(J, R, r, z, phi, zO, alpha), 
                         -R, R, lambda zO: 0, lambda zO: 2*pi)[0]
    return J*integral_H/4/pi/mu0 


# In[4]:




# In[11]:



# In[33]:


# code verification


# ## Computation of magnetic field of diametrically magnetised cylinder using exact equations

# In[59]:


from numpy import sin, cos, pi, tan, arcsin
from cmath import sqrt
#from scipy.special import ellipk as K, ellipe as E, ellipkinc as F, ellipeinc as Einc # ellipk and ellipkinc: Complete and incomplete elliptic 
# integrals of the first kind; ellipe and ellipeinc: Complete and incomplete elliptic integrals of the second kind
from scipy.integrate import quad
import numpy as np



"""
def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x):
        return scipy.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

"""



mu0 = 4*pi*10**-7

def K(m): # Function to compute the complete elliptic integral of the first kind
    func = lambda theta: 1/sqrt((1 - m*sin(theta)**2 + 1e-20))
    
    def real_func(x):
        return np.real(func(x))
    
    def imag_func(x):
        return np.imag(func(x))
    
    real_integral = quad(real_func, 0, pi/2)
    imag_integral = quad(imag_func, 0, pi/2)
    
    return real_integral[0] + 1j*imag_integral[0]


def E(m): # Function to compute the complete elliptic integral of the second kind
    func = lambda theta: sqrt((1 - m*sin(theta)**2 + 1e-20))
    
    def real_func(x):
        return np.real(func(x))
    
    def imag_func(x):
        return np.imag(func(x))
    
    real_integral = quad(real_func, 0, pi/2)
    imag_integral = quad(imag_func, 0, pi/2)
    
    return real_integral[0] + 1j*imag_integral[0]
    
    
    
def F(phi, m): # Function to compute the incomplete elliptic integral of the first kind
    

    func = lambda theta: 1/sqrt(1 - m*sin(theta)**2 + 1e-20)
    
    def real_func(x):
        return np.real(func(x))
    
    def imag_func(x):
        return np.imag(func(x))
    
    real_integral = quad(real_func, 0, phi)
    imag_integral = quad(imag_func, 0, phi)
    
    return real_integral[0] + 1j*imag_integral[0]



def Einc(phi, m): # Function to compute the incomplete elliptic integral of the second kind

    func = lambda theta: sqrt(1 - m*sin(theta)**2 + 1e-20)
    
    def real_func(x):
        return np.real(func(x))
    
    def imag_func(x):
        return np.imag(func(x))
    
    real_integral = quad(real_func, 0, phi)
    imag_integral = quad(imag_func, 0, phi)
    
    return real_integral[0] + 1j*imag_integral[0]

def Piinc(n, phi, m): # Function to compute the incomplete elliptic integral of the third kind
    #main_func_Pi = (1 - n*sin(theta)**2)*(1 - m*sin(theta)**2)**(-1/2)
    
    func = lambda theta: 1/(1 - n*sin(theta)**2 + 1e-20)*1/sqrt(1 - m*sin(theta)**2 + 1e-20)
    
    def real_func(x):
        return np.real(func(x))
    
    def imag_func(x):
        return np.imag(func(x))
    
    real_integral = quad(real_func, 0, phi)
    imag_integral = quad(imag_func, 0, phi)
    
    return real_integral[0] + 1j*imag_integral[0]



def dia_cylinder_Haxial(J, R, h, r, alpha, z):
    
    a = (z - h/2)**2 + R**2 + r**2
    b = (z + h/2)**2 + R**2 + r**2
    c = 2*R*r
    p = 2*c/(c - a)
    u = 2*c/(c - b)
    
    main_1 = (a*K(p) + (c - a)*E(p))/c/sqrt(a - c + 1e-20)
    main_2 = (b*K(u) + (c - b)*E(u))/c/sqrt(b - c + 1e-20)
    
    return J*R*sin(alpha)*(main_1 - main_2)/pi/mu0

def dia_cylinder_Hazimuthal(J, R, h, r, alpha, z):

    a = R**2 + r**2
    b = 2*R*r
    c = h/2 - z
    d = h/2 + z
    t1 = np.float64(0.99999)
    t2 = np.float64(-0.99999)
    zeta = lambda t: sqrt(1 - t**2 + 1e-20)
    eta = lambda t: sqrt(b*(t + 1)/(a + b + c**2 + 1e-20))
    kappa = lambda t: sqrt(b*(t - 1)/(a - b + c**2 + 1e-20))
    lamda = lambda t: sqrt((a - b*t + c**2)/(a + b + c**2 +1e-20))
    nu = lambda t: arcsin(sqrt((t + 1)/2) + 1e-20)
    xi = 2*b/(a + b + c**2)
    sigma = 2*b/(a + b)
    chi = lambda t: arcsin(sqrt((c**2 + a - b*t)/(c**2 + a + b + 1e-20)) + 1e-20)
    psi = (c**2 + a + b)/(c**2 + a - b)
    upsilon = 4*r**2/(c**2 + 4*r**2 + 1e-20)
    
    # auxilliary function
    
    delta = lambda t, a, b, c: -2*c*lamda(t)/(b**2*zeta(t)*sqrt(a - b*t + c**2 + 1e-20))*(-a*zeta(t)*F(nu(t), xi) + (a - b)*zeta(t)*Piinc(sigma, nu(t), xi) +
                                                                          (t + 1)*b*kappa(t)*F(chi(t), psi) + (t + 1)*(-(-a + b - c**2))*kappa(t)*Einc(chi(t), psi))
    
    
    #delta = lambda t, a, b, c: -2*c*1/(b**2*zeta(t)*2)*(-a*zeta(t)*F(nu(t), xi)  +
                                                                          #(t + 1)*b*2*F(chi(t), psi) + (t + 1)*-(-a + b - c**2)*2*Einc(chi(t), psi))
    
    return J*R**2*cos(alpha)/2/pi/mu0*(delta(t2, a, b, c) - delta(t1, a, b, c) +
                                      delta(t2, a, b, d) - delta(t1, a, b, d))
    


# In[67]:



# In[81]:


# Testing functions


# In[82]:
