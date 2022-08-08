#! /usr/bin/python3

# Shear stress across a channel cross section
# Based on Eq. 25 & environs from Parker (1978, gravel-bed)
# Then I'm going to see find the magnitude of lateral stresses.

# Started 31 OCT 2021 by ADW in Potsdam @ Aaron Bufe's place!

import numpy as np
from matplotlib import pyplot as plt
#plt.ion()


##########
# Domain #
##########

# Half widths
#b_flat_half = 14 # [m]
b_half = 20 # [m: full half width]
b_bank_half = 6  # [m]
b_flat_half = b_half - b_bank_half

dy = 0.001
y = np.arange(0,b_flat_half+b_bank_half+1E-6,dy)

##################
# Base constants #
##################
rho = 1000.
rho_s = 2650.
g = 9.805

#############
# Variables #
#############
S = 0.001
D = 0.02 # [m]

##########################
# Defined semi-constants #
##########################
tau_star_c = 0.05 # Close to Wong and Parker (2003)
mu = 0.84 # Internal friction for bank gravels
beta = 0.85 # Lift efficiency


###########################
# Derived variables & co. #
###########################

# Roughness: 6.8 D50 (Clifford)
k_s = 6.8 * D

# Bankfull depth s.t. 1.2 tau_star_c maintained
h_c_bf = 1.2 * tau_star_c * (rho_s - rho)/rho * D / S
# For now
h_c = h_c_bf

# Garry's normalized depth thing
R_k = h_c_bf/k_s

psi_approx = (1/12. * np.log(30*R_k) - 5/72.) * (1 + 1 / (2*np.log(30*R_k) - 17/3.))

B_margin_ratio = b_bank_half / (b_flat_half + b_bank_half)

epsilon = (2*h_c/(2*b_bank_half))**2 # p. 135, line 2; 2*b_bank_half = Bs

theta = 1/psi_approx**0.5 * (1-B_margin_ratio) / (B_margin_ratio * epsilon**0.5)

eta = y / b_bank_half - (1 - B_margin_ratio)/B_margin_ratio

#################
# Channel shape #
#################
h = h_c * np.ones(y.shape)
h2 = h_c / (1 - mu*beta) * ( np.cos( eta*np.arccos(mu*beta)) - mu*beta )
h[y >= b_flat_half] = h2[y >= b_flat_half]
z_bed = h_c - h

#########################################################
# Approx curvature on boundary using central difference #
#########################################################
d2h_dy2 = np.zeros(y.shape)
d2h_dy2[1:-1] = ( h[:-2] - 2*h[1:-1] + h[2:] ) / dy**2
# Will just calculate at the needed point in the future. For now:
d2h_dy2_junction = float( d2h_dy2[y == b_flat_half] )

d2f_deta2__yA = b_bank_half / h_c * d2h_dy2_junction
gamma = -0.5 * d2f_deta2__yA


# I lost a factor of negative two somewhere!
# Need to track down, but for now, 2*epsilon
#### (just removed temporariily, in case I was wrong!)
#### (and perhaps I was!)
#### undone
# Eq. 22.3
d__h_c = 1 - 2*epsilon * ( gamma * (1 + 2*psi_approx) / ( np.sinh(theta) + np.cosh(theta) ) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)
z_bed_effective = h_c - d__h_c * h_c

tau_b_x = d__h_c * rho*g*h_c*S
tau_b_y = -np.diff(tau_b_x)/dy * np.sqrt(2)/2. # What this needs to be to keep uniform stress
tau_b_x_mid = (tau_b_x[:-1] + tau_b_x[1:])/2.

tau_b_angle_deg = np.arctan(tau_b_y / -tau_b_x_mid) * 180/np.pi

# Nope, already defined above
#z_bed = -np.hstack((0, np.cumsum( (tau_b_y / -tau_b_x_mid)*dy )))


# For plotting
y_mid = (y[:-1] + y[1:])/2.
#z_bed_effective_mid = (z_bed_effective[:-1] + z_bed_effective[1:])/2.
z_bed_mid = (z_bed[:-1] + z_bed[1:])/2.
tau_b_tot = (tau_b_x_mid**2 + tau_b_y**2)**.5
#plt.plot(y_mid, tau_b_tot)

plt.plot(y_mid, z_bed_mid, 'k-', label='Bed elevation')
plt.legend()
plt.show() # if running from bash


"""
plt.plot(y_mid, tau_b_angle_deg, 'k-', label='Shear angle')
plt.legend()
plt.show() # if running from bash
 """


# Gravity pulling rock down

"""
# NOT NEEDED -- until gravity portion(below) is used
z_bed_ext = np.zeros(len(y)+1)
z_bed_ext[1:-1] = (z_bed[:-1] + z_bed[1:])/2.
z_bed_ext[0] = 0
z_bed_ext[-1] = (np.max(y)+dy - b_flat_half)**2 / 8
# Hard coded parabola here
theta_angle = np.arctan( np.diff(z_bed_ext)/dy )

d__h_c = 1 - epsilon * ( gamma * (1 + 2*psi_approx) / ( np.sinh(theta) + np.cosh(theta) ) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)
#d__h_c = 1 - epsilon * ( gamma * (1 - 2*psi_approx) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)
#d__h_c = 1 - epsilon * ( gamma * (1 - 2*psi_approx) / ( np.sinh(theta_angle) + np.cosh(theta_angle) ) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)


#/ ( np.sinh(theta) * np.cosh(theta)



tau_y_body = (rho_s - rho)*g*D*np.sin(theta_angle)
tau_y_body_mid = (tau_y_body[:-1] + tau_y_body[1:])/2.
# But not accounting for normal force, friction, repose, etc.

plt.plot(y_mid, tau_b_y - tau_y_body_mid)


plt.legend()
plt.show() # if running from bash


tau_shields = tau_b_x / ((rho_s-rho)*g*D)


"""










