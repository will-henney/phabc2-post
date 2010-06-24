# -*- coding: utf-8 -*-
"""
Reimplementation in python of functions from C2Ray/src/Phab/pdr.f90

09 Apr 2009
"""

import numpy as N

# Module Variables
heat_cosmicray_background = 1.e-27
cosmicrayfactor = 1.0
sigma_dust_vband = 5.e-22
AV_per_tau = 1.086
parsec = 3.085677582e18

class Xray(object):
    Lstar = 2.29e33
    Lcluster = 3.47e33
    rcluster = 0.3 * 3.085677582e18

class Flux(object):
    def __init__(self, radius, QH0=1.e49, fuv_euv_ratio=0.625, xraypars=Xray()):
	rsq = (radius*parsec)**2
	# Ionizing flux in photons/cm^2
	self.euv = QH0 / (4.0*N.pi*rsq)
	# FUV flux in photons/cm^2
	self.fuv = fuv_euv_ratio*self.euv
	# X-ray flux in erg/cm^2 from ionizing star and from YSOs in cluster
	self.xray = xraypars.Lstar/(4*N.pi*rsq) + xraypars.Lcluster/(4*N.pi*(rsq + xraypars.rcluster**2))

def heat(dens, Av, radius, returncomponents=False, radiation="XO", newvalues=True):
    """
    Heating function 
    """
    # Find unatennuated fluxes
    flux = Flux(radius)
    # Cosmic ray heating may optionally be increased by a constant factor
    heat_cosmicray = heat_cosmicray_background * cosmicrayfactor

    def g(x):
	return (x/(1.+x))**2
    if "O" in radiation:
	# Heating due to non-ionizing UV. Important for A_V = 0-5
	heat_fuv = 1.6e-33*flux.fuv*N.exp(-1.9*Av)
	if newvalues:
	    heat_fuv_limit = 3.e-23*(0.03 + (dens/1.e4)) # Limit heating for low densities	else:
	    heat_fuv_limit = 3.e-23*(dens/1.e4)
	heat_fuv = heat_fuv/(1.0 + heat_fuv/heat_fuv_limit)
	# For high density, gas-grain equilibrium at high A_V gives extra heating
	heat_fir = 4.e-27*(flux.euv/1.e12) * N.exp(-0.05*Av) * (dens/(3.e4+dens))**2
    else:
	heat_fuv = 0.0
	heat_fir = 0.0
    if "X" in radiation:
	# Heating due to hard X-rays. Currently, attenuation is ignored
	heat_xray = 6.e-26*(flux.xray/1.e-3)
	if newvalues:
	    # actually, now we have the attenuation back in
	    heat_xray *= N.exp(-0.1*Av/(1.0 + (N.log10(dens)/6.0)**3) )
    else:
	heat_xray = 0.0
    
    heat = heat_fuv + heat_xray + heat_cosmicray + heat_fir

    # Return the total, followed by a dict of the components
    if returncomponents:
	return heat, {"fuv": heat_fuv, "xray": heat_xray, 
		      "cosmicray": heat_cosmicray, "fir": heat_fir}
    else:
	return heat


def cool(dens, temp, thousandK=False):
    """
    The cooling function
    """
    cool = 3.981e-27 * (dens**(-0.4)) * N.sqrt(temp) * N.exp(-72.81*(dens**0.1)/temp)
    if thousandK:
	cool += 2.e-26*(1.e7*N.exp(-1.184e5/(temp+1000)))
    return cool

def coolKI(T):
    """
    Koyama & Inutsuka 2000, as corrected by Vázquez-Semadeni et al 2007
    """
    return 2.e-26*(1.e7*N.exp(-1.184e5/(T+1000)) + 1.4e-2*N.sqrt(T)*N.exp(-92/T))
