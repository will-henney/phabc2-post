"""
Find lines of constant cooling time
"""
import numpy
from heatcool import heat, cool

boltzman = 1.3806503e-16
year = 3.15576e7
gamma = 5./3.

def cooling_time_years(density, temperature):
    """
    P / (g -1) L
    """
    return boltzman*temperature/((gamma-1.0)*density*cool(density, temperature, True)*year)


def thermal_time_years(density, temperature, Av, D):
    """
    P / (g -1) (L - H)
    """
    return boltzman*temperature/(
	(gamma-1.0)*(
	    density*cool(density, temperature, thousandK=True) - heat(density, Av, D, newvalues=True)
	    )*year)
