"""
Find equilibrium temperatures for PDR gas
"""

import numpy
import scipy
from scipy.optimize import fsolve, bisect
from heatcool import heat, cool

def func(x):
    "Heating = Cooling"
    return heat(dens, Av, D, newvalues=True) - dens*cool(dens, x, True)


densarray = numpy.logspace(0.0, 8.0, num=200)
Tmin, Tmax = 0.1, 1.0e9
for Av, D in [
    [0.0, 0.2],
    [3.0, 0.2],
    [30.0, 0.2],
    [0.0, 0.5],
    [3.0, 0.5],
    [30.0, 0.5],
    [0.0, 2.0],
    [3.0, 2.0],
    [30.0, 2.0],
    ]:
    f = open("hc-equilib-Av%.1f-D%.1f.dat" % (Av, D), "w")
    for dens in densarray:
	T = bisect(func, Tmin, Tmax)
	f.write("%.2e "*2 % (dens, T) + "\n")
    f.close()
