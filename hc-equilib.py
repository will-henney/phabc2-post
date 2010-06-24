"""
Find equilibrium temperatures for PDR gas
"""

import numpy
import scipy
from scipy.optimize import fsolve, bisect
from heatcool import heat, cool

def func(x):
    "Heating = Cooling"
    return heat(dens, Av, D, newvalues=False) - dens*cool(dens, x, True)


densarray = numpy.logspace(0.0, 8.0, num=200)
Tmin, Tmax = 0.1, 1.0e9
for Av, D in [
    [0.0, 0.3],
    [1.0, 0.3],
    [2.0, 0.3],
    [5.0, 0.3],
    [10.0, 0.3],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 1.0],
    [5.0, 1.0],
    [10.0, 1.0],
    ]:
    f = open("hc-equilib-Av%.1f-D%.1f.dat" % (Av, D), "w")
    for dens in densarray:
	T = bisect(func, Tmin, Tmax)
	f.write("%.2e "*2 % (dens, T) + "\n")
    f.close()
