"""
Test the cooling times
"""
import numpy
from hctcool import cooling_time_years, thermal_time_years

dmin, dmax = 0.0, 8.0
tmin, tmax = 0.0, 4.5

# dgrid, tgrid = numpy.meshgrid(
#     numpy.logspace(dmin, dmax, num=100), # density 1->10^7
#     numpy.logspace(tmin, tmax, num=100), # temperature 10->10^4
#     )

dgrid = numpy.linspace(dmin, dmax, num=200)
tgrid = numpy.linspace(tmin, tmax, num=200)

# meshgrid creates 2D arrays from the 1D arrays
ddgrid, ttgrid = numpy.meshgrid(10**dgrid, 10**tgrid)


# numpy.set_printoptions(precision=2, linewidth=150)
# print "Density"
# print numpy.log10(dgrid)
# print "Temperature"
# print numpy.log10(tgrid)
# print "Cooling time"
# print numpy.log10(tcgrid)

import pyx, pyxgraph, contour

for Av, D in [
    [0.0, 0.3],
    [1.0, 0.3],
    [2.0, 0.3],
    [5.0, 0.3],
    [10.0, 0.3],
    [50.0, 0.3],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 1.0],
    [5.0, 1.0],
    [10.0, 1.0],
    [50.0, 1.0],
    [0.0, 2.0],
    [1.0, 2.0],
    [2.0, 2.0],
    [5.0, 2.0],
    [10.0, 2.0],
    [50.0, 2.0],
    ]:
    tcgrid = thermal_time_years(ddgrid, ttgrid, Av, D)
    g = pyxgraph.pyxgraph(xlimits=(dmin, dmax), ylimits=(tmin, tmax), width=10, height=10)
    cooltimes = [1./x for x in [1.e2, 1.e3, 1.e4, 3.e4, 1.e5, 1.e6]]
    g.pyxplotcontour(1./tcgrid, dgrid, tgrid,
		     levels=cooltimes, colors="color", color=pyx.color.rgb.blue, labels=True)
    g.pyxplotcontour(-1./tcgrid, dgrid, tgrid,
		      levels=cooltimes, colors="color", color=pyx.color.rgb.red, labels=True)
    g.writePDFfile("hc-tcool-contour-Av%.1f-D%.1f" % (Av, D))
