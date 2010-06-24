import pyx
import numpy as N
import os, sys

pyx.text.set(mode="latex")
pyx.text.preamble("""\usepackage{mathptmx}""")
figwidth = 6
figheight = 6

pmin, pmax = -12.5, -7.5
nmin, nmax = 0.9, 7.3
bmin, bmax = -6.1, -3.1
tmin, tmax = 0.7, 4.7
npix = 50

g = pyx.graph.graphxy(width=figwidth, height=figheight, 
		      x=pyx.graph.axis.linear(min=nmin, max=nmax, title=r'$\log n_\mathrm{gas}$'),
		      y=pyx.graph.axis.linear(min=tmin, max=tmax, title=r'$\log T$'),
		      )

from pyx.style import linestyle, linewidth
from pyx.color import transparency
from pyx.graph.style import line
from pyx.graph import data
from pyx.trafo import rotate, scale
line.defaultlineattrs += [linewidth.thick, transparency(0.2)]
# Free-fall time < 1.e5 years
g.plot(data.function("x(y) = 5.31"), [line([linestyle.dashed])])
xx, yy = g.pos(5.6, 3.8)
g.text(xx, yy, r"\(t_\mathrm{ff} = 10^5\,\mathrm{yr}\)", [scale(0.7), rotate(90)])
# Plot contours of Jeans instability
for MJeans in (0.0, 1.0): 
    g.plot(data.function("y(x) = 0.491 + 0.6667*MJeans + 0.333*(x-3.0)", 
			 context=locals()), 
	   [line([linestyle.solid, linewidth.Thick])])
ang = N.arctan2(0.3333*(nmax-nmin), (tmax-tmin))*180.0/N.pi
x, y = g.pos(1.6, 0.8)
g.text(x, y, r"\(M_J = 10~M_\odot\)", [scale(0.7), rotate(ang)])
x, y = g.pos(3.6, 0.8)
g.text(x, y, r"\(M_J = 1~M_\odot\)", [scale(0.7), rotate(ang)])

# Plot contours of equilibrium T
execdir = os.path.dirname(sys.argv[0])
for Av, D in [
    [0.0, 0.3],
    [2.0, 0.3],
    [10.0, 0.3],
    [0.0, 1.0],
    [2.0, 1.0],
    [10.0, 1.0],
    ]:
    f = open(os.path.join(execdir, "hc-equilib-Av%.1f-D%.1f.dat" % (Av, D)), "r")
    dens, Teq = N.loadtxt(f, unpack=True)
    g.plot(data.values(x=N.log10(dens), y=N.log10(Teq)), [line([linestyle.solid, linewidth.thin])])
x, y = g.pos(5.7, 3.7)
g.text(x, y, r"\(T_\mathrm{eq}(A_V = 0, r = 0.3)\)", [scale(0.6), rotate(-45)])
x, y = g.pos(4.5, 3.0)
g.text(x, y, r"\(T_\mathrm{eq}(A_V = 0, r = 1.0)\)", [scale(0.6), rotate(-30)])
x, y = g.pos(3.9, 2.7)
g.text(x, y, r"\(T_\mathrm{eq}(A_V = 2, r = 0.3)\)", [scale(0.6), rotate(-23)])
x, y = g.pos(3.0, 2.25)
g.text(x, y, r"\(T_\mathrm{eq}(A_V = 2, r = 1.0)\)", [scale(0.6), rotate(-15)])
x, y = g.pos(2.0, 1.65)
g.text(x, y, r"\(T_\mathrm{eq}(A_V = 10, r = 0.3)\)", [scale(0.6), rotate(-5)])
x, y = g.pos(1.0, 1.5)
g.text(x, y, r"\(T_\mathrm{eq}(A_V = 10, r = 1.0)\)", [scale(0.6), rotate(-3)])

# Plot contours of cooling time
cooltimes = [1.e2, 1.e3, 1.e4, 1.e5]
from hctcool import cooling_time_years
from contour import contour_rectgrid
dgrid = N.linspace(nmin, nmax, num=200) # finer grid necessary to give smooth contours
tgrid = N.linspace(tmin, tmax, num=200)
tcgrid = cooling_time_years(*N.meshgrid(10**dgrid, 10**tgrid))
contours = contour_rectgrid(dgrid, tgrid, tcgrid, cooltimes)
for i, contourset in enumerate(contours):
    for contour in contourset:
	g.plot(data.values(x=contour[0], y=contour[1]), [line([linestyle.dotted])])
x, y = g.pos(1.3, 2.7)
g.text(x, y, r"\(t_\mathrm{cool} = 10^4~\mathrm{yr}\)", [scale(0.7), rotate(65)])
x, y = g.pos(3.05, 2.7)
g.text(x, y, r"\(t_\mathrm{cool} = 10^3~\mathrm{yr}\)", [scale(0.7), rotate(65)])
x, y = g.pos(3.0, 4.0)
g.text(x, y, r"\(t_\mathrm{cool} = 10^2~\mathrm{yr}\)", [scale(0.7), rotate(-10)])

# 

g.writePDFfile("mhd-press-globx-key")

