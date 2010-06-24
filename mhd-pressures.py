import pyfits
import numpy as N
import scipy as S
import scipy.stats as stats
from PIL import Image
import pyx

printextra = 1			# change to 1 to print range, median, quartiles, and linear fit

# Atomic data needed to find the temperature
boltzmann_k = 1.3806503e-16
mp = 1.67262158e-24
mu = 1.3

class PlotVariable(object):
    """
    A variable that we might want to plot
    """
    drawzero = 0
    title = 'Variable name, units'
    shorttitle = 'Variable'
    n = None
    min = None
    max = None

    badvalue = 0.0

    def __init__(self, data, drawzero=0):
	self.drawzero = drawzero
	self.data = data
	self.mask = self.data != self.badvalue

    def settitle(self, title, shorttitle=None):
	self.title = title
	if shorttitle: self.shorttitle = shorttitle

    def setminmaxn(self, min=None, max=None, n=None, logarithmic=0):
	"n is number of bins to use for the images"
	self.n = n
	self.min = min
	self.max = max
	self.logarithmic = logarithmic
	if self.logarithmic:
	    self.min = N.log10(self.min)
	    self.max = N.log10(self.max)
	    # this will fail if setminmaxn has already been called
	    self.data = N.log10(self.data)
	self.setmask((self.data > self.min) & (self.data <= self.max))

    def setmask(self, mask):
	"Adds mask to the data mask"
	self.mask = self.mask & mask

# end class PlotVariable

def squared(data):
    return data**2

def correct_widths(data):
    wsq = data**2 - w_ins**2 - w_th**2 - w_fs**2
    wsq[wsq<0.0] = 0.0
    return N.sqrt(wsq)

class Graph(pyx.graph.graphxy):

    epsilon = 1.e-4		# misses of axis extremes
    labelsonx = True
    labelsony = True
    figwidth = 10
    figheight = 10

    def __init__(self, xvar, yvar, weights=None, gamma=1.0, statslevel=1):
	"""
	Create an individual pyx graph of one variable vs another

	Arguments:
		   xvar - x variable (instance of PlotVariable)
		   yvar - y variable (instance of PlotVariable)

	Derived from graphxy, so can be inserted in a canvas
	"""

	# Mask out any outlying points
	# This makes all the difference!
	xvar.data = xvar.data[xvar.mask & yvar.mask]
	yvar.data = yvar.data[xvar.mask & yvar.mask]
	if weights is None:
	    self.weights = weights
	else:
	    self.weights = weights[xvar.mask & yvar.mask]

	if statslevel > 0:
	    for var in (xvar, yvar):
		if printextra: print "Range of %s is [%0.4f, %0.4f]" % (var.shorttitle,
									var.data.min(),
									var.data.max())
		print "Mean +/- sigma of %s is %0.4f +/- %0.4f" % (var.shorttitle,
								   var.data.mean(),
								   var.data.std())
		if printextra: print "Median of %s is %0.4f" % (var.shorttitle,
								N.median(var.data))

	    print "Correlation coefficient of %s with %s is %0.4f" % (xvar.shorttitle,
								      yvar.shorttitle,
								      N.corrcoef(xvar.data, yvar.data)[0,1])
	    a, b = S.polyfit(xvar.data, yvar.data, 1)
	    if printextra: print "Linear fit is %s = %0.4f %s + %0.4f" % (yvar.shorttitle, a,
									  xvar.shorttitle, b)

	    print "Mean +/- sigma of (%s - %s) is %0.4f +/- %0.4f" % (yvar.shorttitle,
								      xvar.shorttitle,
								      (yvar.data-xvar.data).mean(),
								      (yvar.data-xvar.data).std())
	    quart25 = stats.scoreatpercentile(yvar.data - xvar.data, 25)
	    quart50 = stats.scoreatpercentile(yvar.data - xvar.data, 50)
	    quart75 = stats.scoreatpercentile(yvar.data - xvar.data, 75)

	    if printextra: print "Median of (%s - %s) is %0.4f or %0.4f" % (xvar.shorttitle,
									    yvar.shorttitle,
									    N.median(yvar.data - xvar.data),
									    quart50)
	    if printextra: print "1st & 3rd quartiles of (%s - %s) are %0.4f and %0.4f (diff  %0.4f)" % (
		xvar.shorttitle, yvar.shorttitle,
		quart25, quart75, quart75 - quart25)
	    print

	# Make a grid containing histogram of how many pixels there are
	# with each combination of variables
	xygrid, xedges, yedges = N.histogram2d(
	    xvar.data.flatten(), yvar.data.flatten(),
	    bins=[xvar.n, yvar.n], range=[[xvar.min, xvar.max], [yvar.min, yvar.max]],
	    normed=True, weights=self.weights)

#	print 'Max value is %s' % xygrid.max()
#	print 'Mean value is %s' % xygrid.mean()
#	print 'X range is %0.2d to %0.2d' % (xedges[0], xedges[-1])
#	print 'Y range is %0.2d to %0.2d' % (yedges[0], yedges[-1])

	# Turn it into an image, using PIL
	xyim = Image.new(mode='L', size=(xvar.n, yvar.n))
	# Subtract from 1 to get a "negative" image
	xyscale = 1.0 - xygrid.flatten()/xygrid.max()
	if gamma != 1.0: xyscale = xyscale**(1./gamma)
	xyim.putdata(xyscale, scale=255.0)
	# note that transpose does not work in-place!
	xyim = xyim.transpose(Image.ROTATE_90)
#	xyim = xyim.transpose(Image.FLIP_TOP_BOTTOM)

	# Now make a graph to return
	if self.labelsonx:
	    xpainter = pyx.graph.axis.painter.regular()
	else:
	    xpainter = pyx.graph.axis.painter.linked()
	if self.labelsony:
	    ypainter = pyx.graph.axis.painter.regular()
	else:
	    ypainter = pyx.graph.axis.painter.linked()

	if xvar.logarithmic:
	    xaxis = pyx.graph.axis.logarithmic(title=xvar.title, painter=xpainter,
					       min=10**xvar.min+self.epsilon,
					       max=10**xvar.max-self.epsilon)
	else:
	    xaxis = pyx.graph.axis.linear(title=xvar.title, painter=xpainter,
					  min=xvar.min+self.epsilon, max=xvar.max-self.epsilon)
	if yvar.logarithmic:
	    yaxis = pyx.graph.axis.logarithmic(title=yvar.title, painter=ypainter,
					       min=10**yvar.min+self.epsilon,
					       max=10**yvar.max-self.epsilon)
	else:
	    yaxis = pyx.graph.axis.linear(title=yvar.title, painter=ypainter,
					  min=yvar.min+self.epsilon, max=yvar.max-self.epsilon)

	pyx.graph.graphxy.__init__(self, width=self.figwidth, height=self.figheight,
				   x=xaxis, y=yaxis)
	self.insert(pyx.bitmap.bitmap(0, 0, xyim, width=self.figwidth, height=self.figheight))
	# add the zero-velocity line where needed
	self.dolayout()
	zerolinestyle = [pyx.color.transparency(0.5),
			 pyx.style.linestyle.solid,
			 pyx.style.linewidth.Thin]
	if xvar.drawzero: self.stroke(self.xgridpath(0), zerolinestyle)
	if yvar.drawzero: self.stroke(self.ygridpath(0), zerolinestyle)

# End class Graph

import os, sys

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try:
    runid = sys.argv[1]
    itime = int(sys.argv[2])
    varstring = sys.argv[3]
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME VARSTRING" % execname
    exit

uniform = False
if runid.find('krum') >= 0:
    # krum models write out every 10,000 yrs
    age_Myr = float(itime)/100.0
    # special treatment for uniform initial conditions
    uniform = True
else:
    # others write out every 1000 yrs
    age_Myr = float(itime)/1000.0

for id in [
    'dd', 'pp', 'bx', 'by','bz', 'vx', 'vy', 'vz', 'xn',
    ]:
    try:
	locals()[id] = pyfits.open('%s-%s%4.4i.fits' % (runid, id, itime))['PRIMARY'].data
    except:
	locals()[id] = N.zeros([256, 256, 256])

# dd = pyfits.open('%s-%s%4.4i.fits' % (runid, 'dd', itime))['PRIMARY'].data
# pp = pyfits.open('%s-%s%4.4i.fits' % (runid, 'pp', itime))['PRIMARY'].data
# bx = pyfits.open('%s-%s%4.4i.fits' % (runid, 'bx', itime))['PRIMARY'].data
# by = pyfits.open('%s-%s%4.4i.fits' % (runid, 'by', itime))['PRIMARY'].data
# bz = pyfits.open('%s-%s%4.4i.fits' % (runid, 'bz', itime))['PRIMARY'].data
# vx = pyfits.open('%s-%s%4.4i.fits' % (runid, 'vx', itime))['PRIMARY'].data
# vy = pyfits.open('%s-%s%4.4i.fits' % (runid, 'vy', itime))['PRIMARY'].data
# vz = pyfits.open('%s-%s%4.4i.fits' % (runid, 'vz', itime))['PRIMARY'].data

try:
    # Look for a pre-existing temperature file
    te = pyfits.open('%s-%s%4.4i.fits' % (runid, 'te', itime))['PRIMARY'].data
except IOError:
    # generate the temperature ourselves
    print "Generating temperature on the fly"
    te = (mp * mu / dd) * pp / (2.-xn) / boltzmann_k


pb = 0.5*(bx**2 + by**2 + bz**2)
pt = 0.5*dd*(vx**2 + vy**2 + vz**2)
bb = N.sqrt(4.0*N.pi*(bx**2 + by**2 + bz**2))
dn = dd / (1.3*1.67262158e-24)
for a, atext in [(pp, 'P_gas'), (pb, 'P_mag'), (pt, 'P_ram'), (bb, 'B'), (dn, 'n_gas')] :
    print '%s min/mean/max: %.2e/%.2e/%.2e' % (atext, a.min(), a.mean(), a.max())


pyx.text.set(mode="latex")
pyx.text.preamble("""\usepackage{mathptmx}""")
Graph.figwidth = 6
Graph.figheight = 6

offsetB = 6.281724544	# log(B) for n = 1 pcc, va = 1 km/s

pmin, pmax = -14.5, -7.5
nmin, nmax = 0.8, 7.2
bmin, bmax = -6.1, -3.8
tmin, tmax = 0.2, 4.7
npix = 50

if varstring.endswith("mass"):
    weights = dd
    wstring = "Mass-weighted"
elif varstring.endswith("therm"):
    weights = dd*te
    wstring = "Thermal energy-weighted"
else:
    weights = None
    wstring = "Volume-weighted"

if varstring.startswith("pgas-pmag"):
    p1var = PlotVariable(N.log10(pp))
    p1var.setminmaxn(min=pmin, max=pmax, n=npix)
    p1var.settitle(r'$\log P_\mathrm{gas}$', 'log P_gas')

    p2var = PlotVariable(N.log10(pb))
    p2var.setminmaxn(min=pmin, max=pmax, n=npix)
    p2var.settitle(r'$\log P_\mathrm{mag}$', 'log P_mag')

    g = Graph(p1var, p2var, weights=weights, gamma=1.0, statslevel=0)

elif varstring.startswith("pgas-pram"):
    p1var = PlotVariable(N.log10(pp))
    p1var.setminmaxn(min=pmin, max=pmax, n=npix)
    p1var.settitle(r'$\log P_\mathrm{gas}$', 'log P_gas')

    p2var = PlotVariable(N.log10(pt))
    p2var.setminmaxn(min=pmin, max=pmax, n=npix)
    p2var.settitle(r'$\log P_\mathrm{ram}$', 'log P_ram')

    g = Graph(p1var, p2var, weights=weights, gamma=1.0, statslevel=0)

elif varstring.startswith("pram-pmag"):
    p1var = PlotVariable(N.log10(pt))
    p1var.setminmaxn(min=pmin, max=pmax, n=npix)
    p1var.settitle(r'$\log P_\mathrm{ram}$', 'log P_ram')

    p2var = PlotVariable(N.log10(pb))
    p2var.setminmaxn(min=pmin, max=pmax, n=npix)
    p2var.settitle(r'$\log P_\mathrm{mag}$', 'log P_mag')

    g = Graph(p1var, p2var, weights=weights, gamma=1.0, statslevel=0)

elif varstring.startswith("n-B"):
    dnvar = PlotVariable(N.log10(dn))
    dnvar.setminmaxn(min=nmin, max=nmax, n=npix)
    dnvar.settitle(r'$\log n_\mathrm{gas}$', 'log n')

    bbvar = PlotVariable(N.log10(bb))
    bbvar.setminmaxn(min=bmin, max=bmax, n=npix)
    bbvar.settitle(r'$\log \vert B \vert$', 'log |B|')

    if uniform:
	# mask out the uniform ambient values
	bbvar.setmask(abs(bb - bb[0,0,0]) > 0.01*bb[0,0,0])
	dnvar.setmask(abs(dn - dn[0,0,0]) > 0.01*dn[0,0,0])

    g = Graph(dnvar, bbvar, weights=weights, gamma=1.0, statslevel=0)


# V_a = B/sqrt(4 pi rho)
# => log10(B) = log10(V_5) + 5 + 0.5 log10(4 pi mu mp) + 0.5 log10(rho)
#             = log10(V_5) - 6.281724544

elif varstring.startswith("n-T"):
    dnvar = PlotVariable(N.log10(dn))
    dnvar.setminmaxn(min=nmin, max=nmax, n=npix)
    dnvar.settitle(r'$\log n_\mathrm{gas}$', 'log n')

    tevar = PlotVariable(N.log10(te))
    tevar.setminmaxn(min=tmin, max=tmax, n=npix)
    tevar.settitle(r'$\log T$', 'log T')

    g = Graph(dnvar, tevar, weights=weights, gamma=0.2, statslevel=0)


else:
    raise ValueError, "Invalid varstring: %s" % varstring

if varstring.startswith("n-B"):
    g.plot(pyx.graph.data.function("y(x) = 0.5*x-offsetB", context=locals()), # V_a = 1 km/s
	   [pyx.graph.style.line([pyx.style.linestyle.dashed])])
    g.plot(pyx.graph.data.function("y(x) = 0.5*x-offsetB+1", context=locals()), # V_a = 10 km/s
	   [pyx.graph.style.line([pyx.style.linestyle.dashed])])
    dydx = 0.5*(dnvar.max-dnvar.min)/(bbvar.max-bbvar.min)
    ang = 180*N.arctan(dydx)/N.pi
    dn1 = dnvar.min + 0.5
    x, y = g.pos(dn1, 0.5*dn1-offsetB)
    g.text(x+12*pyx.unit.x_pt, y, r"$V_\mathrm{A} = 1~\mathrm{km\ s^{-1}}$",
	   [pyx.trafo.scale(0.7), pyx.trafo.rotate(ang)])
    x, y = g.pos(dn1, 0.5*dn1-offsetB+1)
    g.text(x+12*pyx.unit.x_pt, y, r"$V_\mathrm{A} = 10~\mathrm{km\ s^{-1}}$",
	   [pyx.trafo.scale(0.7), pyx.trafo.rotate(ang)])
elif varstring.startswith("p"):
    g.plot(pyx.graph.data.function("y(x) = x"),
	   [pyx.graph.style.line([pyx.style.linestyle.dashed])])
    ang = 45.0			# assume square graph
    pp1 = p1var.min + 1.0
    x, y = g.pos(pp1, pp1)
    g.text(x+12*pyx.unit.x_pt, y, r"Equal pressures",
	   [pyx.trafo.scale(0.7), pyx.trafo.rotate(ang)])


#####
#
# New part 12 Apr 2009 - a multitude of lines on the thermodynamics graphs
# Modified 11 May 2009 for turbulence runs
#
#####
if varstring.startswith("n-T"):
    from pyx.style import linestyle, linewidth
    from pyx.color import transparency, rgb
    from pyx.graph.style import line
    from pyx.graph import data
    line.defaultlineattrs += [linewidth.thick, transparency(0.5)]
    myred = rgb(0.5, 0.0, 0.0)
    mygreen = rgb(0.0, 0.5, 0.0)
    myblue = rgb(0.0, 0.0, 0.5)
#     # Free-fall time < 1.e5 years
#     g.plot(data.function("x(y) = 5.31"), [line([linestyle.dashed])])
#     # Plot contours of Jeans instability
#     for MJeans in (0.0, 1.0):
#	g.plot(data.function("y(x) = 0.491 + 0.6667*MJeans + 0.333*(x-3.0)",
#			     context=locals()),
#	       [line([linestyle.solid, linewidth.Thick])])
    # Plot contours of equilibrium T
    execdir = os.path.dirname(sys.argv[0])
    for col, Av, D in [
	[myblue, 0.0, 0.2],
	[mygreen, 3.0, 0.2],
	[myred, 30.0, 0.2],
	[myblue, 0.0, 0.5],
	[mygreen, 3.0, 0.5],
	[myred, 30.0, 0.5],
	[myblue, 0.0, 2.0],
	[mygreen, 3.0, 2.0],
	[myred, 30.0, 2.0],
	]:
	f = open(os.path.join(execdir, "hc-equilib-Av%.1f-D%.1f.dat" % (Av, D)), "r")
	dens, Teq = N.loadtxt(f, unpack=True)
	g.plot(data.values(x=N.log10(dens), y=N.log10(Teq)), [line([linestyle.solid, linewidth.thin, col])])
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



gtitle = r'\texttt{%s}, $t = %.3f$~Myr' % (runid, age_Myr)
g.text(0.5*g.figwidth, g.figheight+6*pyx.unit.x_pt, gtitle, [pyx.text.halign.center])
g.text(0.5*g.figwidth, 12*pyx.unit.x_pt, wstring, [pyx.text.halign.center])
g.writePDFfile('%s-%s-%4.4i-%s' % (execname.split('.')[0], runid, itime, varstring))
