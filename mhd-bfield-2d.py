"""
Program to calculate the magnetic lines of force and plot them

Version 0.2: WJH 17 Jan 2008 - extend to central xy cuts of 3d models
Version 0.1: WJH 30 Dec 2007

Method is to first calculate the magnetic vector potential and then
contour it. This is only valid for 2D simulations in slab cartesian
symmetry. Or for symmetry planes of a 3d model, where the field can be
guaranteed to be in the plane. In the latter case, we use the central
xy plane.
"""

import numpy as N
import pyfits
import os, sys

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try: 
    runid = sys.argv[1]
    itime = int(sys.argv[2])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME [VNORM]" % execname
    exit

# Velocity normalization
try: 
    # look for third arg in km/s
    vnorm =  1.e5*float(sys.argv[3])
except IndexError, ValueError:
    # otherwise, assume 30 km/s
    vnorm = 3.e6 

# Check whether model is 2d or 3d
is2d = runid.find("2D") != -1
if is2d:
    zcut = 0
else:
    zcut = None

# Read in data arrays for B field
varlist = ["bx", "by", "vx", "vy", "dd", "xn", "pp"] 
data = {}
for var in varlist: 
    f = pyfits.open('%s-%s%4.4i.fits' % (runid, var, itime))
    if zcut is None:
	zcut = f[0].data.shape[0]/2 # half way along first axis
    data[var] = f[0].data[zcut,:,:]



# Calculate vector potential: A = \int B_y dx - \int B_x dy
# First stab - simply use cumulative sum for the integration
fx = N.cumsum(data['bx'][:,0])[:,N.newaxis] # \int B_x dy
fy = N.cumsum(data['by'], axis=1) # \int B_y dx
vecpot = fy - fx

# Do the graph
import pyxgraph, pyx
pyx.text.preamble(r"""\usepackage{mathptmx}""")
ny, nx = vecpot.shape
figwidth = 10
pixmargin = 4 
figheight = (ny+2*pixmargin)*figwidth/(nx + 2*pixmargin)
g = pyxgraph.pyxgraph(xlimits=(-pixmargin, nx+pixmargin), 
		      ylimits=(-pixmargin, ny+pixmargin),
		      width=figwidth, height=figheight, key=None)
# image of the pressure
graymap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=0.3, pm3d=[3,3,3])
# g.pyxplotarray(data['pp'], colmap=graymap,
g.pyxplotarray(N.log10(data['dd']), colmap=graymap,
	       xpos=0.0, ypos=0.0, width=nx, height=ny, graphcoords=True)


# Plot vectors of B-field
# B = N.sqrt(data['bx']**2 + data['by']**2) # magnitude of B-field
# theta = N.arctan2(data['by'], data['bx']) # angle of B with x-axis

# Subtract off mean velocity of neutral gas
vmean = (data['vx']*data['dd']*data['xn']).sum() / (data['dd']*data['xn']).sum()
data['vx'] -= vmean
print "Mean x-velocity of neutral gas: %i km/s" % (vmean/1.e5)

# Plot vectors of velocity
v = N.sqrt(data['vx']**2 + data['vy']**2)/vnorm # normalized magnitude of velocity
theta = N.arctan2(data['vy'], data['vx']) # angle of v with x-axis
# skip = 8 			# do an arrow every this many pixels
# mx = nx/skip
# my = ny/skip
mx, my = 64, 64*ny/nx			# fixed grid of arrows, independent of resolution 
skip = nx/mx

# we abuse a parametric function below, so we express everything in
# terms of a parameter k
import random
x = lambda k: random.randint(skip/4,3*skip/4) + skip*(int(k)/my)
y = lambda k: random.randint(skip/4,3*skip/4) + skip*(int(k)%my)
# s = lambda k: 3.e4*B[y(k),x(k)]
# s = lambda k: 0.5			# all arrows the same size 
s = lambda k: v[y(k),x(k)]
a = lambda k: theta[y(k),x(k)]*180/N.pi   
arrowfunc = pyx.graph.data.paramfunction(
    "k", 0, mx*my-1, "x, y, size, angle = x(k), y(k), s(k), a(k)",
    points=mx*my, context=locals())
vectorstyles = [
    pyx.color.rgb.green, 
#     pyx.color.transparency(0.6),
    ]
g.plot(arrowfunc, [pyx.graph.style.arrow(arrowsize=0.1,
					 arrowattrs=vectorstyles,
					 lineattrs=vectorstyles,
					 )])


# Plot contours of vector potential
mycolmap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, brightness=0.2)
xx = N.arange(nx) + 0.5
yy = N.arange(ny) + 0.5
# g.pyxplotcontour(vecpot, levels=200, colors='map', colmap=mycolmap)
g.pyxplotcontour(vecpot, xx+1.0, yy+1.0, # add 1 to positions to align with image
		 levels=100, colors='color', color='orange',
		 lw=0.8, 
		 lineattrs=[pyx.color.transparency(0.2)]
		 )
# Plot contour at ionization front
g.pyxplotcontour(data['xn'], xx+0.5, yy+0.5, 
		 levels=[0.1,0.5,0.9], colors='color', color='cyan',
		 lw=0.5, 
		 lineattrs=[pyx.color.transparency(0.4)]
		 )

# g.pyxplot(data=((nx/2,),(ny/2,)), style="points", color="green", ps=5, pt=0)

# Write file
g.writePDFfile('%s-%s-%i' % (execname.split('.')[0], runid, itime))
