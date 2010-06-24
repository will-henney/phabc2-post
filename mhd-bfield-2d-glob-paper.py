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

# Box size
try: 
    # look for fourth arg in pc
    worldwidth =  float(sys.argv[4])
except IndexError, ValueError:
    # otherwise, assume 2 pc
    worldwidth = 2.0		# must be float!!!


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
    if f[0].data.shape[0] % 2 == 1:
	# odd number of planes
	data[var] = f[0].data[zcut,:,:]
    else:
	data[var] = f[0].data[zcut-1:zcut+1,:,:].sum(axis=0)

# Subtract off mean velocity of neutral gas
vmean = (data['vx']*data['dd']*data['xn']).sum() / (data['dd']*data['xn']).sum()
data['vx'] -= vmean
print "Mean x-velocity of neutral gas: %i km/s" % (vmean/1.e5)

vmax = N.sqrt(data['vx']**2 + data['vy']**2).max()/1.e5 # km/s
bmax = N.sqrt(data['bx']**2 + data['by']**2).max()*N.sqrt(4.*N.pi)/1.e-6 # micro G
bmin = N.sqrt(data['bx']**2 + data['by']**2).min()*N.sqrt(4.*N.pi)/1.e-6 # micro G

if runid.endswith('2D1024x'):
    time = float(itime)*100.0
else:
    time = float(itime)*1000.0
    
plottitle = r'''
\shortstack{Model \texttt{%s}, $t = %.3f$ Myr,\\
$V_\mathrm{max} = %.1f$~km/s, $V_\mathrm{glob} = %.1f$~km/s, 
$B = [%.0f\ldots %.0f]~\mu$G}
''' % (runid, 1.e-6*time, vmax, vmean/1.e5, bmin, bmax)


# Calculate vector potential: A = \int B_y dx - \int B_x dy
# First stab - simply use cumulative sum for the integration
# fx = N.cumsum(data['bx'][:,0])[:,N.newaxis] # \int B_x dy
# fy = N.cumsum(data['by'], axis=1) # \int B_y dx

# do it the other way round
by0 = data['by'][0,:]
fy = N.cumsum(by0)[N.newaxis,:]
fx = N.cumsum(data['bx'], axis=0)


vecpot = fy - fx

# Do the graph
import pyxgraph, pyx
pyx.text.preamble(r"""\usepackage{mathpazo}""")
# make sure contours are not jagged
pyx.graph.style.line.defaultlineattrs += [pyx.graph.style.linejoin.round]

# set up colors and other stuff
graymap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=0.3, pm3d=[3,3,3])
orangemap = pyxgraph.ColMapper.ColorMapper("white-yellow-red-black", 
					   exponent=1.6, brightness=0.4)
yellowgreen = pyx.color.rgb(0.8,1.0,0.0)
midblue = pyx.color.rgb(0.3,0.2,1.0)
darkgreen = pyx.color.rgb(0.2,0.5,0.0)
midgreen = pyx.color.rgb(0.3,1.0,0.0)
cyan = pyx.color.rgb(0.0,1.0,1.0)
opaque = pyx.color.transparency(0.1)
opaqueish = pyx.color.transparency(0.2)
seethru = pyx.color.transparency(0.5)

# map colors to graph elements
imagecolmap = orangemap
fieldlinecolor = yellowgreen
ifrontcolor = cyan
velocitycolor = cyan
vectorstyles = [velocitycolor, seethru]

ny, nx = vecpot.shape
figwidth = 10.0
worldheight = worldwidth*ny/nx
figheight = ny*figwidth/nx
g = pyxgraph.pyxgraph(xlimits=(0, worldwidth), 
		      ylimits=(0, worldheight),
		      width=figwidth, height=figheight, 
		      key=None,
		      title=plottitle,
		      xlabel='$x$ (pc)',
		      ylabel='$y$ (pc)',
		      )
# image of the pressure
imagedata = vecpot
# imagedata = data['by']
g.pyxplotarray(imagedata[::-1,:], colmap=imagecolmap,
	       xpos=0.0, ypos=0.0, width=worldwidth, height=worldheight, graphcoords=True)

# # image of the pressure
# imagedata = N.log10(data['dd']/(1.3*1.67262158e-24)) # log10 number density
# imagemin = 1.0
# imagemax = 5.0

# g.pyxplotarray(imagedata[::-1,:], colmap=imagecolmap,
# 	       xpos=0.0, ypos=0.0, width=worldwidth, height=worldheight, graphcoords=True)

# cb = pyxgraph.pyxcolorbar(lut=imagecolmap.generate_lut(), frame=g, pos=(1.03,0.0),
# 			  orientation="vertical2",
# 			  width=12*pyx.unit.x_pt,
# 			  minlabel="$n = 10^{%.1f}$" % imagemin, 
# 			  maxlabel="$n = 10^{%.1f}$" % imagemax,
# 			  textattrs=[pyx.trafo.scale(0.7)])
# g.insert(cb)

# Plot contours of vector potential
mycolmap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, brightness=0.2)
xx = N.arange(nx) + 0.5
yy = N.arange(ny) + 0.5
# add 1 to positions to align with image
g.pyxplotcontour(vecpot, (worldwidth/nx)*(xx+1.0), (worldwidth/nx)*(yy+1.0), 
		 levels=200, colors='color', color=fieldlinecolor,
		 lw=0.5, 
		 lineattrs=[opaqueish]
		 )
# Plot contour at ionization front
# g.pyxplotcontour(data['xn'], (worldwidth/nx)*(xx+0.5), (worldwidth/nx)*(yy+0.5), 
# 		 levels=[0.1,0.5,0.9], colors='color', color=ifrontcolor,
# 		 lw=1.0, 
# 		 lineattrs=[opaqueish]
# 		 )



# Changed to b-field for testing 
v = 0.5*N.ones(data['bx'].shape)
theta = N.arctan2(data['by'], data['bx']) # angle of v with x-axis
# # Plot vectors of velocity
# v = N.sqrt(data['vx']**2 + data['vy']**2)/vnorm # normalized magnitude of velocity
# theta = N.arctan2(data['vy'], data['vx']) # angle of v with x-axis
# mx, my = 64, 64*ny/nx			# fixed grid of arrows, independent of resolution 
mx, my = 128, 128*ny/nx			# fixed grid of arrows, independent of resolution 
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
    "k", 0, mx*my-1, "x, y, size, angle = (worldwidth/nx)*x(k), (worldheight/ny)*y(k), s(k), a(k)",
    points=mx*my, context=locals())
g.plot(arrowfunc, [pyx.graph.style.arrow(arrowsize=0.1,
					 arrowattrs=vectorstyles,
					 lineattrs=vectorstyles,
					 )])


# Write file
g.writePDFfile('%s-%s-%i' % (execname.split('.')[0], runid, itime))
