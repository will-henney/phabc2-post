"""
Program to calculate the magnetic lines of force and plot them

Version 0.3: WJH 05 Apr 2009 - specialized for globule paper
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
# Box size
try: 
    # look for fourth arg in pc
    worldwidth =  float(sys.argv[4])
except IndexError, ValueError:
    # otherwise, assume 20 pc
    worldwidth = 20.0		# must be float!!!

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


vmax = N.sqrt(data['vx']**2 + data['vy']**2).max()/1.e5 # km/s
bmax = N.sqrt(data['bx']**2 + data['by']**2).max()*N.sqrt(4.*N.pi)/1.e-6 # micro G
bmin = N.sqrt(data['bx']**2 + data['by']**2).min()*N.sqrt(4.*N.pi)/1.e-6 # micro G
# plottitle = r'''
# \shortstack{Uniform medium, $t = %.2f$ Myr,\\
# $V_\mathrm{max} = %.1f$~km/s, 
# $B = [%.0f, %.0f]~\mu$G}
# ''' % (float(itime)/100.0, vmax, bmin, bmax)
plottitle=None

# Calculate vector potential: A = \int B_y dx - \int B_x dy
# First stab - simply use cumulative sum for the integration
fx = N.cumsum(data['bx'][:,0])[:,N.newaxis] # \int B_x dy
fy = N.cumsum(data['by'], axis=1) # \int B_y dx
vecpot = fy - fx

# Do the graph
import pyxgraph, pyx
pyx.text.set(mode="latex")
pyx.text.preamble(r"""\usepackage{mathpazo}""")
ny, nx = vecpot.shape
figwidth = 10.0
worldheight = worldwidth*ny/nx
figheight = ny*figwidth/nx
gleft = pyxgraph.pyxgraph(xlimits=(0, worldwidth), 
			  ylimits=(0, worldheight),
			  width=figwidth, height=figheight, 
			  key=None,
			  title=plottitle,
			  xlabel='$x$ (pc)',
			  ylabel='$y$ (pc)',
		      )
gright = pyxgraph.pyxgraph(xlimits=(0, worldwidth), 
			   ylimits=(0, worldheight),
			   width=figwidth, height=figheight, 
			   key=None,
			   title=plottitle,
			   xlabel='$x$ (pc)',
			   ylabel=None, ytexter=False,
			   )
# image of the density
# pm3d=[3,3,3] means linear in each channel
# pm3d=[4,4,4] means x**2 in each channel
# pm3d=[7,7,7] means x**0.5 in each channel
graymap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, invert=0, pm3d=[7,7,7])
# mycolmap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, brightness=0.2)
# mycolmap = pyxgraph.ColMapper.ColorMapper("white-yellow-red-black", exponent=0.6, brightness=0.4)
mycolmap = graymap

# g.pyxplotarray(data['pp'], colmap=graymap,
imagedata = N.log10(data['dd']/(1.3*1.67262158e-24)) # log10 number density
imagemin, imagemax = N.log10(10.0), N.log10(1.5e6)
# WJH 05 Apr 2009 - change back to fixed limits on density for consistency between times
# imagemin, imagemax = N.log10(N.floor(10**imagedata.min())), N.log10(100*N.ceil(10**imagedata.max()/100))
gright.pyxplotarray(imagedata[::-1,:], colmap=mycolmap, minvalue=imagemin, maxvalue=imagemax,
	       xpos=0.0, ypos=0.0, width=worldwidth, height=worldheight, graphcoords=True)

cb = pyxgraph.pyxcolorbar(lut=mycolmap.generate_lut(), frame=gright, pos=(1.1,0.0),
			  orientation="vertical2",
			  minlabel="$\log n = %.2f$" % imagemin, 
			  maxlabel="$\log n = %.2f$" % imagemax 
			  )
gright.insert(cb)


# Plot vectors of B-field
# B = N.sqrt(data['bx']**2 + data['by']**2) # magnitude of B-field
# theta = N.arctan2(data['by'], data['bx']) # angle of B with x-axis

# Subtract off mean velocity of neutral gas
vmean = (data['vx']*data['dd']*data['xn']).sum() / (data['dd']*data['xn']).sum()

# No longer subtract this velocity - WJH 05 Apr 2009
# data['vx'] -= vmean 
print "Mean x-velocity of neutral gas: %i km/s" % (vmean/1.e5)
print "But that was not subtracted!"

# Plot contours of vector potential
xx = N.arange(nx) + 0.5
yy = N.arange(ny) + 0.5
gleft.pyxplotcontour(vecpot, (worldwidth/nx)*(xx+1.0), (worldheight/ny)*(yy+1.0), 
		 # add 1 to positions to align with image
		 levels=50, colors='color', color=pyx.color.gray(0.5),
		 lw=2, 
		 lineattrs=[pyx.color.transparency(0.2)]
		 )

# Plot contour at ionization front
gright.pyxplotcontour(data['xn'], (worldwidth/nx)*(xx+0.5), (worldheight/ny)*(yy+0.5), 
		 levels=[0.5], 
# 		 levels=[0.1, 0.5, 0.9], 
		 colors='color', 
		 color=pyx.color.gray.black,
		 lw=1, 
		 lineattrs=[pyx.color.transparency(0.5), pyx.style.linestyle.dashed]
		 )

# Plot vectors of velocity
v = N.sqrt(data['vx']**2 + data['vy']**2)/vnorm # normalized magnitude of velocity
theta = N.arctan2(data['vy'], data['vx']) # angle of v with x-axis
# skip = 8 			# do an arrow every this many pixels
# mx = nx/skip
# my = ny/skip
# mx, my = 64, 64*ny/nx			# fixed grid of arrows, independent of resolution 
mx, my = 40, 40*ny/nx			# fixed grid of arrows, independent of resolution 
skip = nx/mx

# we abuse a parametric function below, so we express everything in
# terms of a parameter k
import random
# allow skip do be non-integer
x = lambda k: random.randint(int(skip)/4,int(3*skip)/4) + (int(k*skip)/my)
y = lambda k: random.randint(int(skip)/4,int(3*skip)/4) + int(skip*(int(k)%my)) - 1
# s = lambda k: 3.e4*B[y(k),x(k)]
# s = lambda k: 0.5			# all arrows the same size 
s = lambda k: v[min(y(k),ny-1),min(x(k),nx-1)]# *nx/worldwidth 
a = lambda k: theta[min(y(k),ny-1),min(x(k),nx-1)]*180/N.pi   
arrowfunc = pyx.graph.data.paramfunction(
    "k", 0, mx*my-1, "x, y, size, angle = (worldwidth/nx)*x(k), (worldheight/ny)*y(k), s(k), a(k)",
    points=mx*my, context=locals())
vectorstyles = [
#     pyx.color.rgb(1.0,0.6,0.0), # orange
#     pyx.color.rgb(0.3,1.0,0.0), # green
#     pyx.color.rgb(1.0,0.4,0.6), # pink
    pyx.color.gray.black,
    pyx.color.transparency(0.1),
    pyx.style.linewidth.Thick,
    ]
# vectorbgstyles = [
#     pyx.color.rgb(0.0,0.0,0.0), 
#     pyx.color.transparency(0.5),
#     pyx.style.linewidth.THick,
#     ]

# arrowscale=1.2
# s = lambda k: v[y(k),x(k)]*arrowscale
# g.plot(arrowfunc, [pyx.graph.style.arrow(arrowsize=0.15*arrowscale,
# 					 arrowattrs=vectorbgstyles,
# 					 lineattrs=vectorbgstyles,
# 					 )])
# arrowscale=1.0
# s = lambda k: v[y(k),x(k)]*arrowscale
gleft.plot(arrowfunc, [pyx.graph.style.arrow(arrowsize=0.15,
					 arrowattrs=vectorstyles,
					 lineattrs=vectorstyles,
					 )])


# WJH 05 Apr 2009 - show evolutionary time
textmargin = 12*pyx.unit.x_pt
textalignment = [pyx.text.halign.right, pyx.text.valign.bottom] 
gright.text(figwidth - textmargin, textmargin, "%i,000 years" % (itime), textalignment)

c = pyx.canvas.canvas()
c.insert(gleft)
margin = 0.05*figwidth
c.insert(gright, [pyx.trafo.translate(figwidth + margin, 0)])

# g.pyxplot(data=((nx/2,),(ny/2,)), style="points", color="green", ps=5, pt=0)

# Write file
c.writePDFfile('%s-%s-%i' % (execname.split('.')[0], runid, itime))
