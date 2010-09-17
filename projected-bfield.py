from __future__ import division
import os, sys, argparse
# Parse command line arguments
parser = argparse.ArgumentParser(description="""
Graph the projected B field from a certan viewing angle, as calculated by makerotbmaps

Line-of-sight field rotation measure is shown as grayscale, while plane-of-sky field is shown by arrows. 
""")
parser.add_argument("--runid", "-r", type=str, default='Ostar-et', 
                    help='ID for model run')
parser.add_argument("--itime", "-i", type=int, default=280, 
                    help='Integer save time counter')
parser.add_argument("--theta", type=int, default=0, 
                    help='Viewing angle theta in degrees')
parser.add_argument("--phi", type=int, default=0, 
                    help='Viewing angle phi in degrees')
parser.add_argument("--timestep", type=float, default=1000.0, 
                    help='Length of time in years between saves')
parser.add_argument("--boxsize", "-b", type=float, default=4.0, 
                    help="Size of simulation cube in parsecs")
args = parser.parse_args()              # we can now use args.runid, args.itime


import numpy as N

def load_fits(id):
    import pyfits
    return pyfits.open('%s-bmap-%s--rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, id, args.theta, args.phi, args.itime))['PRIMARY'].data

# Load in the required B-component maps
bxi, byi, bzi = [load_fits("i-%s" % (xyz)) for xyz in "xyz"] # ionized
bxn, byn, bzn = [load_fits("n-%s" % (xyz)) for xyz in "xyz"] # neutral


# Do the graph
import pyxgraph, pyx

figwidth, figheight, margin = 10, 10, 2

# pm3d=[3,3,3] means linear in each channel
# pm3d=[4,4,4] means x**2 in each channel
# pm3d=[7,7,7] means x**0.5 in each channel
graymap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, invert=0, pm3d=[7,7,7])
# mycolmap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, brightness=0.2)
# mycolmap = pyxgraph.ColMapper.ColorMapper("white-yellow-red-black", exponent=0.6, brightness=0.4)
mycolmap = graymap

bi_scale, bn_scale = 0.08, 15.0
worldwidth, worldheight = args.boxsize, args.boxsize

## functions needed for vector plots
ny, nx = bxi.shape
mx, my = 64, 64*ny/nx			# fixed grid of arrows, independent of resolution 
skip = nx/mx
# we abuse a parametric function below, so we express everything in
# terms of a parameter k
import random
# allow skip do be non-integer
x = lambda k: random.randint(int(skip)//4,int(3*skip)//4) + (int(k*skip)//my)
y = lambda k: random.randint(int(skip)//4,int(3*skip)//4) + int(skip*(int(k)%my)) - 1
arrow_adjust = 2.0
def s(k):
    jj = min(y(k),ny-1)
    ii = min(x(k),nx-1)
    return arrow_adjust*b[jj,ii]
a = lambda k: theta[min(y(k),ny-1),min(x(k),nx-1)]*180/N.pi   

c = pyx.canvas.canvas()

title_i = r'Rotation measure: \(\int n_\mathrm{e} B_\parallel \, dz\)'
title_n = r'Rotation measure: \(\int n_\mathrm{n} B_\parallel \, dz\)'

for title, bx, by, bz, bscale, xshift in [
    [title_i, bxi, byi, bzi, bi_scale, 0],
    [title_n, bxn, byn, bzn, bn_scale, figwidth + margin],
    ]:
    g = pyxgraph.pyxgraph(
        xlimits = (0, args.boxsize), ylimits = (0, args.boxsize), 
        width = figwidth, height = figheight,
        key = None,
        title = title,
        xlabel=r'$x$ (pc)',
        ylabel=r'$y$ (pc)',
        )

    # Grayscale of z component
    g.pyxplotarray(bz[::-1,:], colmap=mycolmap, minvalue=-bscale, maxvalue=bscale,
                   xpos=0.0, ypos=0.0, width=args.boxsize, height=args.boxsize, graphcoords=True)

    # Plot vectors of plane-of-sky field
    b = N.sqrt(bx**2 + by**2)/bscale # normalized magnitude of B-field
    b.clip(0.0, 1.0, out=b)       # make sure we don't get arrows too big

    theta = N.arctan2(by, bx) # angle of B with x-axis
    arrowfunc = pyx.graph.data.paramfunction(
        "k", 0, mx*my-1, "x, y, size, angle = (worldwidth/nx)*x(k), (worldheight/ny)*y(k), s(k), a(k)",
        points=mx*my, context=locals())
    vectorstyles = [
        pyx.color.rgb.red,
        pyx.color.transparency(0.5),
        pyx.style.linewidth.thick,
        ]

    g.plot(arrowfunc, [pyx.graph.style.arrow(arrowsize=0.15,
                                             arrowattrs=vectorstyles,
                                             lineattrs=vectorstyles,
                                             )])
    c.insert(g, [pyx.trafo.translate(xshift, 0)])

# Write out to PDF file
execroot = os.path.split(sys.argv[0])[-1].split('.')[0]

c.writePDFfile("%s-%s-%i-rot%+3.3i%+3.3i" % (execroot, args.runid, args.itime, args.theta, args.phi))
