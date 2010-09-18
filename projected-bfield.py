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
    return pyfits.open('%s-bmap-%s-rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, id, args.theta, args.phi, args.itime))['PRIMARY'].data

# Load in the required B-component maps
bxi, byi, bzi = [load_fits("i-%s" % (xyz)) for xyz in "xyz"] # ionized
bxn, byn, bzn = [load_fits("n-%s" % (xyz)) for xyz in "xyz"] # neutral
bxm, bym, bzm = [load_fits("m-%s" % (xyz)) for xyz in "xyz"] # neutral


# Do the graph
import pyxgraph, pyx
pyx.text.set(mode="latex")
pyx.text.preamble(r"\usepackage{txfonts,color}\AtBeginDocument{\bfseries\boldmath}")
figwidth, figheight, margin = 10, 10, 0.0

# pm3d=[3,3,3] means linear in each channel
# pm3d=[4,4,4] means x**2 in each channel
# pm3d=[7,7,7] means x**0.5 in each channel
graymap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, invert=0, pm3d=[7,7,7])
# mycolmap = pyxgraph.ColMapper.ColorMapper("pm3d", exponent=1.0, brightness=0.2)
# mycolmap = pyxgraph.ColMapper.ColorMapper("white-yellow-red-black", exponent=0.6, brightness=0.4)
mycolmap = graymap

bi_scale, bn_scale, bm_scale = 0.09, 8.0, 12.0
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

title_i, title_n, title_m = [r'\(B\)-Field weighted by %s gas: \(\int n_\mathrm{%s} \vec{B} \, dz\)' 
                             % (species, sp) 
                             for species, sp in ('ionized', 'ion'), ('neutral', 'neut'), ('molecular', 'mol')
                            ]

def add_text(c, s, fmt="%s"):
    "Put the text label 's' inside the box of canvas 'c'"
    x, y = [6*pyx.unit.x_pt,] * 2
    c.text(x, y, fmt % (s))
    
for title, bx, by, bz, bscale, xshift, yshift in [
    [title_i, bxi, byi, bzi, bi_scale, 0, figwidth + margin],
    [title_n, bxn, byn, bzn, bn_scale, figwidth + margin, figwidth + margin],
    [title_m, bxm, bym, bzm, bm_scale, figwidth + margin, 0],
    ]:
    g = pyxgraph.pyxgraph(
        xlimits = (0, args.boxsize), ylimits = (0, args.boxsize), 
        width = figwidth, height = figheight,
        key = None,
        title = None, xlabel = None, ylabel = None, 
        # title = title,
        # xlabel=r'$x$ (pc)',
        # ylabel=r'$y$ (pc)',
        xtexter = False, ytexter = False,
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
    add_text(g, title)
    c.insert(g, [pyx.trafo.translate(xshift, yshift)])

# Now read in the column densities
from myfitsutils import RGB3Image
from PIL import Image
colm, coln, coli = ['%s-colmap-%s-rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, mni, args.theta, args.phi, args.itime)
                    for mni in "mni"]
maxcolm, maxcoln, maxcoli = [8.e5, 6.e5, 4.e4] # these are funny units - should multiply by cell size
rgbim = RGB3Image(
    redfile=colm, greenfile=coln, bluefile=coli, 
    redmin = 0.0, greenmin = 0.0, bluemin = 0.0,
    redmax = maxcolm, greenmax = maxcoln, bluemax = maxcoli,
    )
rgbim = rgbim.transpose(Image.FLIP_TOP_BOTTOM)
c.insert(pyx.bitmap.bitmap(0, 0, rgbim, width=figwidth))
add_text(c, 
         r"RGB column densities: " + 
         ",".join([r"\color{%s}{\(\int n_\mathrm{%s}\, dz\)}" % (col, sp) 
                   for col, sp in ["red", "mol"], ["green", "neut"], ["blue", "ion"]]),
         fmt=r"\color{white}{%s}")

# Write out to PDF file
execroot = os.path.split(sys.argv[0])[-1].split('.')[0]

c.writePDFfile("%s-%s-%i-rot%+3.3i%+3.3i" % (execroot, args.runid, args.itime, args.theta, args.phi))
