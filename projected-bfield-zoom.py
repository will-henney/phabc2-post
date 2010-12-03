from __future__ import division
import os, sys, argparse
import pyfits
# Parse command line arguments
parser = argparse.ArgumentParser(description="""
Graph the projected B field from a certan viewing angle, as calculated by makerotbmaps

Line-of-sight field rotation measure is shown as grayscale, while plane-of-sky field is shown by arrows. 
""")
parser.add_argument("--runid", "-r", type=str, default='Ostar-et', 
                    help='ID for model run')
parser.add_argument("--itime", "-i", type=int, default=280, 
                    help='Integer save time counter')
parser.add_argument("--compact", action="store_true",
                    help='Use compact labels on panels')
parser.add_argument("--arrows", action="store_true", default=False,
                    help='Plot signed vector integrated B with arrows (otherwise use unsigned vectors)')
parser.add_argument("--theta", type=int, default=0, 
                    help='Viewing angle theta in degrees')
parser.add_argument("--phi", type=int, default=0, 
                    help='Viewing angle phi in degrees')
parser.add_argument("--timestep", type=float, default=1000.0, 
                    help='Length of time in years between saves')
parser.add_argument("--boxsize", "-b", type=float, default=4.0, 
                    help="Size of simulation cube in parsecs")
parser.add_argument("--zoom", "-z", type=float, default=4.0, 
                    help="Zoom factor")
parser.add_argument("--xcenter", "-x", type=int, default=100, 
                    help="X coordinate of center of zoomed area")
parser.add_argument("--ycenter", "-y", type=int, default=50, 
                    help="Y coordinate of center of zoomed area")
args = parser.parse_args()              # we can now use args.runid, args.itime


import numpy as N

def load_fits(id):
    return pyfits.open('%s-bmap-%s-rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, id, args.theta, args.phi, args.itime))['PRIMARY'].data

# Load in the required B-component maps
if args.arrows:
    # int B_i n dz  for i = x, y, z
    XYZ = "xyz"
else:
    # int |B_i| n dz  for i = x, y; int B_z n dz
    XYZ = "UVz"

bxi, byi, bzi = [load_fits("i-%s" % (xyz)) for xyz in XYZ] # ionized
bxn, byn, bzn = [load_fits("n-%s" % (xyz)) for xyz in XYZ] # neutral
bxm, bym, bzm = [load_fits("m-%s" % (xyz)) for xyz in XYZ] # neutral

# Filenames for column densities - used later in the RGB images
colm, coln, coli = ['%s-colmap-%s-rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, mni, args.theta, args.phi, args.itime)
                    for mni in "mni"]

# The actual data for the molecular column densities - used in contours
mcol = pyfits.open(colm)["PRIMARY"].data

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

if args.runid.startswith("Ostar"):
    bi_scale, bn_scale, bm_scale = 0.09, 8.0, 12.0
elif args.runid.startswith("Bstar"):
    bi_scale, bn_scale, bm_scale = 0.002, 1.5, 0.7
else: 
    raise NotImplementedError, "Only Ostar and Bstar allowed so far"

worldwidth, worldheight = args.boxsize/args.zoom, args.boxsize/args.zoom

## functions needed for vector plots
ny, nx = bxi.shape
nx /= args.zoom
ny /= args.zoom
i1 = args.xcenter - nx/2 
i2 = i1 + nx
j1 = args.ycenter - ny/2
j2 = j1 + ny
print "[%d %d %d %d]" % (i1, j1, i2, j2)
bxi = bxi[j1:j2, i1:i2]
byi = byi[j1:j2, i1:i2]
bzi = bzi[j1:j2, i1:i2]
bxn = bxn[j1:j2, i1:i2]
byn = byn[j1:j2, i1:i2]
bzn = bzn[j1:j2, i1:i2]
bxm = bxm[j1:j2, i1:i2]
bym = bym[j1:j2, i1:i2]
bzm = bzm[j1:j2, i1:i2]
mcol = mcol[j1:j2, i1:i2]

# we abuse a parametric function below, so we express everything in
# terms of a parameter k
import random
x = lambda k: int(k)//ny + random.gauss(0.0, 0.25)
y = lambda k: int(k)%nx + random.gauss(0.0, 0.25)
# x = lambda k: int(k)//ny + (random.random() - 0.5) 
# y = lambda k: int(k)%nx + (random.random() - 0.5)
arrow_adjust = 1.0
def s(k):
    jj = min(y(k),ny-1)
    ii = min(x(k),nx-1)
    return arrow_adjust*b[jj,ii]
a = lambda k: theta[min(y(k),ny-1),min(x(k),nx-1)]*180/N.pi   

c = pyx.canvas.canvas()

if args.compact:
    title_i, title_n, title_m = [r'\(\int n_\mathrm{%s} \vec{B} \, dz\)' 
                                 % (sp) for sp in 'ion', 'neut', 'mol']
else:
    title_i, title_n, title_m = [r'\(B\)-Field weighted by %s gas: \(\int n_\mathrm{%s} \vec{B} \, dz\)' 
                                 % (species, sp) 
                                 for species, sp in 
                                 ('ionized', 'ion'), ('neutral', 'neut'), ('molecular', 'mol')]

def add_text(c, s, fmt="%s"):
    "Put the text label 's' inside the box of canvas 'c'"
    x, y = [12*pyx.unit.x_pt,] * 2
    text = pyx.text.text(x, y, fmt % (s))
    tpath = text.bbox().enlarged(3*pyx.unit.x_pt).path()
    c.fill(tpath, [pyx.color.gray(0.5), pyx.color.transparency(0.4)])
    c.insert(text)

# These are used by the contour routine
xx = N.arange(nx) + 0.5
yy = N.arange(ny) + 0.5
xx = (worldwidth/nx)*xx
yy = (worldheight/ny)*yy
# xx = (worldwidth/nx)*(xx+1.0)
# yy = (worldheight/ny)*(yy+1.0)
contourlevels = 3

for title, bx, by, bz, bscale, xshift, yshift in [
    [title_i, bxi, byi, bzi, bi_scale, 0, figwidth + margin],
    [title_n, bxn, byn, bzn, bn_scale, figwidth + margin, figwidth + margin],
    [title_m, bxm, bym, bzm, bm_scale, figwidth + margin, 0],
    ]:
    g = pyxgraph.pyxgraph(
        xlimits = (0, args.boxsize/args.zoom), ylimits = (0, args.boxsize/args.zoom), 
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
                   xpos=0.0, ypos=0.0, width=args.boxsize/args.zoom, height=args.boxsize/args.zoom, graphcoords=True)

    # Plot vectors of plane-of-sky field
    b = N.sqrt(bx**2 + by**2)/bscale # normalized magnitude of B-field
    b.clip(0.0, 1.0, out=b)       # make sure we don't get arrows too big

    theta = N.arctan2(by, bx) # angle of B with x-axis
    arrowfunc = pyx.graph.data.paramfunction(
        "k", 0, nx*ny-1, "x, y, size, angle = (worldwidth/nx)*x(k), (worldheight/ny)*y(k), s(k), a(k)",
        points=nx*ny, context=locals())
    vectorstyles = [
        pyx.color.rgb.red,
        pyx.color.transparency(0.5),
        pyx.style.linewidth.thick,
        ]

    if args.arrows: 
        arrowsize = 0.15
    else:
        arrowsize = 0.0                 # no arrowhead

    g.plot(arrowfunc, [pyx.graph.style.arrow(arrowsize=0.15,
                                             arrowattrs=vectorstyles,
                                             lineattrs=vectorstyles,
                                             )])

    # Add contours of molecular column
    g.pyxplotcontour(mcol, xx, yy,
                     levels=contourlevels, colors='color', color=pyx.color.rgb.blue,
                     lw=2, 
                     lineattrs=[pyx.color.transparency(0.5)]
                     )
    
    add_text(g, title)
    c.insert(g, [pyx.trafo.translate(xshift, yshift)])

# Now do RGB map of the column densities
from myfitsutils import RGB3Image
from PIL import Image
if args.runid.startswith("Ostar"):
    maxcolm, maxcoln, maxcoli = [8.e5, 6.e5, 4.e4] # these are funny units - should multiply by cell size
elif args.runid.startswith("Bstar"):
    maxcolm, maxcoln, maxcoli = [6.e5, 4.e5, 6.e3]

rgbim = RGB3Image(
    redfile=colm, greenfile=coln, bluefile=coli, 
    redmin = 0.0, greenmin = 0.0, bluemin = 0.0,
    redmax = maxcolm, greenmax = maxcoln, bluemax = maxcoli,
    bb=[i1, j1, i2, j2]
    )
print rgbim.red.shape
print rgbim.red.bb
rgbim = rgbim.transpose(Image.FLIP_TOP_BOTTOM)
g = pyxgraph.pyxgraph(
    xlimits = (0, args.boxsize/args.zoom), ylimits = (0, args.boxsize/args.zoom), 
    width = figwidth, height = figheight,
    key = None,
    title = None, xlabel = None, ylabel = None, 
    # title = title,
    # xlabel=r'$x$ (pc)',
    # ylabel=r'$y$ (pc)',
    xtexter = False, ytexter = False,
    )

g.insert(pyx.bitmap.bitmap(0, 0, rgbim, width=figwidth))
# g.pyxplotcontour(mcol, xx, yy, 
#                  levels=contourlevels, colors='color', color=pyx.color.rgb.green,
#                  lw=2, 
#                  lineattrs=[pyx.color.transparency(0.5)]
#                  )

if args.compact:
    title_prefix = ""
else:
    title_prefix = r"RGB column densities: "
add_text(g, 
         title_prefix + 
         ",".join([r"\color{%s}{\(\int n_\mathrm{%s}\, dz\)}" % (col, sp) 
                   for col, sp in ["red", "mol"], ["green", "neut"], ["blue", "ion"]]),
         fmt=r"\color{white}{%s}")

c.insert(g)

# Write out to PDF file
execroot = os.path.split(sys.argv[0])[-1].split('.')[0]

c.writePDFfile("%s-%s-%i-rot%+3.3i%+3.3i-x%3.3i-y%3.3i" 
               % (execroot, args.runid, args.itime, args.theta, args.phi, args.xcenter, args.ycenter))
