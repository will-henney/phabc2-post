from __future__ import division
import os, sys, argparse
import pyfits
import pyxgraph, pyx
import warnings

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
                    help='Plot signed vector-integrated B with arrows (otherwise plot "Stokes"-integrated B with lines)')
parser.add_argument("--scale", type=str, default="runid", choices=["runid", "max"],
                    help='Policy for mapping onto brightness scale')
parser.add_argument("--gamma", type=float, default=1.0, 
                    help="Gamma correction for RGB image. Values > 1 make image brighter")
parser.add_argument("--bscale", type=float, default=1.0, 
                    help="Extra scale factor for the B-field (<1 emphasizes faint structures)")
parser.add_argument("--molecular-bscale", type=float, default=1.0, 
                    help="Extra extra scale factor for the molecular B-field (<1 emphasizes faint structures)")
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
parser.add_argument("--xcenter", "-x", type=int, default=-1, 
                    help="X coordinate of center of zoomed area. If < 0, grid center is used")
parser.add_argument("--ycenter", "-y", type=int, default=-1, 
                    help="Y coordinate of center of zoomed area. If < 0, grid center is used")
parser.add_argument("--vector-density", "-d", type=float, default=1.0, 
                    help="Density of vectors to plot (<= 1.0). A value of 1.0 means plot every point, while with values < 1.0 a fraction of the points are chosen randomly")
parser.add_argument("--jitter", "-j", type=float, default=0.3, 
                    help="Amplitude of noise added to x,y values for vector plotting - prevents overlap of vectors")

args = parser.parse_args()              # we can now use args.runid, args.itime


import numpy as N

def load_fits(id):
    return pyfits.open('%s-bmap-%s-rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, id, args.theta, args.phi, args.itime))['PRIMARY'].data

# Load in the required B-component maps
vectorstyles = [
    pyx.color.rgb.red,
    pyx.color.transparency(0.5),
    pyx.style.linewidth.thick,
    ]
arrowattrs = vectorstyles
if args.arrows:
    # int B_i n dz  for i = x, y, z
    XYZ = "xyz"
    arrowsize = 0.15
    linelength = 0.25*pyx.unit.v_cm
    arrows_id = "-vec"
else:
    # int |B_i| n dz  for i = x, y; int B_z n dz
    XYZ = "UVz"
    arrowsize = None
    linelength = 0.5*pyx.unit.v_cm
    arrows_id = "-abs"
    

bxi, byi, bzi = [load_fits("i-%s" % (xyz)) for xyz in XYZ] # ionized
bxn, byn, bzn = [load_fits("n-%s" % (xyz)) for xyz in XYZ] # neutral
bxm, bym, bzm = [load_fits("m-%s" % (xyz)) for xyz in XYZ] # neutral

# Filenames for column densities - used later in the RGB images
colm, coln, coli = ['%s-colmap-%s-rot%+3.3i%+3.3i-%4.4i.fits' 
                       % (args.runid, mni, args.theta, args.phi, args.itime)
                    for mni in "mni"]

# The actual data for the molecular column densities - used in contours
mcol = pyfits.open(colm)["PRIMARY"].data
ncol = pyfits.open(coln)["PRIMARY"].data
icol = pyfits.open(coli)["PRIMARY"].data

# PyX gives lots of unnecessary warnings - let's turn them off!
warnings.simplefilter("ignore")

# Do the graph
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


if args.scale == "runid":
    if args.runid.startswith("Ostar"):
        bi_scale, bn_scale, bm_scale = 0.09, 8.0, 12.0
    elif args.runid.startswith("Bstar"):
        bi_scale, bn_scale, bm_scale = 0.002, 1.5, 0.7
    elif args.runid.startswith("s"):
        bi_scale, bn_scale, bm_scale = 1.5, 100.0, 100.0
    elif args.runid.startswith("w"):
        bi_scale, bn_scale, bm_scale = 0.5, 30.0, 30.0
    else: 
        raise NotImplementedError, "Unimplemented runid-based scaling for this type of run"
elif args.scale == "max":
    bi_scale = max([b.max() for b in bxi, byi, bzi])
    bn_scale = max([b.max() for b in bxn, byn, bzn])
    bm_scale = max([b.max() for b in bxm, bym, bzm])

# extra fiddle factor
bi_scale *= args.bscale
bn_scale *= args.bscale
bm_scale *= args.bscale*args.molecular_bscale

## functions needed for vector plots
ny, nx = bxi.shape
parsec = 3.085677582e18
# if args.runid.startswith("Ostar") or args.runid.startswith("Bstar"):
#     boxsize_pc = 4.0
# elif args.runid.startswith("s") or args.runid.startswith("w"):
#     boxsize_pc = 1.0
# else:
#     raise NotImplementedError, "Unknown boxsize for this type of run"

cellsize = parsec*args.boxsize/ny

nx //= args.zoom
ny //= args.zoom
# print nx, ny
worldwidth, worldheight = (nx/ny)*args.boxsize/args.zoom, args.boxsize/args.zoom
figheight *= float(ny)/nx

# Default value for [xy]center is -1, which has special meaning: we convert to the center of the input grid
if args.xcenter < 0: 
    args.xcenter = nx/2
if args.ycenter < 0: 
    args.ycenter = ny/2

# Calculate limits of zoomed patch in input arrays
i1 = args.xcenter - nx/2 
i2 = i1 + nx
j1 = args.ycenter - ny/2
j2 = j1 + ny
# print "[%d %d %d %d]" % (i1, j1, i2, j2)
# Crop data to zoomed patch
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
ncol = ncol[j1:j2, i1:i2]
icol = icol[j1:j2, i1:i2]

# 04 Dec 2010 Redo vector field as data arrays, not parametric function
y, x = N.mgrid[0:ny, 0:nx]
x = x.astype(N.float)
y = y.astype(N.float)
# add some jitter to the positions to help with the case of vector_density = 1.0
y += N.random.normal(scale=args.jitter, size=[ny, nx])
x += N.random.normal(scale=args.jitter, size=[ny, nx])

# decimate the points to plot for vector_density < 1
indices = range(ny*nx)
N.random.shuffle(indices)
indices = indices[:int(ny*nx*args.vector_density)]

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
    [title_i, bxi, byi, bzi, bi_scale, 0, figheight + margin],
    [title_n, bxn, byn, bzn, bn_scale, figwidth + margin, figheight + margin],
    [title_m, bxm, bym, bzm, bm_scale, figwidth + margin, 0],
    ]:
    g = pyxgraph.pyxgraph(
        xlimits = (0, worldwidth), ylimits = (0, worldheight), 
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
                   xpos=0.0, ypos=0.0, width=worldwidth, height=worldheight, graphcoords=True)

    # Plot vectors of plane-of-sky field
    b = N.sqrt(bx**2 + by**2)/bscale # normalized magnitude of B-field
    b.clip(0.0, 1.0, out=b)       # make sure we don't get arrows too big

    theta = N.arctan2(by, bx) # angle of B with x-axis
    vectordata = pyx.graph.data.values(
        x = (worldwidth/nx)*x.ravel()[indices], y = (worldheight/ny)*y.ravel()[indices], 
        size = b.ravel()[indices], angle = theta.ravel()[indices]*180.0/N.pi
        )
    g.plot(vectordata, [pyx.graph.style.arrow(arrowsize=arrowsize,
                                              arrowattrs=arrowattrs,
                                              lineattrs=vectorstyles,
                                              linelength=linelength
                                              )])

    # Add contours of molecular column
    g.pyxplotcontour(mcol, xx, yy,
                     levels=contourlevels, colors='color', color=pyx.color.rgb.blue,
                     lw=2, 
                     lineattrs=[pyx.color.transparency(0.5)]
                     )
    g.finish()
    add_text(g, title)
    c.insert(g, [pyx.trafo.translate(xshift, yshift)])

# Now do RGB map of the column densities
from myfitsutils import RGB3Image
from PIL import Image


if args.scale == "runid":
    if args.runid.startswith("Ostar"):
        maxcolm, maxcoln, maxcoli = [8.e5, 6.e5, 4.e4] # these are funny units - should multiply by cell size
    elif args.runid.startswith("Bstar"):
        maxcolm, maxcoln, maxcoli = [6.e5, 4.e5, 6.e3]
    elif args.runid.startswith("s") or args.runid.startswith("w"):
        maxcolm, maxcoln, maxcoli = [6.e6, 2.e6, 2.e5]
elif args.scale == "max":
     maxcolm, maxcoln, maxcoli = None, None, None

rgbim = RGB3Image(
    redfile=colm, greenfile=coln, bluefile=coli, 
    redmin = 0.0, greenmin = 0.0, bluemin = 0.0,
    redmax = maxcolm, greenmax = maxcoln, bluemax = maxcoli,
    bb=[i1, j1, i2, j2], gamma=args.gamma,
    )
# need to recover these here (in case they were None) since the transpose destroys access to the channels
maxcolm, maxcoln, maxcoli = rgbim.red.max, rgbim.green.max, rgbim.blue.max 
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

savefile = "%s-%s-%s-%i-rot%+3.3i%+3.3i-x%3.3i-y%3.3i%s" \
               % (execroot, args.zoom, args.runid, args.itime, args.theta, args.phi, args.xcenter, args.ycenter, arrows_id)

c.writePDFfile(savefile)

print "Written graph to [[%s.pdf]]" % (savefile)

sigma_dust = 5e-22
AV_per_tau = 1.08573620476
microgauss = 1.e-6/N.sqrt(4.0*N.pi)
# print "Maximum column densities [ion/neut/mol] = [%.1e, %.1e, %.1e] mag V-band" % tuple(
#     [m*cellsize*sigma_dust*AV_per_tau for m in maxcoli, maxcoln, maxcolm])
print "Maximum column densities [ion/neut/mol] = [%.1e, %.1e, %.1e] /cm^2" % tuple(
    [m*cellsize for m in maxcoli, maxcoln, maxcolm])
print "Maximum integrated B [ion/neut/mol] = [%.1e, %.1e, %.1e] micro G /cm^2" % tuple(
    # [m*cellsize*sigma_dust*AV_per_tau/microgauss for m in bi_scale, bn_scale, bm_scale]
    [m*cellsize/microgauss for m in bi_scale, bn_scale, bm_scale])
print "Mean B fields [ion/neut/mol] =  [%.1e, %.1e, %.1e] micro G" % tuple(
    [N.mean(N.sqrt(bx**2 + by**2 + bz**2))/N.mean(colden)/microgauss 
     for colden, bx, by, bz in [icol, bxi, byi, bzi], [ncol, bxn, byn, bzn], [mcol, bxm, bym, bzm]])
