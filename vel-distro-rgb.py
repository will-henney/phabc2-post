"""
RGB plots of radial velocity distributions as a function of radius

Color represents ion/molecular state. 

This program originally inspired 23 Aug 2010 by puzzling results from
stats programs, where the molecular gas seems to have negative
velocities

"""
import bivar                            # my bivariate distribution plotting library

import os, sys, argparse

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
parser = argparse.ArgumentParser(description="Radial velocity as a function of radius")
parser.add_argument("--runid", "-r", type=str, default='Bstar-ep', 
                    help='ID for model run')
parser.add_argument("--itime", "-i", type=int, default=1, 
                    help='Integer save time counter')
parser.add_argument("--timestep", type=float, default=None, 
                    help='Length of time in years between saves')
parser.add_argument("--boxsize", "-b", type=float, default=4.0, 
                    help="Size of simulation cube in parsecs")
args = parser.parse_args()              # we can now use args.runid, args.itime


# Heuristic to identify timestep if not specified
if args.timestep is None:
    raise NotImplementedError("Sorry, --timestep must be given explicitly for now")

print args

def load_fits(id):
    import pyfits
    return pyfits.open('%s-%s%4.4i.fits' % (args.runid, id, args.itime))['PRIMARY'].data

# Load in the required data cubes
vx, vy, vz, xn, AV = [load_fits(id) for id in ["vx", "vy", "vz", "xn", "AV"]]

import numpy
# make cube of unit radius vector
nx, ny, nz = vx.shape                          # grid shape
x, y, z, = numpy.ogrid[0:nx, 0:ny, 0:nz]       # grid coordinates
xc, yc, zc = [n/2 - 0.5 for n in [nx, ny, nz]] # grid center
# grid of scalar radius from center (in units of cell size)
rr = numpy.sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2)
# components of unit radius vector from center
ux, uy, uz = (x - xc) / rr, (y - yc) / rr, (z - zc) / rr

# rescale radius to be in parsecs
rr *= args.boxsize/nx

# radial velocity is dot product of vector velocity with unit radius vector
vr = vx*ux + vy*uy + vz*uz

# rescale velocities to be in km/s
vr /= 1.e5

import molfrac
xmol = molfrac.molfrac(AV)

# List of weights will produce RGB pdf images
weights = [
    xn*xmol,                            # molecular
    xn*(1.0-xmol),                      # neutral
    1.0-xn                              # ionized
    ]


# Now make the plots
import pyx
pyx.text.set(mode="latex")
# pyx.text.preamble("""\usepackage{mathptmx}""")
bivar.printextra = 0
bivar.Graph.figwidth = 6
bivar.Graph.figheight = 6
bivar.PlotVariable.n = 50                               # size of pdf images

rrvar = bivar.PlotVariable(rr)          # radius
rmax = rr.max()
rrvar.setminmaxn(min=0.0, max=rmax)
rrvar.settitle(r'Radius, \(r\), pc', 'Radius')

vrvar = bivar.PlotVariable(vr, drawzero=1)          # radial velocity
vrmin, vrmax = -20.0, 20.0
vrvar.setminmaxn(min=vrmin, max=vrmax)
vrvar.settitle(r'Radial velocity, \(V_r\), km s\(^{-1}\)', 'V_r')

g = bivar.Graph(rrvar, vrvar, weights=weights, gamma=2.0, statslevel=1)

# add a title
gtitle = r'\texttt{%s}, $t = %.3f$~Myr' % (args.runid, args.itime*args.timestep/1.e6)
g.text(0.5*g.figwidth, g.figheight+6*pyx.unit.x_pt, gtitle, [pyx.text.halign.center])

# and save
g.writePDFfile('%s-%s-%4.4i' % (execname.split('.')[0], args.runid, args.itime))
