"""
Performs circular shift on the datacubes of all variables and writes them out with a new name.
Currently, only does shift by half of the grid size along each axis. Motivation is to test that 
the periodic boundary conditions are working properly. 

v 0.1: WJH 16 Feb 2010 
"""

import pyfits, numpy
import os, sys

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try:
    runid = str(sys.argv[1]) 	# name of run
    itime = int(sys.argv[2])	# save time
    newrunid = str(sys.argv[3]) # new name for cshifted cubes
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME NEWRUNID" % execname
    exit

varlist = ["bx", "by", "bz", "vx", "vy", "vz", "dd", "xn", "xi", "pp"]

for var in varlist:
    prefix = '%s-%s%4.4i' % (runid, var, itime)
    # read in variable at t=itime
    data = pyfits.open(prefix + '.fits')[0].data

    # circular shift by half cube along each axis
    nz, ny, nx = data.shape
    data = numpy.roll(numpy.roll(numpy.roll(data, nx/2, axis=2), ny/2, axis=1), nz/2, axis=0)

    # write out the new fits file
    hdu = pyfits.PrimaryHDU()
    hdu.data = data
    hdu.writeto('%s-%s%4.4i.fits' % (newrunid, var, itime), clobber=True)

