import pyfits
import numpy as N

import os, sys

# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
except IndexError, ValueError:
    print "Usage: %s RUNID" % execname
    exit


def get_fits_data(varid):
    return pyfits.open('%s-%s%4.4i.fits' % (runid, varid, itime))['PRIMARY'].data

maxtime = 1000			# maximum possible time to look for
for itime in range(maxtime): 
    dd = get_fits_data('dd')    
    pp = get_fits_data('pp')
    bx = get_fits_data('bx')
    by = get_fits_data('by')
    bz = get_fits_data('bz')
    vx = get_fits_data('vx')
    vy = get_fits_data('vy')
    vz = get_fits_data('vz')
    xi = get_fits_data('xi')

    eb = 0.5*(bx**2 + by**2 + bz**2)
    pt = 0.5*dd*(vx**2 + vy**2 + vz**2)
    bb = N.sqrt(4.0*N.pi*(bx**2 + by**2 + bz**2))
    dn = dd / (1.3*1.67262158e-24)


