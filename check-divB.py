"""
Calculate Div B for PhabC2 data, writing results to a FITS file. 

Uses same discretization that Fabio uses. 
"""
import pyfits
import os, sys
import numpy as N

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try:
    runid = str(sys.argv[1])
    itime = int(sys.argv[2])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME" % execname
    exit

bx, by, bz = [ pyfits.open('%s-%s%4.4i.fits' % (runid, varid, itime))[0].data
	       for varid in 'bx', 'by', 'bz' ]

if len(bx.shape) == 2: 
    # promote 2d array to 3d array
    bx = N.array([bx, bx, bx])
    by = N.array([by, by, by])
    bz = N.array([bz, bz, bz])

print bx.shape
divbx = N.zeros(bx.shape)
divby = N.zeros(bx.shape)
divbz = N.zeros(bx.shape)


# Note: This only works if the BCs are periodic
# Remember, x-axis is the last one in pyfits!

# WJH 08 May 2010 - rw
# Bx(i+1,j,k) - Bx(i-1,j,k)
divbx = N.roll(bx, -1, axis=2) - N.roll(bx, 1, axis=2)
# By(i,j+1,k) - By(i,j-1,k)
divby = N.roll(by, -1, axis=1) - N.roll(by, 1, axis=1)
# Bz(i,j,k+1) - Bz(i,j,k-1)
divbz = N.roll(bz, -1, axis=0) - N.roll(bz, 1, axis=0)

# divbx[:,:,1:-1] = bx[:,:,2:] -  bx[:,:,:-2]
# divbx[:,:,0] = bx[:,:,1] - bx[:,:,-1]
# divbx[:,:,-1] = bx[:,:,0] - bx[:,:,-2]
# divby[:,1:-1,:] = by[:,2:,:] -  by[:,:-2,:]
# divby[:,0,:] = by[:,1,:] - by[:,-1,:]
# divby[:,-1,:] = by[:,0,:] - by[:,-2,:]
# divbz[1:-1,:,:] = bz[2:,:,:] -  bz[:-2,:,:]
# divbz[0,:,:] = bz[1,:,:] - bz[-1,:,:]
# divbz[-1,:,:] = bz[0,:,:] - bz[-2,:,:]

divb = divbx + divby + divbz

b = N.sqrt(bx**2 + by**2 + bz**2)

# normalised divergence
divb_b = divb/b

hdu = pyfits.PrimaryHDU()
hdu.data = divb
hdu.writeto('%s-%s%4.4i.fits' % (runid, 'db', itime), clobber=True)
hdu.data = divb_b
hdu.writeto('%s-%s%4.4i.fits' % (runid, 'nd', itime), clobber=True)


print 'h Div B : mean || %3.3e min %3.3e max %3.3e' % (
    N.abs(divb).mean(), divb.min(), divb.max())

print 'h Div B / |B|: mean || %3.3e min %3.3e max %3.3e' % (
    N.abs(divb_b).mean(), divb_b.min(), divb_b.max())
