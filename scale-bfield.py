"""
Copies a set of datacubes to a new name while multiplying all B components by SCALE
"""
import os, sys, pyfits

# Parse command line arguments
try:
    runid = str(sys.argv[1])
    itime = int(sys.argv[2])
    newrunid = str(sys.argv[3])
    scale = float(sys.argv[4])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME NEWRUNID SCALE" % execname
    exit

bvars = ["bx", "by", "bz"]
othervars = ["vx", "vy", "vz", "dd", "xn", "xi", "pp"]

for var in bvars + othervars:
    prefix = '%s-%s%4.4i' % (runid, var, itime)
    # read in variable at t=itime
    data = pyfits.open(prefix + '.fits')[0].data
    # scale the B-field by desired factor 
    if var in bvars:
        data *= scale
    # write out the new fits file
    hdu = pyfits.PrimaryHDU()
    hdu.data = data
    hdu.writeto('%s-%s%4.4i.fits' % (newrunid, var, itime), clobber=True)

