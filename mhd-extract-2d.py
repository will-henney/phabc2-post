import pyfits
import os, sys

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try:
    runid = sys.argv[1]
    kslice = int(sys.argv[2])
    extraid = sys.argv[3]
except IndexError, ValueError:
    print "Usage: %s RUNID KSLICE EXTRAID" % execname
    exit

varlist = ["bx", "by", "bz", "vx", "vy", "vz", "dd", "xn", "xi", "pp"]
itime = 0			# use true initial conditions

for var in varlist:
    # read in variable at t=0
    f = pyfits.open('%s-%s%4.4i.fits' % (runid, var, itime))
    # take slice
    f[0].data = f[0].data[kslice,:,:] 
    # write out to new file with extra ID
    f.writeto('%s-%s-%s%4.4i.fits' % (runid, extraid, var, itime))

