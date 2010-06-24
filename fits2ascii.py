import pyfits
import os, sys

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try:
    runid = str(sys.argv[1]) 	# name of run
    itime = int(sys.argv[2])	# save time
except IndexError, ValueError:
    print "Usage: %s RUNID" % execname
    exit

varlist = ["bx", "by", "bz", "vx", "vy", "vz", "dd", "xn", "xi", "pp"]

for var in varlist:
    prefix = '%s-%s%4.4i' % (runid, var, itime)
    # read in variable at t=itime
    data = pyfits.open(prefix + '.fits')[0].data
    f = file(prefix + '.dat', 'w') # output file for array
    fmt = data.ndim*"%s " + "\n"
    f.write(fmt % data.shape) # write array dimensions 
    data.tofile(f, sep="\n", format="%7e") # write array
    f.close()

