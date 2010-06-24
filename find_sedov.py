"""
Find at what time the temperature suddenly took off
"""

import sys, glob
import pyfits
from scipy.ndimage import extrema

runid = sys.argv[1]
pattern = "%s-te????.fits" % (runid)
filelist = glob.glob(pattern)

for file in filelist:
    Tarray = pyfits.open(file)[0].data
    Tmin, Tmax, Tminloc, Tmaxloc = extrema(Tarray)
    # switch to fortran array indexing
    Tmaxloc = [ i+1 for i in Tmaxloc[::-1] ]
    print "%s: max T = %.2e at %s" % (file, Tmax, Tmaxloc)

