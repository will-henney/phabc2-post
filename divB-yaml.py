"""
Calculate Div B for PhabC2 data, writing results for each time to
a separate YAML file.

Uses same discretization that Fabio uses. 
"""
from __future__ import with_statement   # needed for python 2.5
import pyfits
import os, sys
import numpy as N
import yaml

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
try:
    runid = str(sys.argv[1])
    itime1 = int(sys.argv[2])
    itime2 = int(sys.argv[3])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME1 ITIME2" % execname
    exit



for itime in range(itime1, itime2+1):
    try:
        bx, by, bz = [ pyfits.open('%s-%s%4.4i.fits' % (runid, varid, itime))[0].data
                       for varid in 'bx', 'by', 'bz' ]
    except NameError:
        continue                        # if a file is not found, just carry on with the next

    if len(bx.shape) == 2: 
        # promote 2d array to 3d array
        bx = N.array([bx, bx, bx])
        by = N.array([by, by, by])
        bz = N.array([bz, bz, bz])

    print itime, bx.shape
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
    divb = divbx + divby + divbz
    b = N.sqrt(bx**2 + by**2 + bz**2)
    divb_b = divb/b

    # make a dictionary of the statistics we want to save
    data = dict( dbmean=N.abs(divb).mean(), dbmin=divb.min(), dbmax=divb.max(),
                 ndbmean=N.abs(divb_b).mean(), ndbmin=divb_b.min(), ndbmax=divb_b.max() )

    # convert to standard floats since yaml does not recognise most numpy types
    for key, value in data.items(): data[key] = float(value)

    # add in the details of which run and which time
    data.update(runid=runid, itime=itime) 

    with open("divB-%s%4.4i.yaml" % (runid, itime), "w") as f:
        # use verbose more readable format
        f.write(yaml.dump(data, default_flow_style=False))

        

