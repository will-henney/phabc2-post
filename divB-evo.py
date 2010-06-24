"""
Calculate Div B for PhabC2 data, writing results to a FITS file. 

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


columns = ["itime", "dbmean", "dbmin", "dbmax", "ndbmean", "ndbmin", "ndbmax"]
fmts = ["%4.4i"] + 6*["%.3e"]              # formats for the above columns
datatable = []
for itime in range(itime1, itime2+1):
    # get this time's data from the yaml file
    yamlfile = "divB-%s%4.4i.yaml" % (runid, itime)
    try:
        datadict = yaml.load(open(yamlfile))
        # format to string before putting in table
        datatable.append([fmt % (datadict[key]) for key, fmt in zip(columns, fmts)])
    except IOError:
        print "YAML file %s not found - skipping..." % (yamlfile)
        continue
    except KeyError:
        print "Could not read all values from %s - skipping..." % (yamlfile)
        continue

# Now we have made the table, write it to a file
with open("divB-evo-%s.dat" % (runid), "w") as f:
    f.write("# " + "\t".join(columns) + "\n")
    for row in datatable:
        f.write("\t".join(row) + "\n")

        
